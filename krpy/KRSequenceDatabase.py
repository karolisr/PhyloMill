#!/usr/bin/env python
# -*- coding: utf-8 -*-


from __future__ import print_function
from __future__ import unicode_literals


class KRSequenceDatabase:


    ############################################################################


    _DB_SCRIPT = '''

    CREATE TABLE taxonomies(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        taxonomy TEXT NOT NULL UNIQUE
        );

    CREATE TABLE ncbi_tax_ids(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        org_id INTEGER NOT NULL REFERENCES organisms(id) ON DELETE CASCADE,
        ncbi_tax_id INTEGER NOT NULL
        );

    CREATE TABLE organisms(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        active INTEGER NOT NULL,
        taxonomy_id INTEGER REFERENCES taxonomies(id) ON DELETE RESTRICT,
        genus TEXT NOT NULL,
        species TEXT,
        subspecies TEXT,
        variety TEXT,
        hybrid TEXT,
        other TEXT,
        authority TEXT,
        common_name TEXT,
        synonymy_check_done INTEGER NOT NULL
        );

    CREATE TABLE blacklist(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        ncbi_gi INTEGER,
        ncbi_version TEXT,
        internal_reference TEXT,
        notes TEXT
        );

    CREATE TABLE records(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        org_id INTEGER NOT NULL REFERENCES organisms(id) ON DELETE CASCADE,
        active INTEGER NOT NULL,
        ncbi_gi INTEGER,
        ncbi_version TEXT,
        internal_reference TEXT,
        description TEXT
        );

    CREATE TABLE record_ancestry(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        rec_id INTEGER NOT NULL REFERENCES records(id) ON DELETE CASCADE,
        parent_rec_id INTEGER NOT NULL REFERENCES records(id) ON DELETE RESTRICT
        );

    CREATE TABLE record_action_history(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        rec_id INTEGER NOT NULL REFERENCES records(id) ON DELETE CASCADE,
        action TEXT NOT NULL
        );

    CREATE TABLE record_features(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        type TEXT NOT NULL,
        rec_id INTEGER NOT NULL REFERENCES records(id) ON DELETE CASCADE,
        location TEXT NOT NULL
        );

    CREATE TABLE record_feature_qualifiers(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        type TEXT NOT NULL,
        rec_feat_id INTEGER NOT NULL REFERENCES record_features(id) ON DELETE CASCADE,
        qualifier TEXT NOT NULL
        );

    CREATE TABLE record_annotation_types(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        type TEXT NOT NULL UNIQUE
        );

    CREATE TABLE record_annotation_values(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        value TEXT NOT NULL UNIQUE
        );

    CREATE TABLE record_annotations(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        rec_id INTEGER NOT NULL REFERENCES records(id) ON DELETE CASCADE,
        rec_ann_type_id INTEGER NOT NULL REFERENCES record_annotation_types(id) ON DELETE CASCADE,
        rec_ann_value_id INTEGER NOT NULL REFERENCES record_annotation_values(id) ON DELETE CASCADE
        );

    CREATE TABLE sequence_alphabets(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        alphabet TEXT NOT NULL UNIQUE
        );

    INSERT INTO sequence_alphabets VALUES (NULL, 'DNA');
    INSERT INTO sequence_alphabets VALUES (NULL, 'RNA');
    INSERT INTO sequence_alphabets VALUES (NULL, 'AA');

    CREATE TABLE sequences(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        rec_id INTEGER NOT NULL REFERENCES records(id) ON DELETE CASCADE,
        seq_alpha_id INTEGER NOT NULL REFERENCES sequence_alphabets(id) ON DELETE RESTRICT,
        sequence TEXT NOT NULL
        );

    CREATE TABLE sequence_representations(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        seq_id INTEGER NOT NULL REFERENCES sequences(id) ON DELETE CASCADE,
        rec_id INTEGER UNIQUE REFERENCES records(id) ON DELETE SET NULL,
        aln_id INTEGER REFERENCES alignments(id) ON DELETE CASCADE,
        representation SEQREP NOT NULL
        );

    CREATE TABLE alignments(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        name TEXT NOT NULL,
        rec_id INTEGER UNIQUE REFERENCES records(id) ON DELETE CASCADE,
        description TEXT
        );

    '''

    # CREATE TABLE record_features(
    #     id INTEGER PRIMARY KEY AUTOINCREMENT,
    #     rec_id INTEGER NOT NULL REFERENCES records(id) ON DELETE CASCADE,
    #     rec_feat_type_id INTEGER NOT NULL REFERENCES record_feature_types(id) ON DELETE CASCADE,
    #     location TEXT NOT NULL
    #     );

    # CREATE TABLE record_feature_qualifiers(
    #     id INTEGER PRIMARY KEY AUTOINCREMENT,
    #     rec_feat_id INTEGER NOT NULL REFERENCES record_features(id) ON DELETE CASCADE,
    #     rec_feat_qual_type_id INTEGER NOT NULL REFERENCES record_feature_qualifier_types(id) ON DELETE CASCADE,
    #     qualifier TEXT NOT NULL
    #     );

    # CREATE TABLE record_feature_types(
    #     id INTEGER PRIMARY KEY AUTOINCREMENT,
    #     type TEXT NOT NULL UNIQUE
    #     );

    # CREATE TABLE record_feature_qualifier_types(
    #     id INTEGER PRIMARY KEY AUTOINCREMENT,
    #     type TEXT NOT NULL UNIQUE
    #     );

    ############################################################################


    class _KRSeqEdits(list):
        pass


    def produce_seq_edits(self, s1, s2):
        from krpy import krstring
        return self._KRSeqEdits(krstring.produce_edits(s1, s2))


    def apply_seq_edits(self, e, s):
        from krpy import krstring
        return krstring.apply_edits(e, s)


    def _adapt_seq_edits(self, e):
        from krpy import krstring
        return krstring.edits_to_string(e)


    def _convert_seq_edits(self, e):
        from krpy import krstring
        return self._KRSeqEdits(krstring.string_to_edits(e))


    ############################################################################


    def __init__(self, db_file):

        import os
        import sqlite3

        self._DB_FILE = db_file

        if not os.path.exists(self._DB_FILE):
            self._db_prepare()

        sqlite3.register_adapter(self._KRSeqEdits, self._adapt_seq_edits)
        sqlite3.register_converter(b'SEQREP', self._convert_seq_edits)

        self._DB_CONN = sqlite3.connect(self._DB_FILE,
            detect_types=sqlite3.PARSE_DECLTYPES, isolation_level=None)
        self._DB_CONN.row_factory = sqlite3.Row
        self._DB_CONN.text_factory = str

        self._DB_CURSOR = self._DB_CONN.cursor()

        self._DB_CONN.execute('PRAGMA foreign_keys = ON;')
        # self._DB_CONN.execute('PRAGMA cache_size = 100000;')
        # self._DB_CONN.execute('PRAGMA max_page_count = 100000;')
        self._DB_CONN.execute('PRAGMA synchronous = OFF;')
        self._DB_CONN.execute('PRAGMA journal_mode = MEMORY;')


    ############################################################################


    def _db_prepare(self):

        import sys
        import sqlite3

        db = sqlite3.connect(self._DB_FILE)

        try:
            db.executescript(self._DB_SCRIPT)
        except sqlite3.Error as error:
            print(error, file=sys.stderr)
            db.rollback()
            raise

        db.commit()
        db.close()


    def db_select(self, table_name_list, column_list, where_dict=None, join_rules_str=None, order_by_column_list=None):

        import sys
        import sqlite3

        columns_str = ', '.join(column_list)

        if where_dict:
            where_str = ' WHERE ' + ' AND '.join(
                [str(x) + ' IS :' + str(x.replace('.', '_')) for x in where_dict.keys()])

            # This solves the issue where sqlite3 interpolator freaks out if
            # there are periods in value placeholders
            new_where_dict = dict()
            for key in where_dict.keys():
                new_where_dict[key.replace('.', '_')] = where_dict[key]
            where_dict = new_where_dict

        else:
            where_str = ''
            where_dict = {}

        order_by_str = ''
        if order_by_column_list:
            order_by_str = ' ORDER BY ' + ', '.join(order_by_column_list)

        tables_str = ''

        if len(table_name_list) == 1:
            tables_str = table_name_list[0]
        else:
            tables_str = ' LEFT OUTER JOIN '.join(table_name_list)
            tables_str = tables_str + ' ON ' + join_rules_str

        select_str = str('SELECT ' + columns_str + ' FROM ' + tables_str + where_str + order_by_str)

        # print(select_str)

        try:
            self._DB_CURSOR.execute(select_str, where_dict)
        except sqlite3.Error as error:
            print(error, file=sys.stderr)
            raise

        results = self._DB_CURSOR.fetchall()

        return results


    def db_get_row_ids(self, table_name, where_dict=None):

        results = self.db_select(
            table_name_list=[table_name],
            column_list=['id'],
            where_dict=where_dict)

        # row_id = None

        # if len(results) == 1:
        #     row_id = results[0][b'id']
        # elif len(results) > 1:
        #     row_id = -1

        row_ids = None

        if len(results) > 0:
            row_ids = [x[b'id'] for x in results]

        return row_ids


    def db_insert(self, table_name, values_dict, check_exists=True):

        import sys
        import sqlite3

        row_id = None
        already_in_db = False

        columns_str = ', '.join(values_dict.keys())
        values_str = ':' + ', :'.join(values_dict.keys())

        exists = False
        if check_exists:
            exists = self.db_get_row_ids(table_name, values_dict)

        if check_exists and exists:
            row_id = exists[0]
            already_in_db = True

        else:
            insert_str = str('INSERT INTO ' + table_name + ' (' + columns_str +
                ') ' + 'VALUES (' + values_str + ');')

            try:
                self._DB_CURSOR.execute(insert_str, values_dict)
            except sqlite3.Error as error:
                print(error, file=sys.stderr)
                self._DB_CONN.rollback()
                raise

            row_id = self._DB_CURSOR.lastrowid

        return (row_id, already_in_db)


    def db_update(self, table_name, values_dict, where_dict=None):

        import sys
        import sqlite3

        values_str_list = [str(x) + '=:' + str(x) for x in values_dict.keys()]
        values_str = ', '.join(values_str_list)

        where_str = ''
        if where_dict:
            where_str = ' WHERE ' + ' AND '.join(
                [str(x) + ' IS :' + str(x) + '_where_string' for x in where_dict.keys()])

            # This solves the issue where sqlite3 interpolator freaks out if
            # there are periods in value placeholders
            new_where_dict = dict()
            for key in where_dict.keys():
                new_where_dict[str(key) + str('_where_string')] = where_dict[key]
            where_dict = new_where_dict

        else:
            where_str = ''
            where_dict = {}

        update_str = str(
            'UPDATE ' + table_name + ' SET ' + values_str + where_str + ';')

        update_dict = dict(values_dict.items() + where_dict.items())

        try:
            self._DB_CURSOR.execute(update_str, update_dict)
        except sqlite3.Error as error:
            print(error, file=sys.stderr)
            self._DB_CONN.rollback()
            raise

        return self._DB_CURSOR.lastrowid


    def db_delete(self, table_name, where_dict=None):

        import sys
        import sqlite3

        where_str = ''
        if where_dict:
            where_str = ' WHERE ' + ' AND '.join(
                [str(x) + ' IS :' + str(x) for x in where_dict.keys()])
        else:
            where_str = ''
            where_dict = {}

        delete_string = str(
            'DELETE FROM  ' + table_name + where_str + ';')

        try:
            self._DB_CURSOR.execute(delete_string, where_dict)
        except sqlite3.Error as error:
            print(error, file=sys.stderr)
            self._DB_CONN.rollback()
            raise


    ############################################################################


    def get_cursor(self):
        return self._DB_CURSOR


    def save(self):
        self._DB_CONN.commit()


    def close(self):
        self._DB_CONN.close()


    ############################################################################


    def bio_alphabet_to_string(self, bio_alphabet):

        import Bio

        sequence_alphabet_str = ''

        if (
            isinstance(bio_alphabet, Bio.Alphabet.IUPAC.IUPACAmbiguousDNA) or
            isinstance(bio_alphabet, Bio.Alphabet.IUPAC.IUPACUnambiguousDNA) or
            isinstance(bio_alphabet, Bio.Alphabet.IUPAC.ExtendedIUPACDNA)
            ):
            sequence_alphabet_str = 'DNA'
        elif (
            isinstance(bio_alphabet, Bio.Alphabet.IUPAC.IUPACAmbiguousRNA) or
            isinstance(bio_alphabet, Bio.Alphabet.IUPAC.IUPACUnambiguousRNA)
            ):
            sequence_alphabet_str = 'RNA'
        elif (
            isinstance(bio_alphabet, Bio.Alphabet.IUPAC.ExtendedIUPACProtein) or
            isinstance(bio_alphabet, Bio.Alphabet.IUPAC.IUPACProtein)
            ):
            sequence_alphabet_str = 'AA'

        return sequence_alphabet_str


    def string_to_bio_alphabet(self, str_alphabet):

        import Bio

        bio_alphabet = None

        if str_alphabet == 'DNA':
            bio_alphabet = Bio.Alphabet.IUPAC.IUPACAmbiguousDNA()
        elif str_alphabet == 'RNA':
            bio_alphabet = Bio.Alphabet.IUPAC.IUPACAmbiguousRNA()
        elif str_alphabet == 'AA':
            bio_alphabet = Bio.Alphabet.IUPAC.IUPACProtein()

        return bio_alphabet


    def _sequence_alphabet(self, bio_seq_object):

        alphabet = bio_seq_object.alphabet
        return self.bio_alphabet_to_string(bio_alphabet=alphabet)


    ############################################################################


    def set_active(self, table_name, where_dict):

        self.db_update(
            table_name=table_name,
            values_dict={'active': 1},
            where_dict=where_dict)


    def set_inactive(self, table_name, where_dict):

        self.db_update(
            table_name=table_name,
            values_dict={'active': 0},
            where_dict=where_dict)


    def add_organism(self, organism_dict, taxonomy_list=None,
        ncbi_tax_id_list=None, synonymy_check_done=False, active=True):

        taxonomy_id = None
        taxonomy_str = None

        if taxonomy_list:
            taxonomy_str = ','.join(taxonomy_list)
            values_dict = {'taxonomy': taxonomy_str}
            taxonomy_id = self.db_insert('taxonomies', values_dict)[0]

        for key in organism_dict.keys():
            if str(organism_dict[key]) == b'':
                organism_dict[key] = None

        if synonymy_check_done:
            synonymy_check_done = 1
        else:
            synonymy_check_done = 0

        if active:
            active = 1
        else:
            active = 0

        where_dict = {
            # 'taxonomy_id': taxonomy_id,
            'genus': organism_dict['genus'],
            'species': organism_dict['species'],
            'subspecies': organism_dict['subspecies'],
            'variety': organism_dict['variety'],
            'hybrid': organism_dict['hybrid'],
            'other': organism_dict['other']
            # 'authority': organism_dict['authority'],
            # 'synonymy_check_done': synonymy_check_done
        }

        values_dict = {
            'active': active,
            'taxonomy_id': taxonomy_id,
            'genus': organism_dict['genus'],
            'species': organism_dict['species'],
            'subspecies': organism_dict['subspecies'],
            'variety': organism_dict['variety'],
            'hybrid': organism_dict['hybrid'],
            'other': organism_dict['other'],
            'authority': organism_dict['authority'],
            'common_name': organism_dict['common_name'],
            'synonymy_check_done': synonymy_check_done
        }

        already_in_db = False

        org_id = None
        org_ids = self.db_get_row_ids('organisms', where_dict=where_dict)

        if org_ids and len(org_ids) == 1:

            org_id = org_ids[0]

            org_dict = self.db_select(
                table_name_list=['organisms'],
                column_list=['taxonomy_id', 'authority', 'synonymy_check_done'],
                where_dict={'id': org_id},
                join_rules_str=None,
                order_by_column_list=None)[0]

            if not values_dict[b'taxonomy_id']:
                values_dict[b'taxonomy_id'] = org_dict[b'taxonomy_id']
            if not values_dict[b'authority']:
                values_dict[b'authority'] = org_dict[b'authority']
            if not values_dict[b'synonymy_check_done']:
                values_dict[b'synonymy_check_done'] = org_dict[b'synonymy_check_done']

            self.db_update('organisms',
                values_dict=values_dict,
                where_dict={'id': org_id})
            already_in_db = True

        else:
            org_id = self.db_insert('organisms', values_dict)[0]

        if ncbi_tax_id_list:
            for ncbi_tax_id in ncbi_tax_id_list:

                # CREATE TABLE ncbi_tax_ids(
                #     id INTEGER PRIMARY KEY AUTOINCREMENT,
                #     org_id INTEGER NOT NULL REFERENCES organisms(id) ON DELETE CASCADE,
                #     ncbi_tax_id INTEGER NOT NULL

                values_dict = {
                    'org_id': org_id,
                    'ncbi_tax_id': ncbi_tax_id
                    }

                self.db_insert('ncbi_tax_ids', values_dict)

        return (org_id, already_in_db)


    def delete_organisms(self, where_dict):

        table_name = 'organisms'

        self.db_delete(
            table_name=table_name,
            where_dict=where_dict)


    def delete_orphaned_organisms(self):

        delete_string = 'DELETE FROM organisms WHERE id NOT IN \
            (SELECT org_id FROM records WHERE org_id IS NOT NULL)'

        try:
            self._DB_CURSOR.execute(delete_string)
        except sqlite3.Error as error:
            print(error, file=sys.stderr)
            self._DB_CONN.rollback()
            raise


    def delete_orphaned_taxonomies(self):

        delete_string = 'DELETE FROM taxonomies WHERE id NOT IN \
            (SELECT taxonomy_id FROM organisms WHERE taxonomy_id IS NOT NULL)'

        try:
            self._DB_CURSOR.execute(delete_string)
        except sqlite3.Error as error:
            print(error, file=sys.stderr)
            self._DB_CONN.rollback()
            raise


    def get_organisms(self, where_dict=None):

        table_name_list = ['organisms', 'taxonomies']
        column_list = ['organisms.id', 'genus', 'species',
                       'subspecies', 'variety', 'hybrid', 'other', 'authority',
                       'taxonomy', 'common_name', 'synonymy_check_done']
        order_by_column_list = ['genus', 'species']

        orgs = self.db_select(
            table_name_list=table_name_list,
            column_list=column_list,
            where_dict=where_dict,
            join_rules_str='organisms.taxonomy_id=taxonomies.id',
            order_by_column_list=order_by_column_list)

        results = list()

        for org in orgs:

            org = dict(org)

            ncbi_tax_ids = self.db_select(
                table_name_list=['ncbi_tax_ids'],
                column_list=['ncbi_tax_id'],
                where_dict={'org_id': org[b'id']},
                join_rules_str=None,
                order_by_column_list=None)

            # print(ncbi_tax_ids)

            # for x in ncbi_tax_ids:
            #     print(x[b'ncbi_tax_id'])

            # print([x[b'ncbi_tax_id'] for x in ncbi_tax_ids])

            org[b'ncbi_tax_ids'] = [x[b'ncbi_tax_id'] for x in ncbi_tax_ids]

            results.append(org)

        return results


    def add_sequence(self, rec_id, sequence_str, sequence_alphabet_str):

        sequence_str = sequence_str.upper()
        sequence_alphabet_str = sequence_alphabet_str.upper()

        values_dict = {'alphabet': sequence_alphabet_str}
        seq_alpha_id = self.db_get_row_ids('sequence_alphabets',
            values_dict)[0]

        values_dict = {
            'rec_id': rec_id,
            'seq_alpha_id': seq_alpha_id,
            'sequence': sequence_str
        }

        row_id = self.db_insert('sequences', values_dict)

        return row_id


    def add_sequence_representation(
        self, seq_id, repr_list, rec_id=None, aln_id=None):

        values_dict = {
            'aln_id': aln_id,
            'rec_id': rec_id,
            'seq_id': seq_id,
            'representation': repr_list
        }

        row_id = self.db_insert('sequence_representations', values_dict,
            check_exists=False)

        return row_id


    def add_record_feature(self, rec_id, type_str, location_str):

        # values_dict = {'type': type_str}
        # rec_feat_type_id = self.db_insert('record_feature_types',
        #     values_dict)[0]

        values_dict = {
            'type': type_str,
            'rec_id': rec_id,
            # 'rec_feat_type_id': rec_feat_type_id,
            'location': location_str
        }

        row_id = self.db_insert('record_features', values_dict, check_exists=False)

        return row_id


    def add_record_feature_qualifier(self, rec_feat_id, type_str, qualifier_str):

        # values_dict = {'type': type_str}
        # rec_feat_qual_type_id = self.db_insert(
        #     'record_feature_qualifier_types',
        #     values_dict)[0]

        values_dict = {
            'type': type_str,
            'rec_feat_id': rec_feat_id,
            # 'rec_feat_qual_type_id': rec_feat_qual_type_id,
            'qualifier': qualifier_str
        }

        row_id = self.db_insert('record_feature_qualifiers', values_dict, check_exists=False)

        return row_id


    def _add_record_annotation(self, rec_id, type_str, annotation_str):

        values_dict = {'type': type_str}
        rec_ann_type_id = self.db_insert('record_annotation_types',
            values_dict)[0]

        values_dict = {'value': annotation_str}
        rec_ann_value_id = self.db_insert('record_annotation_values',
            values_dict)[0]

        values_dict = {
            'rec_id': rec_id,
            'rec_ann_type_id': rec_ann_type_id,
            'rec_ann_value_id': rec_ann_value_id
        }

        row_id = self.db_insert('record_annotations', values_dict)

        return row_id


    def add_record_annotation(
        self,
        record_reference,
        type_str,
        annotation_str,
        record_reference_type='gi'  # gi version internal raw
        ):

        where_dict_key = ''

        if record_reference_type == 'gi':
            where_dict_key = 'ncbi_gi'
        elif record_reference_type == 'version':
            where_dict_key = 'ncbi_version'
        elif record_reference_type == 'internal':
            where_dict_key = 'internal_reference'
        elif record_reference_type == 'raw':
            where_dict_key = 'id'

        where_dict = {where_dict_key: record_reference}

        rec_id = self.db_get_row_ids(
            table_name='records',
            where_dict=where_dict)[0]

        row_id = self._add_record_annotation(
            rec_id=rec_id,
            type_str=type_str,
            annotation_str=annotation_str)

        return row_id



    def _add_record_action(self, rec_id, action_str):

        values_dict = {
            'rec_id': rec_id,
            'action': action_str
        }

        row_id = self.db_insert('record_action_history', values_dict)

        return row_id


    def in_db(self,
        record_reference,
        record_reference_type='gi'  # gi version internal
        ):

        where_dict_key = ''

        if record_reference_type == 'gi':
            where_dict_key = 'ncbi_gi'
        elif record_reference_type == 'version':
            where_dict_key = 'ncbi_version'
        elif record_reference_type == 'internal':
            where_dict_key = 'internal_reference'

        where_dict = {where_dict_key: record_reference}

        results = self.db_get_row_ids(
            table_name='records',
            where_dict=where_dict)

        if not results:
            return False
        else:
            return True


    def get_all_record_ids(self,
        record_reference_type='gi'  # gi version internal
        ):

        where_dict_key = ''

        if record_reference_type == 'gi':
            where_dict_key = 'ncbi_gi'
        elif record_reference_type == 'version':
            where_dict_key = 'ncbi_version'
        elif record_reference_type == 'internal':
            where_dict_key = 'internal_reference'

        # where_dict = {where_dict_key: record_reference}

        results = self.db_select(
            table_name_list=['records'],
            column_list=[where_dict_key],
            where_dict=None,
            join_rules_str=None,
            order_by_column_list=None)

        results = [x[0] for x in results]

        return results


    def in_blacklist(self,
        record_reference,
        record_reference_type='gi'  # gi version internal
        ):

        where_dict_key = ''

        if record_reference_type == 'gi':
            where_dict_key = 'ncbi_gi'
        elif record_reference_type == 'version':
            where_dict_key = 'ncbi_version'
        elif record_reference_type == 'internal':
            where_dict_key = 'internal_reference'

        where_dict = {where_dict_key: record_reference}

        results = self.db_get_row_ids(
            table_name='blacklist',
            where_dict=where_dict)

        if not results:
            return False
        else:
            return True


    def add_record_to_blacklist(self, ncbi_gi=None, ncbi_version=None,
        internal_reference=None, notes=None):

        row_id = None

        if ncbi_gi or ncbi_version or internal_reference:

            values_dict = {
                'ncbi_gi': ncbi_gi,
                'ncbi_version': ncbi_version,
                'internal_reference': internal_reference,
                'notes': notes
            }

            row_id = self.db_insert('blacklist', values_dict)

        return row_id


    def delete_records(self, where_dict, blacklist=False, blacklist_notes=None):

        table_name = 'records'
        column_list = ['ncbi_gi', 'ncbi_version', 'internal_reference']

        rec_ids_list = self.db_select(
            table_name_list=[table_name],
            column_list=column_list,
            where_dict=where_dict,
            join_rules_str=None)

        self.db_delete(
            table_name=table_name,
            where_dict=where_dict)

        if blacklist:
            for rec_ids in rec_ids_list:
                self.add_record_to_blacklist(
                    ncbi_gi=rec_ids[b'ncbi_gi'],
                    ncbi_version=rec_ids[b'ncbi_version'],
                    internal_reference=rec_ids[b'internal_reference'],
                    notes=blacklist_notes)


    def delete_record(self,
        record_reference,
        record_reference_type='gi',  # gi version internal
        blacklist=False,
        blacklist_notes=None
        ):

        where_dict_key = ''

        if record_reference_type == 'gi':
            where_dict_key = 'ncbi_gi'
        elif record_reference_type == 'version':
            where_dict_key = 'ncbi_version'
        elif record_reference_type == 'internal':
            where_dict_key = 'internal_reference'

        where_dict = {where_dict_key: record_reference}

        self.delete_records(where_dict=where_dict, blacklist=blacklist, blacklist_notes=None)


    def update_records(self, values_dict, where_dict):

        self.db_update(
            table_name='records',
            values_dict=values_dict,
            where_dict=where_dict)


    # def add_record(self,
    #     org,  # string or org_dict
    #     ncbi_gi,
    #     ncbi_version,
    #     internal_reference,
    #     description,
    #     sequence_str,
    #     # seq_rep_id=None,
    #     # aln_id=None,
    #     taxonomy_list=None,
    #     ncbi_tax_id=None,
    #     action_str='New record'):

    #     from krpy import krbionames

    #     if isinstance(org, basestring):
    #         org = krbionames.parse_organism_name(org, sep=' ', ncbi_authority=False)

    #     pass


    def add_record(self, org_id, ncbi_gi, ncbi_version, internal_reference,
        description, sequence_str, sequence_alphabet_str,
        parent_rec_id_list=None, active=True, action_str='New record'):

        # id INTEGER PRIMARY KEY AUTOINCREMENT,
        # org_id INTEGER NOT NULL REFERENCES organisms(id) ON DELETE CASCADE,
        # active INTEGER NOT NULL,
        # ncbi_gi INTEGER,
        # ncbi_version TEXT,
        # internal_reference TEXT,
        # description TEXT


        # Record

        if active:
            active = 1
        else:
            active = 0

        values_dict = {
            'org_id': org_id,
            'active': active,
            'ncbi_gi': ncbi_gi,
            'ncbi_version': ncbi_version,
            'internal_reference': internal_reference,
            'description': description
        }

        rec_id = self.db_insert('records', values_dict)[0]

        self._add_record_action(
            rec_id=rec_id,
            action_str=action_str)

        # Sequence

        seq_id = self.add_sequence(
            rec_id=rec_id,
            sequence_str=sequence_str,
            sequence_alphabet_str=sequence_alphabet_str)[0]

        repr_list = self.produce_seq_edits(sequence_str, sequence_str)

        seq_rep_id = self.add_sequence_representation(
            rec_id=rec_id,
            seq_id=seq_id,
            repr_list=repr_list)[0]

        # Ancestral records

        if parent_rec_id_list:

            # id INTEGER PRIMARY KEY AUTOINCREMENT,
            # rec_id INTEGER NOT NULL REFERENCES records(id) ON DELETE CASCADE,
            # parent_rec_id INTEGER NOT NULL REFERENCES records(id) ON DELETE RESTRICT

            for parent_rec_id in parent_rec_id_list:

                values_dict = {
                    'rec_id': rec_id,
                    'parent_rec_id': parent_rec_id
                    }

                self.db_insert('record_ancestry', values_dict)

        return rec_id


    def add_genbank_record(self, record, action_str='Genbank record'):

        from krpy import krbionames
        from krpy import krncbi

        # Organism

        org = record.annotations['organism']
        org_dict = krbionames.parse_organism_name(org, sep=' ',
            ncbi_authority=False)
        taxonomy_list_temp = record.annotations['taxonomy']
        taxonomy_list = [('name=' + x) for x in taxonomy_list_temp]

        ncbi_tax_id = int(krncbi.get_ncbi_tax_id_for_record(record))

        org_id = self.add_organism(
            organism_dict=org_dict,
            taxonomy_list=taxonomy_list,
            ncbi_tax_id_list=[ncbi_tax_id])[0]

        # Record

        active = 1

        ncbi_gi = record.annotations['gi']

        ncbi_version = record.id

        description = record.description

        values_dict = {
            'org_id': org_id,
            'active': active,
            'ncbi_gi': ncbi_gi,
            'ncbi_version': ncbi_version,
            'description': description
        }

        row_id = self.db_insert('records', values_dict)

        self._add_record_action(
            rec_id=row_id[0],
            action_str=action_str)

        # Sequence

        sequence_str = str(record.seq)
        sequence_alphabet_str = self._sequence_alphabet(record.seq)

        seq_id = self.add_sequence(
            rec_id=row_id[0],
            sequence_str=sequence_str,
            sequence_alphabet_str=sequence_alphabet_str)[0]

        repr_list = self.produce_seq_edits(sequence_str, sequence_str)

        seq_rep_id = self.add_sequence_representation(
            rec_id=row_id[0],
            seq_id=seq_id,
            repr_list=repr_list)[0]

        # Features

        for feature in record.features:

            rec_feat_id = self.add_record_feature(
                rec_id=row_id[0],
                type_str=feature.type,
                location_str=str(feature.location))

            # Feature qualifiers

            for qualifier in feature.qualifiers:

                fqual = None
                if isinstance(feature.qualifiers[qualifier], list):
                    fqual = feature.qualifiers[qualifier][0]
                else:
                    fqual = feature.qualifiers[qualifier]

                rec_feat_qual_id = self.add_record_feature_qualifier(
                    rec_feat_id=rec_feat_id[0],
                    type_str=qualifier,
                    qualifier_str=str(fqual))

        return row_id


    # def _align_sequence_reps(self, seq_rep_id_list, program, options='', program_executable=''):

    #     from Bio import SeqRecord
    #     from krpy import kralign

    #     temp_records = list()

    #     for seq_rep_id in seq_rep_id_list:
    #         seq = self._get_sequence_from_representation(seq_rep_id=seq_rep_id)
    #         seq_record = SeqRecord.SeqRecord(seq, id=str(seq_rep_id), name='', description='')
    #         temp_records.append(seq_record)

    #     aln = kralign.align(
    #         records=temp_records,
    #         program=program,
    #         options='',
    #         program_executable='')

    #     result_seq_rep_id_list = list()

    #     for s in aln:
    #         old_seq_rep_id = int(s.id)

    #         table_name_list = ['sequence_representations']
    #         column_list = ['seq_id']
    #         where_dict = {'id': old_seq_rep_id}
    #         join_rules_str = None

    #         results = self.db_select(
    #             table_name_list=table_name_list,
    #             column_list=column_list,
    #             where_dict=where_dict,
    #             join_rules_str=join_rules_str)[0]

    #         seq_id = results[b'seq_id']
    #         seq = self._get_sequence(sequence_id=seq_id)
    #         aln_seq = str(s.seq).upper()

    #         repr_list = self.produce_seq_edits(str(seq), aln_seq)

    #         new_seq_rep_id = self.add_sequence_representation(
    #             seq_id=seq_id, repr_list=repr_list)[0]

    #         result_seq_rep_id_list.append(new_seq_rep_id)

    #     return result_seq_rep_id_list


    def add_alignment(self, name, seq_rep_id_list, description=None,rec_id=None):

        values_dict = {
            'rec_id': rec_id,
            'name': name,
            'description': description
        }

        row_id = self.db_insert(
            'alignments', values_dict, check_exists=False)[0]

        for seq_rep_id in seq_rep_id_list:

            values_dict = {
                'aln_id': row_id
            }

            where_dict = {
                'id': seq_rep_id
            }

            self.db_update(
                table_name='sequence_representations',
                values_dict=values_dict,
                where_dict=where_dict)


    def get_alignment(self, alignment_id):

        # CREATE TABLE sequence_representations(
        #     id INTEGER PRIMARY KEY AUTOINCREMENT,
        #     seq_id INTEGER NOT NULL REFERENCES sequences(id) ON DELETE CASCADE,
        #     rec_id INTEGER UNIQUE REFERENCES records(id) ON DELETE SET NULL,
        #     aln_id INTEGER REFERENCES alignments(id) ON DELETE CASCADE,
        #     representation SEQREP NOT NULL
        #     );

        # CREATE TABLE alignments(
        #     id INTEGER PRIMARY KEY AUTOINCREMENT,
        #     name TEXT NOT NULL,
        #     rec_id INTEGER UNIQUE REFERENCES records(id) ON DELETE CASCADE,
        #     description TEXT
        #     );

        from Bio.Align import MultipleSeqAlignment
        from Bio import Seq
        from Bio.SeqRecord import SeqRecord

        aln = None

        # table_name_list = ['alignments']
        # column_list = ['name', 'rec_id', 'description']
        # where_dict = {'id': alignment_id}
        # alignment_field_list = self.db_select(
        #     table_name_list=table_name_list,
        #     column_list=column_list,
        #     where_dict=where_dict,
        #     join_rules_str=None)[0]

        table_name_list = ['sequence_representations']
        column_list = ['id', 'seq_id', 'rec_id', 'aln_id', 'representation']
        where_dict = {'aln_id': alignment_id}
        seq_reps = self.db_select(
            table_name_list=table_name_list,
            column_list=column_list,
            where_dict=where_dict,
            join_rules_str=None)

        seq_record_list = list()

        for seq_rep in seq_reps:

            rec_id = self.db_select(
                table_name_list=['sequences'],
                column_list=['rec_id'],
                where_dict={'id': seq_rep[b'seq_id']},
                join_rules_str=None,
                order_by_column_list=None)[0][b'rec_id']

            seq_rep_id = seq_rep[b'id']
            seq = self._get_sequence_from_representation(seq_rep_id=seq_rep_id)

            seq_record = SeqRecord(
                seq=seq,
                id=str(rec_id),
                name='',
                description='')

            seq_record_list.append(seq_record)

            # print(seq_record)

        aln = MultipleSeqAlignment(seq_record_list)

        return aln


    # def align_records(
    #     self,
    #     record_reference_list,
    #     program,
    #     aln_name,
    #     options='',
    #     program_executable='',
    #     record_reference_type='gi',  # gi version internal
    #     description=None
    #     ):

    #     seq_rep_id_list = list()

    #     for r in record_reference_list:
    #         seq_rep_id = self.get_seq_rep_id_for_record(
    #             record_reference=r,
    #             record_reference_type=record_reference_type)
    #         seq_rep_id_list.append(seq_rep_id)

    #     aln_seq_rep_id_list = self._align_sequence_reps(
    #         seq_rep_id_list=seq_rep_id_list,
    #         program=program,
    #         options=options,
    #         program_executable=program_executable)

    #     row_id = self.add_alignment(
    #         name=aln_name,
    #         seq_rep_id_list=aln_seq_rep_id_list,
    #         description=description)

    #     return row_id


    ############################################################################


    def _get_sequence(self, sequence_id):

        from Bio import Seq

        table_name_list = ['sequences', 'sequence_alphabets']
        column_list = ['alphabet', 'sequence']
        where_dict = {'sequences.id': sequence_id}
        join_rules_str = 'sequences.seq_alpha_id=sequence_alphabets.id'

        results = self.db_select(
            table_name_list=table_name_list,
            column_list=column_list,
            where_dict=where_dict,
            join_rules_str=join_rules_str)[0]

        seq_str = results[b'sequence']
        seq_alphabet = self.string_to_bio_alphabet(str_alphabet=results[b'alphabet'])

        seq = Seq.Seq(data=seq_str, alphabet=seq_alphabet)

        return(seq)


    def _get_sequence_from_representation(self, seq_rep_id):

        from Bio import Seq

        table_name_list = ['sequence_representations']
        column_list = ['seq_id', 'representation']
        where_dict = {'id': seq_rep_id}
        join_rules_str = None

        results = self.db_select(
            table_name_list=table_name_list,
            column_list=column_list,
            where_dict=where_dict,
            join_rules_str=join_rules_str)[0]

        sequence_id = results[b'seq_id']
        seq = self._get_sequence(sequence_id)

        sequence_edits = results[b'representation']

        new_seq_str = self.apply_seq_edits(e=sequence_edits, s=str(seq))
        new_seq = Seq.Seq(data=new_seq_str, alphabet=seq.alphabet)

        return new_seq


    def get_seq_rep_id_for_record(self,
        record_reference,
        record_reference_type='gi'  # gi version internal raw
        ):

        where_dict_key = ''

        if record_reference_type == 'gi':
            where_dict_key = 'ncbi_gi'
        elif record_reference_type == 'version':
            where_dict_key = 'ncbi_version'
        elif record_reference_type == 'internal':
            where_dict_key = 'internal_reference'
        elif record_reference_type == 'raw':
            where_dict_key = 'id'

        table_name = 'records'
        where_dict = {where_dict_key: record_reference}
        rec_id = self.db_get_row_ids(
            table_name=table_name,
            where_dict=where_dict)[0]

        table_name = 'sequence_representations'
        where_dict = {'rec_id': rec_id}
        seq_rep_id = self.db_get_row_ids(
            table_name=table_name,
            where_dict=where_dict)[0]

        return seq_rep_id


    def get_sequence_for_record(self,
        record_reference,
        record_reference_type='gi'  # gi version internal raw
        ):

        seq_rep_id = self.get_seq_rep_id_for_record(
            record_reference=record_reference,
            record_reference_type=record_reference_type)

        seq = self._get_sequence_from_representation(seq_rep_id=seq_rep_id)

        return seq


    def get_sequences_for_records(self,
        record_reference_list,
        record_reference_type='gi'  # gi version internal raw
        ):

        seq_list = list()

        for r in record_reference_list:
            seq = self.get_sequence_for_record(
                record_reference=r,
                record_reference_type=record_reference_type)

            seq_list.append({'reference': r, 'seq': seq})

        return seq_list


    def get_record(
        self,
        record_reference,
        record_reference_type='gi'  # gi version internal raw
        ):

        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature
        from krpy import krbionames
        from krpy import krseq

        where_dict_key = ''

        if record_reference_type == 'gi':
            where_dict_key = 'ncbi_gi'
        elif record_reference_type == 'version':
            where_dict_key = 'ncbi_version'
        elif record_reference_type == 'internal':
            where_dict_key = 'internal_reference'
        elif record_reference_type == 'raw':
            where_dict_key = 'id'

        where_dict = {where_dict_key: record_reference}
        results = self.db_select(
            table_name_list=['records'],
            column_list=['active', 'id', 'org_id', 'ncbi_gi', 'ncbi_version', 'internal_reference',
                         'description'],
            where_dict=where_dict,
            join_rules_str=None)[0]

        where_dict = {b'organisms.id': results[b'org_id']}

        org_dict = self.get_organisms(where_dict=where_dict)[0]
        org_flat = krbionames.flatten_organism_name(
                parsed_name=org_dict, sep=' ')

        seq_rep_id = self.get_seq_rep_id_for_record(
            record_reference=record_reference,
            record_reference_type=record_reference_type)

        seq = self._get_sequence_from_representation(seq_rep_id=seq_rep_id)

        features_temp = self.db_select(
            # table_name_list=['record_features', 'record_feature_types'],
            table_name_list=['record_features'],
            column_list=['record_features.id', 'location', 'type'],
            where_dict={'rec_id': results[b'id']},
            # join_rules_str='record_features.rec_feat_type_id=record_feature_types.id',
            order_by_column_list=None)

        features = list()

        for feat_raw in features_temp:

            location_string = feat_raw[b'location']
            location = krseq.location_from_string(
                location_string=location_string)

            qualifiers_temp = self.db_select(
                # table_name_list=['record_feature_qualifiers', 'record_feature_qualifier_types'],
                table_name_list=['record_feature_qualifiers'],
                column_list=['qualifier', 'type'],
                where_dict={'rec_feat_id': feat_raw[b'id']},
                # join_rules_str='record_feature_qualifiers.rec_feat_qual_type_id=record_feature_qualifier_types.id',
                order_by_column_list=None)

            qualifiers = dict()
            for qualifier in qualifiers_temp:
                qualifiers[qualifier[b'type']] = [qualifier[b'qualifier']]

            feat = SeqFeature(
                location=location,
                type=feat_raw[b'type'],
                qualifiers=qualifiers
                )

            features.append(feat)

        rec_id = results[b'ncbi_version']
        if not rec_id:
            rec_id = ''

        rec_name = results[b'ncbi_version']
        if not rec_name:
            rec_name = ''

        rec_description = results[b'description']
        if not rec_description:
            rec_description = ''

        record = SeqRecord(
            seq=seq,
            id=rec_id,
            name=rec_name,
            description=rec_description,
            features=features)

        record.annotations[b'gi'] = str(results[b'ncbi_gi'])
        record.annotations[b'organism'] = str(org_flat)
        record.annotations[b'common_name'] = str(org_dict['common_name'])
        record.annotations[b'lineage'] = str(org_dict['taxonomy'])
        record.annotations[b'internal_reference'] = str(results[b'internal_reference'])
        record.annotations[b'kr_seq_db_org_id'] = str(results[b'org_id'])
        record.annotations[b'kr_seq_db_id'] = str(results[b'id'])
        record.annotations[b'kr_seq_db_active'] = str(results[b'active'])

        return record


    def get_records(
        self,
        record_reference_list,
        record_reference_type='gi',  # gi version internal raw
        active=True,
        inactive=False
        ):

        if (not active) and (not inactive):
            return None

        record_list = list()

        for ref in record_reference_list:
            record = self.get_record(
                record_reference=ref,
                record_reference_type=record_reference_type  # gi version internal raw
                )

            active_state = int(record.annotations[b'kr_seq_db_active'])

            accept = False

            if active and inactive:
                accept = True
            elif active and active_state==1:
                accept = True
            elif inactive and active_state==0:
                accept = True

            if accept:
                record_list.append(record)

        return record_list


    def get_all_records(self, active=True, inactive=False):

        # record_reference_list = self.db_select(
        #     table_name_list=['records'],
        #     column_list=['id'],
        #     where_dict=None,
        #     join_rules_str=None,
        #     order_by_column_list=None)

        # record_reference_list = [x[b'id'] for x in record_reference_list]

        record_reference_list = self.db_get_row_ids(
            table_name = 'records',
            where_dict = None)

        records = self.get_records(
            record_reference_list=record_reference_list,
            record_reference_type='raw',
            active=active,
            inactive=inactive)

        return records


    def get_records_with_annotations(self, annotation_type, annotation, active=True, inactive=False):

        sql_string = '''
            SELECT rec_id
            FROM record_annotations
            LEFT OUTER JOIN record_annotation_types ON record_annotation_types.id=rec_ann_type_id
            LEFT OUTER JOIN record_annotation_values ON record_annotation_values.id=rec_ann_value_id
            WHERE type IS ? AND value IS ?;
            '''

        try:
            self._DB_CURSOR.execute(sql_string, [annotation_type, annotation])
        except sqlite3.Error as error:
            print(error, file=sys.stderr)
            self._DB_CONN.rollback()
            raise

        results = self._DB_CURSOR.fetchall()

        rec_id_list = list()

        for result in results:
            rec_id_list.append(result[b'rec_id'])

        records = self.get_records(
            record_reference_list=rec_id_list,
            record_reference_type='raw'  # gi version internal raw
            )

        return records

    def get_parent_rec_ids(self, rec_id):

        # CREATE TABLE record_ancestry(
        #     id INTEGER PRIMARY KEY AUTOINCREMENT,
        #     rec_id INTEGER NOT NULL REFERENCES records(id) ON DELETE CASCADE,
        #     parent_rec_id INTEGER NOT NULL REFERENCES records(id) ON DELETE RESTRICT
        #     );

        results = self.db_select(
            table_name_list=['record_ancestry'],
            column_list=['parent_rec_id'],
            where_dict={'rec_id': rec_id},
            join_rules_str=None,
            order_by_column_list=['parent_rec_id'])

        parent_rec_ids = None

        if len(results) > 0:
            parent_rec_ids = [x[b'parent_rec_id'] for x in results]

        return parent_rec_ids


    ############################################################################


# if __name__ == '__main__':

#     pass

#     import os
#     from krpy import krbioio

#     # db_file = 'sqlite3_database'
#     db_file = '/Users/karolis/Desktop/phylo-prj/db.sqlite3'

#     # if os.path.exists(db_file):
#     #     os.remove(db_file)

#     seq_db = KRSequenceDatabase(db_file)

#     ##########################################################################

#     x = seq_db.db_update(
#         table_name='organisms',
#         values_dict={'subspecies': 'karolis'},
#         where_dict={'genus': 'Solanum'})

#     print(x)

#     ##########################################################################

#     # gb_file = 'records_NADH2.gb'

#     # gb_records = krbioio.read_sequence_file(
#     #     file_path=gb_file,
#     #     file_format='gb',
#     #     ret_type='dict',
#     #     key='gi'
#     # )

#     ##########################################################################

#     # seq_db.save()
#     # seq_db.close()
