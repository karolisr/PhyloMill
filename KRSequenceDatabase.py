#!/usr/bin/env python
# -*- coding: utf-8 -*-


from __future__ import print_function
from __future__ import unicode_literals


class KRSequenceDatabase:


    ############################################################################


    _DB_SCRIPT = '''

    PRAGMA foreign_keys = ON;

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
        taxonomy_id INTEGER REFERENCES taxonomies(id) ON DELETE RESTRICT,
        genus TEXT NOT NULL,
        species TEXT,
        subspecies TEXT,
        variety TEXT,
        hybrid TEXT,
        other TEXT,
        authority TEXT,
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

    CREATE TABLE record_feature_types(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        type TEXT NOT NULL UNIQUE
        );

    CREATE TABLE record_features(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        rec_id INTEGER NOT NULL REFERENCES records(id) ON DELETE CASCADE,
        rec_feat_type_id INTEGER NOT NULL REFERENCES record_feature_types(id) ON DELETE CASCADE,
        location TEXT NOT NULL
        );

    CREATE TABLE record_feature_qualifier_types(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        type TEXT NOT NULL UNIQUE
        );

    CREATE TABLE record_feature_qualifiers(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        rec_feat_id INTEGER NOT NULL REFERENCES record_features(id) ON DELETE CASCADE,
        rec_feat_qual_type_id INTEGER NOT NULL REFERENCES record_feature_qualifier_types(id) ON DELETE CASCADE,
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

    ############################################################################


    class _KRSeqEdits(list):
        pass


    def _produce_seq_edits(self, s1, s2):
        from krpy import krstring
        return self._KRSeqEdits(krstring.produce_edits(s1, s2))


    def _apply_seq_edits(self, e, s):
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
            detect_types=sqlite3.PARSE_DECLTYPES)
        self._DB_CONN.row_factory = sqlite3.Row
        self._DB_CONN.text_factory = str

        self._DB_CONN.execute('PRAGMA foreign_keys = ON;')

        self._DB_CURSOR = self._DB_CONN.cursor()


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


    def _db_select(self, table_name_list, column_list, where_dict=None, join_rules_str=None, order_by_column_list=None):

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


    def _db_get_row_id(self, table_name, values_dict):

        results = self._db_select(
            table_name_list=[table_name],
            column_list=['id'],
            where_dict=values_dict)

        row_id = None

        if len(results) == 1:
            row_id = results[0][b'id']
        elif len(results) > 1:
            row_id = -1

        return row_id


    def _db_insert(self, table_name, values_dict, check_exists=True):

        import sys
        import sqlite3

        row_id = None
        already_in_db = False

        columns_str = ', '.join(values_dict.keys())
        values_str = ':' + ', :'.join(values_dict.keys())

        exists = False
        if check_exists:
            exists = self._db_get_row_id(table_name, values_dict)

        if check_exists and exists:
            row_id = exists
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


    def _db_update(self, table_name, values_dict, where_dict=None):

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


    def _db_delete(self, table_name, where_dict=None):

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


    def _bio_alphabet_to_string(self, bio_alphabet):

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


    def _string_to_bio_alphabet(self, str_alphabet):

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
        return self._bio_alphabet_to_string(bio_alphabet=alphabet)


    ############################################################################


    def add_organism(self, organism_dict, taxonomy_list=None, ncbi_tax_id_list=None, synonymy_check_done=False):

        taxonomy_id = None
        taxonomy_str = None

        if taxonomy_list:
            taxonomy_str = ','.join(taxonomy_list)
            values_dict = {'taxonomy': taxonomy_str}
            taxonomy_id = self._db_insert('taxonomies', values_dict)[0]

        for key in organism_dict.keys():
            if str(organism_dict[key]) == b'':
                organism_dict[key] = None

        if synonymy_check_done:
            synonymy_check_done = 1
        else:
            synonymy_check_done = 0

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
            'taxonomy_id': taxonomy_id,
            'genus': organism_dict['genus'],
            'species': organism_dict['species'],
            'subspecies': organism_dict['subspecies'],
            'variety': organism_dict['variety'],
            'hybrid': organism_dict['hybrid'],
            'other': organism_dict['other'],
            'authority': organism_dict['authority'],
            'synonymy_check_done': synonymy_check_done
        }

        already_in_db = False

        org_id = self._db_get_row_id('organisms', values_dict=where_dict)

        if org_id and org_id >= 0:

            org_dict = self._db_select(
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

            self._db_update('organisms',
                values_dict=values_dict,
                where_dict={'id': org_id})
            already_in_db = True

        else:
            org_id = self._db_insert('organisms', values_dict)[0]

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

                self._db_insert('ncbi_tax_ids', values_dict)

        return (org_id, already_in_db)


    # def update_organism(self, where_dict, organism_dict, taxonomy_list=None, ncbi_tax_id_list=None, synonymy_check_done=False):

    #     taxonomy_id = None
    #     taxonomy_str = None

    #     if taxonomy_list:

    #         taxonomy_str = ','.join(taxonomy_list)
    #         values_dict = {'taxonomy': taxonomy_str}
    #         taxonomy_id = self._db_insert('taxonomies', values_dict)[0]

    #     for key in organism_dict.keys():
    #         if str(organism_dict[key]) == b'':
    #             organism_dict[key] = None

    #     if synonymy_check_done:
    #         synonymy_check_done = 1
    #     else:
    #         synonymy_check_done = 0

    #     values_dict = {
    #         'taxonomy_id': taxonomy_id,
    #         # 'ncbi_tax_id': ncbi_tax_id,
    #         'genus': organism_dict['genus'],
    #         'species': organism_dict['species'],
    #         'subspecies': organism_dict['subspecies'],
    #         'variety': organism_dict['variety'],
    #         'hybrid': organism_dict['hybrid'],
    #         'other': organism_dict['other'],
    #         'authority': organism_dict['authority'],
    #         'synonymy_check_done': synonymy_check_done
    #     }

    #     org_id = self._db_get_row_id(
    #         table_name='organisms', values_dict=where_dict)

    #     self._db_update(
    #         table_name='organisms',
    #         values_dict=values_dict,
    #         where_dict=where_dict)

    #     if ncbi_tax_id_list:

    #         for ncbi_tax_id in ncbi_tax_id_list:

    #             # CREATE TABLE ncbi_tax_ids(
    #             #     id INTEGER PRIMARY KEY AUTOINCREMENT,
    #             #     org_id INTEGER NOT NULL REFERENCES organisms(id) ON DELETE CASCADE,
    #             #     ncbi_tax_id INTEGER NOT NULL

    #             values_dict = {
    #                 'org_id': org_id,
    #                 'ncbi_tax_id': ncbi_tax_id
    #                 }

    #             self._db_insert('ncbi_tax_ids', values_dict)


    def delete_organisms(self, where_dict):

        table_name = 'organisms'

        self._db_delete(
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
                       'taxonomy', 'synonymy_check_done']
        order_by_column_list = ['genus', 'species']

        orgs = self._db_select(
            table_name_list=table_name_list,
            column_list=column_list,
            where_dict=where_dict,
            join_rules_str='organisms.taxonomy_id=taxonomies.id',
            order_by_column_list=order_by_column_list)

        results = list()

        for org in orgs:

            org = dict(org)

            ncbi_tax_ids = self._db_select(
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


    def _add_sequence(self, rec_id, sequence_str, sequence_alphabet_str):

        sequence_str = sequence_str.upper()
        sequence_alphabet_str = sequence_alphabet_str.upper()

        values_dict = {'alphabet': sequence_alphabet_str}
        seq_alpha_id = self._db_get_row_id('sequence_alphabets',
            values_dict)

        values_dict = {
            'rec_id': rec_id,
            'seq_alpha_id': seq_alpha_id,
            'sequence': sequence_str
        }

        row_id = self._db_insert('sequences', values_dict)

        return row_id


    def _add_sequence_representation(self, seq_id, repr_list, rec_id=None, aln_id=None):

        values_dict = {
            'aln_id': aln_id,
            'rec_id': rec_id,
            'seq_id': seq_id,
            'representation': repr_list
        }

        row_id = self._db_insert('sequence_representations', values_dict,
            check_exists=False)

        return row_id


    def add_record_feature(self, rec_id, type_str, location_str):

        values_dict = {'type': type_str}
        rec_feat_type_id = self._db_insert('record_feature_types',
            values_dict)[0]

        values_dict = {
            'rec_id': rec_id,
            'rec_feat_type_id': rec_feat_type_id,
            'location': location_str
        }

        row_id = self._db_insert('record_features', values_dict)

        return row_id


    def add_record_feature_qualifier(self, rec_feat_id, type_str, qualifier_str):

        values_dict = {'type': type_str}
        rec_feat_qual_type_id = self._db_insert(
            'record_feature_qualifier_types',
            values_dict)[0]

        values_dict = {
            'rec_feat_id': rec_feat_id,
            'rec_feat_qual_type_id': rec_feat_qual_type_id,
            'qualifier': qualifier_str
        }

        row_id = self._db_insert('record_feature_qualifiers', values_dict)

        return row_id


    def _add_record_annotation(self, rec_id, type_str, annotation_str):

        values_dict = {'type': type_str}
        rec_ann_type_id = self._db_insert('record_annotation_types',
            values_dict)[0]

        values_dict = {'value': annotation_str}
        rec_ann_value_id = self._db_insert('record_annotation_values',
            values_dict)[0]

        values_dict = {
            'rec_id': rec_id,
            'rec_ann_type_id': rec_ann_type_id,
            'rec_ann_value_id': rec_ann_value_id
        }

        row_id = self._db_insert('record_annotations', values_dict)

        return row_id


    def add_record_annotation(
        self,
        record_reference,
        type_str,
        annotation_str,
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

        rec_id = self._db_get_row_id(
            table_name='records',
            values_dict=where_dict)

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

        row_id = self._db_insert('record_action_history', values_dict)

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

        results = self._db_get_row_id(
            table_name='records',
            values_dict=where_dict)

        if not results:
            return False
        else:
            return True


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

        results = self._db_get_row_id(
            table_name='blacklist',
            values_dict=where_dict)

        if not results:
            return False
        else:
            return True


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

            row_id = self._db_insert('blacklist', values_dict)

        return row_id


    def delete_records(self, where_dict, blacklist=False, blacklist_notes=None):

        table_name = 'records'
        column_list = ['ncbi_gi', 'ncbi_version', 'internal_reference']

        rec_ids_list = self._db_select(
            table_name_list=[table_name],
            column_list=column_list,
            where_dict=where_dict,
            join_rules_str=None)

        self._db_delete(
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

        self._db_update(
            table_name='records',
            values_dict=values_dict,
            where_dict=where_dict)


    def add_genbank_record(self, record, action_str='Genbank record'):

        from krpy import krbionames
        from krpy import krncbi

        # Organism

        org = record.annotations['organism']
        org_dict = krbionames.parse_organism_name(org, sep=' ',
            ncbi_authority=False)
        taxonomy_list = record.annotations['taxonomy']
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

        row_id = self._db_insert('records', values_dict)

        self._add_record_action(
            rec_id=row_id[0],
            action_str=action_str)

        # Sequence

        sequence_str = str(record.seq)
        sequence_alphabet_str = self._sequence_alphabet(record.seq)

        seq_id = self._add_sequence(
            rec_id=row_id[0],
            sequence_str=sequence_str,
            sequence_alphabet_str=sequence_alphabet_str)[0]

        repr_list = self._produce_seq_edits(sequence_str, sequence_str)

        seq_rep_id = self._add_sequence_representation(
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

                rec_feat_qual_id = self.add_record_feature_qualifier(
                    rec_feat_id=rec_feat_id[0],
                    type_str=qualifier,
                    qualifier_str=str(feature.qualifiers[qualifier][0]))

        return row_id


    def _align_sequence_reps(self, seq_rep_id_list, program, options='', program_executable=''):

        from Bio import SeqRecord
        from krpy import kralign

        temp_records = list()

        for seq_rep_id in seq_rep_id_list:
            seq = self._get_sequence_from_representation(seq_rep_id=seq_rep_id)
            seq_record = SeqRecord.SeqRecord(seq, id=str(seq_rep_id), name='', description='')
            temp_records.append(seq_record)

        aln = kralign.align(
            records=temp_records,
            program=program,
            options='',
            program_executable='')

        result_seq_rep_id_list = list()

        for s in aln:
            old_seq_rep_id = int(s.id)

            table_name_list = ['sequence_representations']
            column_list = ['seq_id']
            where_dict = {'id': old_seq_rep_id}
            join_rules_str = None

            results = self._db_select(
                table_name_list=table_name_list,
                column_list=column_list,
                where_dict=where_dict,
                join_rules_str=join_rules_str)[0]

            seq_id = results[b'seq_id']
            seq = self._get_sequence(sequence_id=seq_id)
            aln_seq = str(s.seq).upper()

            repr_list = self._produce_seq_edits(str(seq), aln_seq)

            new_seq_rep_id = self._add_sequence_representation(
                seq_id=seq_id, repr_list=repr_list)[0]

            result_seq_rep_id_list.append(new_seq_rep_id)

        return result_seq_rep_id_list


    def _add_alignment(self, name, seq_rep_id_list, description=None, rec_id=None):

        values_dict = {
            'rec_id': rec_id,
            'name': name,
            'description': description
        }

        row_id = self._db_insert('alignments', values_dict, check_exists=False)

        for seq_rep_id in seq_rep_id_list:

            values_dict = {
                'aln_id': row_id[0]
            }

            where_dict = {
                'id': seq_rep_id
            }

            self._db_update(
                table_name='sequence_representations',
                values_dict=values_dict,
                where_dict=where_dict)


    def align_records(
        self,
        record_reference_list,
        program,
        aln_name,
        options='',
        program_executable='',
        record_reference_type='gi',  # gi version internal
        description=None
        ):

        seq_rep_id_list = list()

        for r in record_reference_list:
            seq_rep_id = self._get_seq_rep_id_for_record(
                record_reference=r,
                record_reference_type=record_reference_type)
            seq_rep_id_list.append(seq_rep_id)

        aln_seq_rep_id_list = self._align_sequence_reps(
            seq_rep_id_list=seq_rep_id_list,
            program=program,
            options=options,
            program_executable=program_executable)

        row_id = self._add_alignment(
            name=aln_name,
            seq_rep_id_list=aln_seq_rep_id_list,
            description=description)

        return row_id


    ############################################################################


    def _get_sequence(self, sequence_id):

        from Bio import Seq

        table_name_list = ['sequences', 'sequence_alphabets']
        column_list = ['alphabet', 'sequence']
        where_dict = {'sequences.id': sequence_id}
        join_rules_str = 'sequences.seq_alpha_id=sequence_alphabets.id'

        results = self._db_select(
            table_name_list=table_name_list,
            column_list=column_list,
            where_dict=where_dict,
            join_rules_str=join_rules_str)[0]

        seq_str = results[b'sequence']
        seq_alphabet = self._string_to_bio_alphabet(str_alphabet=results[b'alphabet'])

        seq = Seq.Seq(data=seq_str, alphabet=seq_alphabet)

        return(seq)


    def _get_sequence_from_representation(self, seq_rep_id):

        from Bio import Seq

        table_name_list = ['sequence_representations']
        column_list = ['seq_id', 'representation']
        where_dict = {'id': seq_rep_id}
        join_rules_str = None

        results = self._db_select(
            table_name_list=table_name_list,
            column_list=column_list,
            where_dict=where_dict,
            join_rules_str=join_rules_str)[0]

        sequence_id = results[b'seq_id']
        seq = self._get_sequence(sequence_id)

        sequence_edits = results[b'representation']

        new_seq_str = self._apply_seq_edits(e=sequence_edits, s=str(seq))
        new_seq = Seq.Seq(data=new_seq_str, alphabet=seq.alphabet)

        return(new_seq)


    def _get_seq_rep_id_for_record(self,
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

        table_name = 'records'
        where_dict = {where_dict_key: record_reference}
        rec_id = self._db_get_row_id(
            table_name=table_name,
            values_dict=where_dict)

        table_name = 'sequence_representations'
        where_dict = {'rec_id': rec_id}
        seq_rep_id = self._db_get_row_id(
            table_name=table_name,
            values_dict=where_dict)

        return seq_rep_id


    def get_sequence_for_record(self,
        record_reference,
        record_reference_type='gi'  # gi version internal
        ):

        seq_rep_id = self._get_seq_rep_id_for_record(
            record_reference=record_reference,
            record_reference_type=record_reference_type)

        seq = self._get_sequence_from_representation(seq_rep_id=seq_rep_id)

        return(seq)


    def get_sequences_for_records(self,
        record_reference_list,
        record_reference_type='gi'  # gi version internal
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
        record_reference_type='gi'  # gi version internal
        ):

        from Bio.SeqRecord import SeqRecord
        from krpy import krbionames

        where_dict_key = ''

        if record_reference_type == 'gi':
            where_dict_key = 'ncbi_gi'
        elif record_reference_type == 'version':
            where_dict_key = 'ncbi_version'
        elif record_reference_type == 'internal':
            where_dict_key = 'internal_reference'

        where_dict = {where_dict_key: record_reference}
        results = self._db_select(
            table_name_list=['records'],
            column_list=['org_id', 'ncbi_gi', 'ncbi_version', 'internal_reference',
                         'description'],
            where_dict=where_dict,
            join_rules_str=None)[0]

        where_dict = {b'organisms.id': results[b'org_id']}

        org_dict = self.get_organisms(where_dict=where_dict)[0]
        org_flat = krbionames.flatten_organism_name(
                parsed_name=org_dict, sep=' ')

        seq_rep_id = self._get_seq_rep_id_for_record(
            record_reference=record_reference,
            record_reference_type=record_reference_type)

        seq = self._get_sequence_from_representation(seq_rep_id=seq_rep_id)

        record = SeqRecord(
            seq=seq,
            id=results[b'ncbi_version'],
            name=results[b'ncbi_version'],
            description=results[b'description'])

        record.annotations[b'gi'] = str(results[b'ncbi_gi'])
        record.annotations[b'organism'] = str(org_flat)

        return record


    def get_records(
        self,
        record_reference_list,
        record_reference_type='gi'  # gi version internal
        ):

        record_list = list()

        for ref in record_reference_list:
            record = self.get_record(
                record_reference=ref,
                record_reference_type=record_reference_type  # gi version internal
                )

            record_list.append(record)

        return record_list


    ############################################################################


if __name__ == '__main__':

    import os
    from krpy import krbioio

    db_file = 'sqlite3_database'

    if os.path.exists(db_file):
        os.remove(db_file)

    seq_db = KRSequenceDatabase(db_file)

    ############################################################################

    # gb_file = 'records_NADH2.gb'

    # gb_records = krbioio.read_sequence_file(
    #     file_path=gb_file,
    #     file_format='gb',
    #     ret_type='dict',
    #     key='gi'
    # )

    ############################################################################

    seq_db.save()
    seq_db.close()
