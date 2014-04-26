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

    CREATE TABLE organisms(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        taxonomy_id INTEGER REFERENCES taxonomies(id),
        ncbi_tax_id INTEGER,
        genus TEXT NOT NULL,
        species TEXT,
        subspecies TEXT,
        variety TEXT,
        hybrid TEXT,
        other TEXT,
        authority TEXT
        );

    CREATE TABLE blacklist(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        ncbi_gi INTEGER,
        ncbi_version TEXT,
        internal_reference TEXT
        );

    CREATE TABLE records(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        org_id INTEGER NOT NULL REFERENCES organisms(id),
        seq_rep_id INTEGER REFERENCES sequence_representations(id),
        aln_id INTEGER REFERENCES alignments(id),
        active INTEGER NOT NULL,
        ncbi_gi INTEGER,
        ncbi_version TEXT,
        internal_reference TEXT,
        description TEXT
        );

    CREATE TABLE record_ancestry(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        rec_id INTEGER NOT NULL REFERENCES records(id),
        parent_rec_id INTEGER NOT NULL REFERENCES records(id)
        );

    CREATE TABLE record_action_history(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        rec_id INTEGER NOT NULL REFERENCES records(id),
        action TEXT NOT NULL
        );

    CREATE TABLE record_feature_types(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        type TEXT NOT NULL UNIQUE
        );

    CREATE TABLE record_features(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        rec_id INTEGER NOT NULL REFERENCES records(id),
        rec_feat_type_id INTEGER NOT NULL REFERENCES record_feature_types(id),
        location TEXT NOT NULL
        );

    CREATE TABLE record_feature_qualifier_types(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        type TEXT NOT NULL UNIQUE
        );

    CREATE TABLE record_feature_qualifiers(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        rec_feat_id INTEGER NOT NULL REFERENCES record_features(id),
        rec_feat_qual_type_id INTEGER NOT NULL REFERENCES record_feature_qualifier_types(id),
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
        rec_id INTEGER NOT NULL REFERENCES records(id),
        rec_ann_type_id INTEGER NOT NULL REFERENCES record_annotation_types(id),
        rec_ann_value_id INTEGER NOT NULL REFERENCES record_annotation_values(id)
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
        seq_alpha_id INTEGER NOT NULL REFERENCES sequence_alphabets(id),
        sequence TEXT NOT NULL
        );

    CREATE TABLE sequence_representations(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        seq_id INTEGER NOT NULL REFERENCES sequences(id),
        representation SEQREP NOT NULL
        );

    CREATE TABLE alignments(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        name TEXT NOT NULL,
        description TEXT
        );

    CREATE TABLE alignment_sequence_representations(
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        aln_id INTEGER NOT NULL REFERENCES alignments(id),
        seq_rep_id INTEGER NOT NULL REFERENCES sequence_representations(id)
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


    def _db_select(self, table_name_list, column_list, where_dict=None, join_rules_str=None):

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

        tables_str = ''

        if len(table_name_list) == 1:
            tables_str = table_name_list[0]
        else:
            tables_str = ' INNER JOIN '.join(table_name_list)
            tables_str = tables_str + ' ON ' + join_rules_str

        select_str = str('SELECT ' + columns_str + ' FROM ' + tables_str + where_str)

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


    ############################################################################


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


    def add_organism(self, organism_dict, taxonomy_list=None, ncbi_tax_id=None):

        # id INTEGER PRIMARY KEY AUTOINCREMENT,
        # taxonomy_id INTEGER REFERENCES taxonomies(id),
        # ncbi_tax_id INTEGER,
        # genus TEXT NOT NULL,
        # species TEXT,
        # subspecies TEXT,
        # variety TEXT,
        # hybrid TEXT,
        # other TEXT,
        # authority TEXT

        taxonomy_id = None
        taxonomy_str = None

        if taxonomy_list:

            taxonomy_str = ','.join(taxonomy_list)
            values_dict = {'taxonomy': taxonomy_str}
            taxonomy_id = self._db_insert('taxonomies', values_dict)[0]

        for key in organism_dict.keys():
            if str(organism_dict[key]) == b'':
                organism_dict[key] = None

        values_dict = {
            'taxonomy_id': taxonomy_id,
            'ncbi_tax_id': ncbi_tax_id,
            'genus': organism_dict['genus'],
            'species': organism_dict['species'],
            'subspecies': organism_dict['subspecies'],
            'variety': organism_dict['variety'],
            'hybrid': organism_dict['cross'],
            'other': organism_dict['other'],
            'authority': organism_dict['authority']
        }

        row_id = self._db_insert('organisms', values_dict)

        return row_id


    def add_sequence(self, sequence_str, sequence_alphabet_str):

        # id INTEGER PRIMARY KEY AUTOINCREMENT,
        # seq_alpha_id INTEGER NOT NULL REFERENCES sequence_alphabets(id),
        # sequence TEXT NOT NULL

        sequence_str = sequence_str.upper()
        sequence_alphabet_str = sequence_alphabet_str.upper()

        values_dict = {'alphabet': sequence_alphabet_str}
        seq_alpha_id = self._db_get_row_id('sequence_alphabets',
            values_dict)

        values_dict = {
            'seq_alpha_id': seq_alpha_id,
            'sequence': sequence_str
        }

        row_id = self._db_insert('sequences', values_dict)

        return row_id


    def add_sequence_representation(self, seq_id, repr_list):

        # id INTEGER PRIMARY KEY AUTOINCREMENT,
        # seq_id INTEGER NOT NULL REFERENCES sequences(id),
        # aln_id INTEGER REFERENCES alignments(id),
        # representation SEQREP NOT NULL

        values_dict = {
            'seq_id': seq_id,
            'representation': repr_list
        }

        row_id = self._db_insert('sequence_representations', values_dict,
            check_exists=False)

        return row_id


    def add_record_feature(self, rec_id, type_str, location_str):

        # id INTEGER PRIMARY KEY AUTOINCREMENT,
        # rec_id INTEGER NOT NULL REFERENCES records(id),
        # rec_feat_type_id INTEGER NOT NULL REFERENCES record_feature_types(id),
        # location TEXT NOT NULL

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

        # id INTEGER PRIMARY KEY AUTOINCREMENT,
        # rec_feat_id INTEGER NOT NULL REFERENCES record_features(id),
        # rec_feat_qual_type_id INTEGER NOT NULL REFERENCES record_feature_qualifier_types(id),
        # qualifier TEXT NOT NULL

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


    def add_record_annotation(self, rec_id, type_str, annotation_str):

        # id INTEGER PRIMARY KEY AUTOINCREMENT,
        # rec_id INTEGER NOT NULL REFERENCES records(id),
        # rec_ann_type_id INTEGER NOT NULL REFERENCES record_annotation_types(id),
        # rec_ann_value_id INTEGER NOT NULL REFERENCES record_annotation_values(id)

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


    def _add_record_action(self, rec_id, action_str):

        # id INTEGER PRIMARY KEY AUTOINCREMENT,
        # rec_id INTEGER NOT NULL REFERENCES records(id),
        # action TEXT NOT NULL

        values_dict = {
            'rec_id': rec_id,
            'action': action_str
        }

        row_id = self._db_insert('record_action_history', values_dict)

        return row_id


    def add_record(self,
        org,  # string or org_dict
        ncbi_gi,
        ncbi_version,
        internal_reference,
        description,
        sequence_str,
        # seq_rep_id=None,
        # aln_id=None,
        taxonomy_list=None,
        ncbi_tax_id=None,
        action_str='New record'):

        from krpy import krbionames

        if isinstance(org, basestring):
            org = krbionames.parse_organism_name(org, sep=' ', ncbi_authority=False)

        # id INTEGER PRIMARY KEY AUTOINCREMENT,
        # org_id INTEGER NOT NULL REFERENCES organisms(id),
        # seq_rep_id INTEGER REFERENCES sequence_representations(id),
        # aln_id INTEGER REFERENCES alignments(id),
        # active INTEGER NOT NULL,
        # ncbi_gi INTEGER,
        # ncbi_version TEXT,
        # internal_reference TEXT,
        # description TEXT

        pass


    def add_genbank_record(self, record, action_str='Genbank record'):

        from krpy import krbionames
        from krpy import krncbi

        # id INTEGER PRIMARY KEY AUTOINCREMENT,
        # org_id INTEGER NOT NULL REFERENCES organisms(id),
        # seq_rep_id INTEGER REFERENCES sequence_representations(id),
        # aln_id INTEGER REFERENCES alignments(id),
        # active INTEGER NOT NULL,
        # ncbi_gi INTEGER,
        # ncbi_version TEXT,
        # internal_reference TEXT,
        # description TEXT

        # Organism

        org = record.annotations['organism']
        org_dict = krbionames.parse_organism_name(org, sep=' ',
            ncbi_authority=False)
        taxonomy_list = record.annotations['taxonomy']
        ncbi_tax_id = int(krncbi.get_ncbi_tax_id(record))

        org_id = self.add_organism(
            organism_dict=org_dict,
            taxonomy_list=taxonomy_list,
            ncbi_tax_id=ncbi_tax_id)[0]

        # Sequence

        sequence_str = str(record.seq)
        sequence_alphabet_str = self._sequence_alphabet(record.seq)

        seq_id = self.add_sequence(
            sequence_str=sequence_str,
            sequence_alphabet_str=sequence_alphabet_str)[0]

        repr_list = self._produce_seq_edits(sequence_str, sequence_str)

        seq_rep_id = self.add_sequence_representation(
            seq_id=seq_id,
            repr_list=repr_list)[0]

        # Record

        active = 1

        ncbi_gi = record.annotations['gi']

        ncbi_version = record.id

        description = record.description

        values_dict = {
            'org_id': org_id,
            'seq_rep_id': seq_rep_id,
            'active': active,
            'ncbi_gi': ncbi_gi,
            'ncbi_version': ncbi_version,
            'description': description
        }

        row_id = self._db_insert('records', values_dict)

        self._add_record_action(
            rec_id=row_id[0],
            action_str=action_str)

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

            results = seq_db._db_select(
                table_name_list=table_name_list,
                column_list=column_list,
                where_dict=where_dict,
                join_rules_str=join_rules_str)[0]

            seq_id = results[b'seq_id']
            seq = self._get_sequence(sequence_id=seq_id)
            aln_seq = str(s.seq).upper()

            repr_list = self._produce_seq_edits(str(seq), aln_seq)

            ####

            # print('--- --- --- ---')
            # print(seq[0:130])
            # print(aln_seq[0:-1])
            # repr_list_str = self._adapt_seq_edits(e=repr_list)
            # print(repr_list_str)
            # repr_list_new = self._convert_seq_edits(e=repr_list_str)
            # print(repr_list)
            # print(repr_list_new)
            # print((self._apply_seq_edits(e=repr_list_new, s=str(seq)))[100:230])
            # print(repr_list)

            ###

            new_seq_rep_id = self.add_sequence_representation(
                seq_id=seq_id, repr_list=repr_list)[0]

            result_seq_rep_id_list.append(new_seq_rep_id)

        return result_seq_rep_id_list


    def _add_alignment(self, name, seq_rep_id_list, description=None):

        # CREATE TABLE sequences(
        #     id INTEGER PRIMARY KEY AUTOINCREMENT,
        #     seq_alpha_id INTEGER NOT NULL REFERENCES sequence_alphabets(id),
        #     sequence TEXT NOT NULL
        #     );

        # CREATE TABLE sequence_representations(
        #     id INTEGER PRIMARY KEY AUTOINCREMENT,
        #     seq_id INTEGER NOT NULL REFERENCES sequences(id),
        #     representation SEQREP NOT NULL
        #     );

        # CREATE TABLE alignments(
        #     id INTEGER PRIMARY KEY AUTOINCREMENT,
        #     name TEXT NOT NULL,
        #     description TEXT
        #     );

        # CREATE TABLE alignment_sequence_representations(
        #     id INTEGER PRIMARY KEY AUTOINCREMENT,
        #     aln_id INTEGER NOT NULL REFERENCES alignments(id),
        #     seq_rep_id INTEGER NOT NULL REFERENCES sequence_representations(id)
        #     );

        values_dict = {
            'name': name,
            'description': description
        }

        row_id = self._db_insert('alignments', values_dict, check_exists=False)

        for seq_rep_id in seq_rep_id_list:

            values_dict = {
                'aln_id': row_id[0],
                'seq_rep_id': seq_rep_id
            }

            self._db_insert(
                'alignment_sequence_representations',
                values_dict,
                check_exists=False)

        return row_id


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

        # SELECT alphabet, sequence FROM sequences INNER JOIN sequence_alphabets ON sequences.seq_alpha_id=sequence_alphabets.id;

        table_name_list = ['sequences', 'sequence_alphabets']
        column_list = ['alphabet', 'sequence']
        where_dict = {'sequences.id': sequence_id}
        join_rules_str = 'sequences.seq_alpha_id=sequence_alphabets.id'

        results = seq_db._db_select(
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

        results = seq_db._db_select(
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

        table_name_list = ['records']
        column_list = ['seq_rep_id']
        where_dict = {where_dict_key: record_reference}
        join_rules_str = None

        results = seq_db._db_select(
            table_name_list=table_name_list,
            column_list=column_list,
            where_dict=where_dict,
            join_rules_str=join_rules_str)[0]

        seq_rep_id = results[b'seq_rep_id']

        return(seq_rep_id)


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



    ############################################################################


if __name__ == '__main__':

    import os
    from krpy import krbioio

    db_file = 'sqlite3_database'

    if os.path.exists(db_file):
        os.remove(db_file)

    seq_db = KRSequenceDatabase(db_file)

    # print(seq_db._DB_CONN.execute('PRAGMA foreign_keys;').fetchall())

    # gb_file = 'records_small.gb'
    # gb_file = 'records_matK.gb'
    gb_file = 'records_NADH2.gb'

    gb_records = krbioio.read_sequence_file(
        file_path=gb_file,
        file_format='gb',
        ret_type='dict',
        key='gi'
    )

    for record in gb_records.values():
        # print(record.id, record.seq[0:100])
        seq_db.add_genbank_record(record=record)

    ############################################################################

    # table_name_list = ['records']
    # column_list = ['id', 'ncbi_gi', 'ncbi_version']
    # where_dict = {'ncbi_gi': 331689645}
    # where_dict = {'active': 1}

    # table_name_list = ['sequences', 'sequence_alphabets']
    # column_list = ['alphabet', 'sequence']
    # where_dict = {}
    # join_rules_str = 'sequences.seq_alpha_id=sequence_alphabets.id'

    # results = seq_db._db_select(
    #     table_name_list=table_name_list,
    #     column_list=column_list,
    #     where_dict=where_dict,
    #     join_rules_str=join_rules_str)

    # print(results)

    # row_id = seq_db._db_get_row_id(
    #     table_name=table_name_list[0],
    #     values_dict=where_dict)

    # print(row_id)

    ############################################################################

    # seq = seq_db._get_sequence(sequence_id=2)
    # print(seq, seq.alphabet)

    ############################################################################

    # new_seq = seq_db._get_sequence_from_representation(seq_rep_id=1)
    # print(new_seq, new_seq.alphabet)

    ############################################################################

    # seq = seq_db.get_sequence_for_record(
    #     record_reference=331689645,
    #     record_reference_type='gi'  # gi version internal
    #     )

    # print(seq, seq.alphabet)

    ############################################################################

    # seq_list = seq_db.get_sequences_for_records(
    #     record_reference_list=[331689645, 14599541],
    #     record_reference_type='gi'  # gi version internal
    #     )

    # for seq_dict in seq_list:
    #     print(seq_dict['reference'], seq_dict['seq'])

    ############################################################################

    # aln = seq_db._align_sequence_reps(
    #     seq_rep_id_list=[1,2,3,4,5,6,7,8],
    #     program='mafft')

    # print(aln)

    # seq_db._add_alignment(name='aln', seq_rep_id_list=aln, description=None)

    ############################################################################

    aln_id = seq_db.align_records(
            record_reference_list=[472278893, 156100661, 74272471, 559797281],
            program='mafft',
            aln_name='test_aln',
            options='',
            program_executable='',
            record_reference_type='gi',  # gi version internal
            description='test_description'
            )

    print(aln_id)

    ############################################################################

    seq_db.save()
    seq_db.close()
