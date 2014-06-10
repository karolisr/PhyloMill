#!/usr/bin/env python
# -*- coding: utf-8 -*-


from __future__ import print_function


################################################################################


if __name__ == '__main__':

    import os
    import sys
    import argparse
    import inspect
    import shutil
    import ConfigParser

    from subprocess import call

    from krpy import KRSequenceDatabase
    from krpy.workflow_functions import wf_phylo as wf
    from krpy import krncbi
    from krpy import krother
    from krpy import krbioio
    from krpy import krbionames
    from krpy import krio
    from krpy import krcl
    from krpy import kralign

    from krpy.krother import write_log

    PS = os.path.sep

    ############################################################################

    PARSER = argparse.ArgumentParser()

    PARSER.add_argument(
        '-p',
        '--project_dir',
        type=unicode,
        help="Prepares clean project directory if it doesn't exist. Otherwise \
              sets the project directory.")

    PARSER.add_argument(
        '-c',
        '--commands',
        type=unicode,
        help='Commands to run.')

    PARSER.add_argument(
        '--gi',
        type=int,
        help='NCBI record GI.')

    PARSER.add_argument(
        '--flat',
        action='store_true',
        help='')

    PARSER.add_argument(
        '--raw',
        action='store_true',
        help='')

    PARSER.add_argument(
        '--locus',
        type=unicode,
        help='')

    ARGS = PARSER.parse_args()

    ############################################################################

    COMMANDS = None
    if ARGS.commands:
        COMMANDS = set([x.strip() for x in ARGS.commands.split(',')])

    ############################################################################

    GI = None
    if ARGS.gi:
        GI = ARGS.gi

    LOCUS = None
    if ARGS.locus:
        LOCUS = ARGS.locus

    FLAT = None
    if ARGS.flat:
        FLAT = ARGS.flat

    RAW = None
    if ARGS.raw:
        RAW = ARGS.raw

    ############################################################################

    DB = None
    PRJ_DIR_PATH = None

    # Prepare clean project directory
    if ARGS.project_dir:
        PRJ_DIR_PATH = ARGS.project_dir.rstrip(PS)
        PRJ_DIR_PATH = PRJ_DIR_PATH + PS
    else:
        write_log('Project directory is required.', LFP, newlines_before=1,
            newlines_after=0)
        sys.exit(0)

    # Log file path
    LFP = PRJ_DIR_PATH + 'log.txt'

    if ARGS.project_dir:

        if os.path.exists(PRJ_DIR_PATH):

            msg = 'Using project directory at ' + PRJ_DIR_PATH.rstrip(PS)
            write_log(msg=msg, log_file_path=LFP, append=True,
                newlines_before=1, newlines_after=0)
            DB = KRSequenceDatabase.KRSequenceDatabase(
                PRJ_DIR_PATH + 'db.sqlite3')

        else:

            # Script filename
            script_file_path = inspect.getfile(inspect.currentframe())

            # Script directory path
            script_dir_path = os.path.dirname(os.path.abspath(script_file_path))

            prj_template_dir_path = script_dir_path.strip('workflows') + \
                'data' + PS + 'phylo_prj_template'

            shutil.copytree(prj_template_dir_path, PRJ_DIR_PATH, symlinks=False,
                            ignore=None)

            print('')
            msg = 'Creating project directory at ' + PRJ_DIR_PATH.rstrip(PS)
            write_log(msg=msg, log_file_path=LFP, append=False,
                newlines_before=0, newlines_after=0)

            KRSequenceDatabase.KRSequenceDatabase(PRJ_DIR_PATH + 'db.sqlite3')

            wd = os.getcwd()
            os.chdir(PRJ_DIR_PATH+'organism_name_files')
            call(['./get_ncbi_data.sh'], stdout=open(os.devnull, 'wb'))
            os.chdir(wd)

            sys.exit(0)

    ############################################################################

    DNLD_DIR_PATH = PRJ_DIR_PATH + 'downloaded_files' + PS
    OUT_DIR_PATH = PRJ_DIR_PATH + 'output' + PS
    ORG_LOC_DIR_PATH = OUT_DIR_PATH + 'flatten' + PS
    SRCH_DIR_PATH = PRJ_DIR_PATH + 'search_strategies' + PS
    ORGN_DIR_PATH = PRJ_DIR_PATH + 'organism_name_files' + PS
    TEMP_DIR_PATH = PRJ_DIR_PATH + 'temporary_files' + PS

    CFG_FILE_PATH = PRJ_DIR_PATH + 'config'

    ############################################################################

    if not COMMANDS:
        write_log('No commands given.', LFP, newlines_before=1,
            newlines_after=0)
        sys.exit(0)

    ############################################################################

    # Constants
    FLAT_ID = 0.97

    ############################################################################

    # Get configuration information
    CFG = ConfigParser.SafeConfigParser(allow_no_value=True)
    CFG.optionxform=str
    CFG.read(CFG_FILE_PATH)

    EMAIL = CFG.get('General', 'email')
    MAX_SEQ_LENGTH = CFG.getint('General', 'max_seq_length')

    # Flatten options
    FLAT_ALN_PROG = CFG.get('Flatten', 'align_program')
    FLAT_ALN_PROG_EXE = CFG.get('General', FLAT_ALN_PROG + '_executable')
    FLAT_ALN_PROG_OPTIONS = CFG.get('Flatten', 'align_program_options')
    FLAT_RESOLVE_AMBIGUITIES = CFG.getboolean('Flatten', 'resolve_ambiguities')

    # Hacks
    hacks_items = CFG.items('Hacks')
    hacks_temp = list()
    for hack in hacks_items:
        h = hack[0].split('_')
        if h[-1] != 'data':
            if CFG.getboolean('Hacks', '_'.join(h)):
                hacks_temp.append('_'.join(h))

    HACKS = dict()

    for hack in hacks_temp:
        if CFG.has_option('Hacks', hack + '_data'):
            HACKS[hack] = CFG.get('Hacks', hack + '_data')
        else:
            HACKS[hack] = None

    # Taxa
    tax_temp = CFG.items('Taxa')
    tax_temp = [x[0] for x in tax_temp]
    TAX_IDS = list()
    for tax in tax_temp:
        if tax.isalpha():
            tax_id = list(krncbi.esearch(tax, 'taxonomy', EMAIL))[0]
            TAX_IDS.append(str(tax_id))
            msg = 'NCBI taxonomy ID for ' + tax + ' is ' + str(tax_id)
            write_log(msg, LFP, newlines_before=0, newlines_after=0)
        else:
            TAX_IDS.append(tax)

    # Organism name resolution
    SYN = list()
    syn_temp = None
    if CFG.has_option('Organism Names', 'synonymize'):
        syn_temp = CFG.get('Organism Names', 'synonymize')
    if syn_temp:
        SYN = syn_temp.split(',')
        SYN = [x.strip().lower() for x in SYN]
    SKIP_PREV_SYN = CFG.getboolean('Organism Names',
        'skip_previously_synonymized')
    REM_HYBRIDS = CFG.getboolean('Organism Names', 'remove_hybrids')
    REM_NONSPECIFIC = CFG.getboolean('Organism Names', 'remove_nonspecific')

    # Loci
    LOCI = dict()
    loci_temp = CFG.items('Loci')
    loci_temp = [x[0] for x in loci_temp]
    for l in loci_temp:
        locus_ss_file_path = SRCH_DIR_PATH + l
        if os.path.exists(locus_ss_file_path):
            locus_cfg = ConfigParser.SafeConfigParser(allow_no_value=True)
            locus_cfg.optionxform=str
            locus_cfg.read(locus_ss_file_path)

            l_dict = dict()
            l_strategies = list()

            l_dict['name'] = l

            for l_sec in locus_cfg.sections():
                if l_sec == l:
                    l_db = locus_cfg.get(l, 'database')
                    l_query = locus_cfg.get(l, 'query')
                    # l_short_name = locus_cfg.get(l, 'short_name')
                    l_dict['database'] = l_db
                    l_dict['query'] = l_query
                    # l_dict['short_name'] = l_short_name
                    # l_bad_gis = locus_cfg.get(l_sec, 'bad_gis')
                    # l_bad_gis = l_bad_gis.split(',')
                    # l_bad_gis = [int(x.strip()) for x in l_bad_gis]
                    # l_dict['bad_gis'] = l_bad_gis
                    # l_dict['min_length'] = 0
                    # if locus_cfg.has_option(section=l_sec, option='min_length'):
                    #     l_dict['min_length'] = locus_cfg.getint(l_sec, 'min_length')
                    # l_dict['max_length'] = 0
                    # if locus_cfg.has_option(section=l_sec, option='max_length'):
                    #     l_dict['max_length'] = locus_cfg.getint(l_sec, 'max_length')
                # elif l_sec == 'blast':
                #     l_blast_dbs = locus_cfg.get(l_sec, 'databases')
                #     l_blast_dbs = l_blast_dbs.split(',')
                #     l_blast_dbs = [x.strip() for x in l_blast_dbs]
                #     l_dict['blast_dbs'] = l_blast_dbs
                else:
                    l_lrp = locus_cfg.getint(l_sec, 'locus_relative_position')
                    l_ft = locus_cfg.get(l_sec, 'feature_type')
                    l_ql = locus_cfg.get(l_sec, 'qualifier_label')
                    l_qv = locus_cfg.get(l_sec, 'qualifier_value')
                    l_re = locus_cfg.getboolean(l_sec, 'regex')
                    l_svm = locus_cfg.getboolean(l_sec, 'strict_value_match')
                    l_ml = locus_cfg.getint(l_sec, 'min_length')
                    l_el = locus_cfg.getint(l_sec, 'extra_length')

                    l_strat_dict = {
                    'locus_relative_position': l_lrp,
                    'feature_type': l_ft,
                    'qualifier_label': l_ql,
                    'qualifier_value': l_qv,
                    'regex': l_re,
                    'strict_value_match': l_svm,
                    'min_length': l_ml,
                    'extra_length': l_el
                    }

                    l_strategies.append(l_strat_dict)

            l_dict['strategies'] = l_strategies

            LOCI[l] = l_dict
        else:
            write_log('There is no search strategy file for locus: ' + l, LFP)
            sys.exit(0)

    ############################################################################

    # Search genbank
    if 'search' in COMMANDS:

        for locus_name in LOCI.keys():

            # ln_split = locus_name.split('_')
            # if len(ln_split) > 1 and ln_split[1] == 'blast':

            #     wf.blast_search(
            #         kr_seq_db_object=DB,
            #         loci=LOCI,
            #         locus_name=locus_name,
            #         ncbi_tax_ids=TAX_IDS,
            #         max_seq_length=MAX_SEQ_LENGTH,
            #         email=EMAIL,
            #         log_file_path=LFP,
            #         dnld_dir_path=DNLD_DIR_PATH,
            #         temp_dir=TEMP_DIR_PATH)

            # else:

            wf.regular_search(
                kr_seq_db_object=DB,
                log_file_path=LFP,
                email=EMAIL,
                loci=LOCI,
                locus_name=locus_name,
                ncbi_tax_ids=TAX_IDS,
                max_seq_length=MAX_SEQ_LENGTH,
                dnld_dir_path=DNLD_DIR_PATH)

    ############################################################################

    # Resolve organism names
    if 'resolve_org_names' in COMMANDS:

        ########################################################################

        msg = 'Resolving organism names.'
        write_log(msg, LFP, newlines_before=1, newlines_after=0)

        ########################################################################

        ncbi_names_table = None
        synonymy_table = None

        if len(SYN) > 0:

            msg = 'Loading synonymy table and NCBI organism name list.'
            write_log(msg, LFP, newlines_before=1, newlines_after=0)

            ncbi_names_table = krio.read_table_file(
                path=ORGN_DIR_PATH + 'ncbi_tax_names',
                has_headers=False,
                headers=('tax_id', 'name_txt', 'unique_name', 'name_class'),
                delimiter='\t|',
                quotechar=None,
                stripchar='"',
                rettype='dict')

            synonymy_table = krio.read_table_file(
                path=ORGN_DIR_PATH + 'synonymy.csv',
                has_headers=True, headers=None, delimiter=',')

        ########################################################################

        taxid_blacklist_set = krio.read_table_file(
            path=ORGN_DIR_PATH + 'taxid_blacklist.tsv',
            has_headers=False,
            headers=None,
            delimiter=',',
            quotechar='"',
            rettype='set')

        ########################################################################

        record_taxon_mappings_list = krio.read_table_file(
            path=ORGN_DIR_PATH + 'record_taxon_mappings.tsv',
            has_headers=False,
            headers=('id', 'taxon'),
            delimiter='\t',
            quotechar=None,
            stripchar='"',
            rettype='dict')
        record_taxon_mappings_dict = dict()
        for t in record_taxon_mappings_list:
            record_taxon_mappings_dict[t['id']]=t['taxon']

        ########################################################################

        taxid_taxon_mappings_list = krio.read_table_file(
            path=ORGN_DIR_PATH + 'taxid_taxon_mappings.tsv',
            has_headers=False,
            headers=('taxid', 'taxon'),
            delimiter='\t',
            quotechar=None,
            stripchar='"',
            rettype='dict')
        taxid_taxon_mappings_dict = dict()
        for t in taxid_taxon_mappings_list:
            taxid_taxon_mappings_dict[t['taxid']]=t['taxon']

        ########################################################################

        taxonomy_cache = dict()

        ########################################################################

        if 'tgrc' in HACKS.keys():

            msg = 'Using "Tomato Genetics Resource Center" to resolve organism names.'
            write_log(msg, LFP, newlines_before=1, newlines_after=0)

            wf.rename_tgrc_organisms(
                kr_seq_db_object=DB,
                taxonomy_cache=taxonomy_cache,
                log_file_path=LFP,
                email=EMAIL
                )

        ########################################################################

        msg = 'Checking for record -> taxon mappings.'
        write_log(msg, LFP, newlines_before=1, newlines_after=0)

        wf.rename_organisms_with_record_taxon_mappings(
            kr_seq_db_object=DB,
            record_taxon_mappings_dict=record_taxon_mappings_dict,
            taxonomy_cache=taxonomy_cache,
            log_file_path=LFP,
            email=EMAIL
            )

        ########################################################################

        msg = 'Checking for tax_id -> taxon mappings and synonymy.'
        write_log(msg, LFP, newlines_before=1, newlines_after=0)

        wf.rename_organisms_using_taxids(
            kr_seq_db_object=DB,
            taxid_blacklist_set=taxid_blacklist_set,
            taxid_taxon_mappings_dict=taxid_taxon_mappings_dict,
            taxonomy_cache=taxonomy_cache,
            log_file_path=LFP,
            email=EMAIL,
            rem_hybrids=REM_HYBRIDS,
            rem_nonspecific=REM_NONSPECIFIC,
            skip_prev_syn=SKIP_PREV_SYN,
            tax_groups_to_syn=SYN,
            authority_alternates_file = ORGN_DIR_PATH + \
                'authority_alternates.dat',
            ncbi_names_table=ncbi_names_table,
            synonymy_table=synonymy_table)

        ########################################################################

        krcl.clear_line()
        msg = 'Organism name check: done.'
        write_log(msg, LFP, newlines_before=1, newlines_after=0)

    ############################################################################

    # Extract loci
    if 'extract_loci' in COMMANDS:

        msg = 'Extracting loci.'
        write_log(msg, LFP, newlines_before=1, newlines_after=0)

        for locus_name in LOCI.keys():

            ln_split = locus_name.split('_')
            if len(ln_split) > 1 and ln_split[1] == 'blast':
                continue

            msg = locus_name
            write_log(msg, LFP, newlines_before=1, newlines_after=0)

            locus_dict = LOCI[locus_name]

            records = DB.get_records_with_annotations(
                annotation_type='locus',
                annotation=locus_name,
                active=True,
                inactive=False)

            acc_rej_gi_dict = wf.extract_loci(
                locus_dict=locus_dict,
                records=records,
                log_file_path=LFP,
                kr_seq_db_object=DB,
                temp_dir=TEMP_DIR_PATH)

            ###
            # acc_gi_list = acc_rej_gi_dict['accept']
            # for acc in acc_gi_list:
            #     print(acc)
            ###

            rej_gi_list = acc_rej_gi_dict['reject']
            delete_note = 'failed sequence similarity test'
            rej_gi_list = [[int(x[1]), delete_note] for x in rej_gi_list]

            no_feature_gi_list = acc_rej_gi_dict['no_feature']
            delete_note = 'no locus annotation'
            no_feature_gi_list = [[x, delete_note] for x in no_feature_gi_list]

            bad_list = rej_gi_list + no_feature_gi_list

            print()

            for bad in bad_list:

                bad_gi = bad[0]
                delete_note = bad[1]

                msg = 'inactivating:' + \
                ' gi:' + str(bad_gi) + ' note:' + delete_note
                write_log(msg, LFP)

                where_dict = {'ncbi_gi': bad_gi}

                blacklist_notes = delete_note

                rec_ids = DB.db_get_row_ids(
                    'records',
                    where_dict=where_dict)

                DB.set_inactive(
                    table_name='records',
                    where_dict=where_dict)

                if rec_ids:

                    for rec_id in rec_ids:

                        rec = DB.get_record(
                            record_reference=rec_id,
                            record_reference_type='raw'
                            )

                        DB.add_record_to_blacklist(
                            ncbi_gi=int(rec.annotations['gi']),
                            ncbi_version=rec.id,
                            internal_reference=rec.annotations['internal_reference'],
                            notes=blacklist_notes)

                # DB.delete_records(
                #     where_dict=where_dict,
                #     blacklist=True,
                #     blacklist_notes=blacklist_notes)

                # DB.delete_orphaned_organisms()
                # DB.delete_orphaned_taxonomies()
                DB.save()

    ############################################################################

    # Export active records
    if 'export_active_records' in COMMANDS:

        msg = 'Exporting active records.'
        write_log(msg, LFP, newlines_before=1, newlines_after=0)

        timestamp = krother.timestamp()
        timestamp = timestamp.replace('-', '_')
        timestamp = timestamp.replace(':', '_')
        timestamp = timestamp.replace(' ', '_')

        gb_file_name = ''

        records = None
        if (not LOCUS) and (FLAT or RAW):
            records = DB.get_all_records(active=True, inactive=False)
            gb_file_name = 'active_records' + '_' + timestamp + '.gb'

        elif LOCUS and FLAT:
            records = DB.get_records_with_annotations(
                annotation_type='locus_flat',
                annotation=LOCUS,
                active=True,
                inactive=False)

            for rec in records:
                rec.description = rec.annotations['organism']

            gb_file_name = 'active_records' + '_' + 'locus_flat' + '_' + LOCUS + '_' + timestamp + '.gb'

        elif LOCUS and RAW:
            records = DB.get_records_with_annotations(
                annotation_type='locus',
                annotation=LOCUS,
                active=True,
                inactive=False)
            gb_file_name = 'active_records' + '_' + 'locus_raw' + '_' + LOCUS + '_' + timestamp + '.gb'

        ar_file_path = OUT_DIR_PATH + gb_file_name

        krbioio.write_sequence_file(
            records=records,
            file_path=ar_file_path,
            file_format='genbank')

    ############################################################################

    # Blacklist gi
    if 'blacklist_gi' in COMMANDS:

        if not GI:
            write_log('No GI given.', LFP, newlines_before=1,
                newlines_after=0)
            sys.exit(0)

        record = DB.get_record(
            record_reference=GI,
            record_reference_type='gi'  # gi version internal raw
            )

        del_rec_id = int(record.annotations['kr_seq_db_id'])

        blacklist_notes = 'user_deleted'

        msg = 'inactivating:' + \
        ' gi:' + str(GI) + ' note:' + blacklist_notes
        write_log(msg, LFP)

        where_dict = {'ncbi_gi': GI}

        DB.set_inactive(
            table_name='records',
            where_dict=where_dict)

        DB.add_record_to_blacklist(
            ncbi_gi=GI,
            ncbi_version=record.id,
            internal_reference=record.annotations['internal_reference'],
            notes=blacklist_notes)

        where_dict = {'parent_rec_id': del_rec_id}

        DB.db_delete(
            table_name='record_ancestry',
            where_dict=where_dict)

        DB.save()

    ############################################################################

    # Whitelist gi
    if 'whitelist_gi' in COMMANDS:

        if not GI:
            write_log('No GI given.', LFP, newlines_before=1,
                newlines_after=0)
            sys.exit(0)

        whitelist_notes = 'user_activated'

        msg = 'activating:' + \
        ' gi:' + str(GI) + ' note:' + whitelist_notes
        write_log(msg, LFP)

        where_dict = {'ncbi_gi': GI}

        DB.set_active(
            table_name='records',
            where_dict=where_dict)

        DB.db_delete(
            table_name='blacklist',
            where_dict=where_dict)

        DB.save()

    ############################################################################

    # Produce one locus per organism
    if 'flatten' in COMMANDS:

        msg = 'Producing one locus per organism.'
        write_log(msg, LFP, newlines_before=1, newlines_after=0)

        for locus_name in LOCI.keys():

            msg = locus_name
            write_log(msg, LFP, newlines_before=1, newlines_after=0)

            locus_dir_path = ORG_LOC_DIR_PATH + locus_name + PS
            krio.prepare_directory(locus_dir_path)

            locus_dict = LOCI[locus_name]

            # Get all unflattened records for current locus
            records = DB.get_records_with_annotations(
                annotation_type='locus',
                annotation=locus_name,
                active=True,
                inactive=False)

            records = sorted(records, key=lambda x: len(x.seq), reverse=True)
            reference_records = records[0:min(100, len(records))]

            org_locus_dict = dict()
            for record in records:
                org_id = record.annotations['kr_seq_db_org_id']
                if org_id not in org_locus_dict.keys():
                    org_locus_dict[org_id] = list()
                org_locus_dict[org_id].append(record)

            # Get all flattened records for current locus
            records_flat = DB.get_records_with_annotations(
                annotation_type='locus_flat',
                annotation=locus_name,
                active=True,
                inactive=False)

            org_locus_flat_dict = dict()
            for record_flat in records_flat:
                org_id = record_flat.annotations['kr_seq_db_org_id']
                if org_id not in org_locus_flat_dict.keys():
                    org_locus_flat_dict[org_id] = list()
                org_locus_flat_dict[org_id].append(record_flat)

            ####################################################################

            org_count = len(org_locus_dict.keys())
            for i, org_id in enumerate(org_locus_dict.keys()):

                krcl.print_progress(
                    current=i+1, total=org_count, length=0,
                    prefix=krother.timestamp() + ' - ',
                    postfix='',  # ' - ' + org_flat,
                    show_bar=False)

                where_dict = {'organisms.id': int(org_id)}
                org_dict = DB.get_organisms(where_dict=where_dict)[0]
                org_flat = krbionames.flatten_organism_name(
                    parsed_name=org_dict, sep=' ')

                aln_name = locus_name + '__' + org_flat.replace(' ', '_')

                aln_file_path = locus_dir_path + aln_name + '.phy'
                seq_file_path = locus_dir_path + aln_name + '.fasta'

                org_records = org_locus_dict[org_id]

                ################################################################

                deleted_rec_ids = False
                new_rec_ids = False

                seq_match = True
                aln_match = True

                seq_file_exists = False
                aln_file_exists = False

                ################################################################

                # Check if there already is a flat record for this organism
                if org_id in org_locus_flat_dict:
                    # print('flat_record_exists', org_flat)
                    flat_rec_list = org_locus_flat_dict[org_id]
                    flat_rec = flat_rec_list[0]
                    flat_rec_id = int(flat_rec.annotations['kr_seq_db_id'])

                    db_aln = None

                    if len(flat_rec_list) > 1:
                        print('More than one flat record!')
                    else:

                        parent_rec_ids = DB.get_parent_rec_ids(rec_id=flat_rec_id)

                        rec_ids_in_db = [int(x.annotations['kr_seq_db_id']) for x in org_records]
                        rec_ids_in_db = sorted(rec_ids_in_db, reverse=False)

                        # Check if there is a locus alignment/sequence file
                        seq_file_exists = os.path.exists(seq_file_path)
                        aln_file_exists = os.path.exists(aln_file_path)

                        existing_record_list = list()

                        if seq_file_exists:
                            existing_record_list = krbioio.read_sequence_file(
                                file_path=seq_file_path,
                                file_format='fasta',
                                ret_type='list',
                                key='gi')

                            # os.remove(seq_file_path)

                        elif aln_file_exists:
                            existing_record_list = krbioio.read_alignment_file(
                                file_path=aln_file_path,
                                file_format='phylip-relaxed')

                            # os.remove(aln_file_path)

                        gis_in_file = list()
                        if existing_record_list:
                            gis_in_file = [int(x.id) for x in existing_record_list]

                        rec_id_gi_dict = dict()
                        rec_ids_in_file = list()
                        for gi_in_file in gis_in_file:
                            rec_id_in_file = DB.db_get_row_ids(
                                table_name='records',
                                where_dict={'ncbi_gi': gi_in_file})[0]
                            rec_ids_in_file.append(rec_id_in_file)
                            rec_id_gi_dict[str(rec_id_in_file)] = gi_in_file
                        rec_ids_in_file = sorted(rec_ids_in_file, reverse=False)

                        ########################################################

                        # Check if there are any deleted records
                        if seq_file_exists or aln_file_exists:
                            if set(rec_ids_in_file) == set(parent_rec_ids):
                                # print('no records were deleted')
                                pass
                            else:
                                deleted_rec_ids = list(set(rec_ids_in_file) - set(parent_rec_ids))
                                # if not (len(deleted_rec_ids) > 0):
                                deleted_rec_ids = deleted_rec_ids + list(set(parent_rec_ids) - set(rec_ids_in_file))
                                # print(len(deleted_rec_ids), 'deleted_rec_ids:', deleted_rec_ids)

                        # Check if there are any new records
                        if set(rec_ids_in_db) == set(parent_rec_ids):
                            # print('no new records')
                            pass
                        else:
                            new_rec_ids = list(set(rec_ids_in_db) - set(parent_rec_ids))
                            # print(len(new_rec_ids), 'new_rec_ids:', new_rec_ids)

                        # print(len(parent_rec_ids), 'parent_rec_ids:', parent_rec_ids)
                        # print(len(rec_ids_in_db), 'rec_ids_in_db:', rec_ids_in_db)
                        # print(len(rec_ids_in_file), 'rec_ids_in_file:', rec_ids_in_file)
                        # print('gis_in_file', gis_in_file)

                        ########################################################

                        # Check if there is an alignment associated with this
                        # record
                        aln_id = DB.db_get_row_ids(
                            table_name='alignments',
                            where_dict={'rec_id': flat_rec_id})

                        if aln_id:
                            aln_id = aln_id[0]
                            db_aln = DB.get_alignment(
                                alignment_id=aln_id)

                        # print('aln_id', aln_id)

                        # If there is a locus alignment/sequence file,
                        # get all the sequences from this file and compare to
                        # the sequences in the database
                        if (seq_file_exists or aln_file_exists) and existing_record_list:
                            if db_aln:
                                for a in db_aln:
                                    if str(a.id) in rec_id_gi_dict.keys():
                                        gi = str(rec_id_gi_dict[str(a.id)])
                                        for er in existing_record_list:
                                            if er.id == gi:
                                                er_seq_str = str(er.seq).upper()
                                                db_seq_str = str(a.seq).upper()
                                                if er_seq_str != db_seq_str:
                                                    aln_match = False
                                # print('aln_match', aln_match)
                            else:
                                er = existing_record_list[0]
                                er_seq_str = str(er.seq).upper()
                                db_seq_str = str(flat_rec.seq).upper()
                                if er_seq_str != db_seq_str:
                                    seq_match = False
                                # print('seq_match', seq_match)

                    ############################################################

                    if deleted_rec_ids:
                        for del_rec_id in deleted_rec_ids:
                            del_rec = None

                            for r in org_records:
                                if int(r.annotations['kr_seq_db_id']) == del_rec_id:
                                    del_rec = r
                                    break

                            if del_rec:
                                org_records.remove(del_rec)

                                blacklist_notes = 'user_deleted'

                                ncbi_gi=int(del_rec.annotations['gi'])

                                #### TODO: USED IN blacklist_gi SHOULD BE REFACTORED ####

                                msg = 'inactivating:' + \
                                ' gi:' + str(ncbi_gi) + ' note:' + blacklist_notes
                                write_log(msg, LFP)

                                where_dict = {'id': del_rec_id}

                                DB.set_inactive(
                                    table_name='records',
                                    where_dict=where_dict)

                                DB.add_record_to_blacklist(
                                    ncbi_gi=ncbi_gi,
                                    ncbi_version=del_rec.id,
                                    internal_reference=del_rec.annotations['internal_reference'],
                                    notes=blacklist_notes)

                                where_dict = {'parent_rec_id': del_rec_id}

                                DB.db_delete(
                                    table_name='record_ancestry',
                                    where_dict=where_dict)

                                ####

                    rec_ids_in_db = [int(x.annotations['kr_seq_db_id']) for x in org_records]
                    rec_ids_in_db = sorted(rec_ids_in_db, reverse=False)
                    parent_rec_id_list = rec_ids_in_db

                    new_seq_str = None
                    flat_locus_produced = False
                    if (len(org_records) > 1) and ((not aln_match) or new_rec_ids or deleted_rec_ids):

                        if seq_file_exists:
                            os.remove(seq_file_path)
                        if aln_file_exists:
                            os.remove(aln_file_path)

                        new_aln = existing_record_list
                        if new_rec_ids or deleted_rec_ids:
                            records_to_align = list(existing_record_list)
                            for r in org_records:
                                if new_rec_ids and int(r.annotations['kr_seq_db_id']) in new_rec_ids:
                                    trimmed_rec = wf.trim_record_to_locus(
                                        record=r,
                                        locus_name=locus_name)
                                    records_to_align.append(trimmed_rec)

                            new_aln = wf.flatten_locus(
                                records=records_to_align,
                                reference_records=reference_records,
                                locus_dict=locus_dict,
                                log_file_path=LFP,
                                already_trimmed=True,
                                aln_program=FLAT_ALN_PROG,
                                aln_program_executable=FLAT_ALN_PROG_EXE,
                                aln_options=FLAT_ALN_PROG_OPTIONS,
                                min_locus_sequence_identity=FLAT_ID)

                        wf.update_record_alignment(
                            rec_id=flat_rec_id,
                            new_aln=new_aln,
                            aln_name=aln_name,
                            kr_seq_db_object=DB)

                        consensus = kralign.consensus(
                            alignment=new_aln,
                            threshold=0.4,
                            unknown='N',
                            resolve_ambiguities=FLAT_RESOLVE_AMBIGUITIES)

                        new_seq_str = str(consensus[0]).upper()

                        flat_locus_produced = True

                        krbioio.write_alignment_file(
                            alignment=new_aln,
                            file_path=aln_file_path,
                            file_format='phylip-relaxed')

                    elif (len(org_records) == 1):

                        DB.db_delete(
                            table_name='alignments',
                            where_dict={'rec_id': flat_rec_id})

                        new_rec = None

                        if not seq_match:
                            er = existing_record_list[0]
                            new_seq_str = str(er.seq).upper()
                            new_rec = er

                        elif deleted_rec_ids:
                            trimmed_rec = wf.trim_record_to_locus(
                                record=org_records[0],
                                locus_name=locus_name)
                            if trimmed_rec:
                                new_seq_str = str(trimmed_rec.seq)
                                new_rec = trimmed_rec
                                db_seq_str = str(flat_rec.seq).upper()
                                if new_seq_str.upper() != db_seq_str:
                                    seq_match = False

                        if not seq_match:

                            flat_locus_produced = True

                            if seq_file_exists:
                                os.remove(seq_file_path)
                            if aln_file_exists:
                                os.remove(aln_file_path)

                            krbioio.write_sequence_file(
                                records=[new_rec],
                                file_path=seq_file_path,
                                file_format='fasta')

                    if flat_locus_produced:

                        where_dict = {'rec_id': flat_rec_id}
                        DB.db_delete(
                            table_name='record_ancestry',
                            where_dict=where_dict)

                        for parent_rec_id in parent_rec_id_list:
                            values_dict = {
                                'rec_id': flat_rec_id,
                                'parent_rec_id': parent_rec_id
                                }
                            DB.db_insert('record_ancestry', values_dict)

                        DB.db_delete(
                            table_name='sequences',
                            where_dict={'rec_id': flat_rec_id})

                        sequence_alphabet_str = DB.bio_alphabet_to_string(
                            bio_alphabet=flat_rec.seq.alphabet)

                        seq_id = DB.add_sequence(
                            rec_id=flat_rec_id,
                            sequence_str=new_seq_str,
                            sequence_alphabet_str=sequence_alphabet_str)[0]

                        repr_list = DB.produce_seq_edits(new_seq_str, new_seq_str)

                        seq_rep_id = DB.add_sequence_representation(
                            rec_id=flat_rec_id,
                            seq_id=seq_id,
                            repr_list=repr_list)[0]

                    if len(org_records) == 0:

                        if seq_file_exists:
                            os.remove(seq_file_path)
                        if aln_file_exists:
                            os.remove(aln_file_path)

                        DB.db_delete(
                            table_name='records',
                            where_dict={'id': flat_rec_id})

                    ############################################################

                org_records_count = len(org_records)
                if (org_records_count > 0) and (org_id not in org_locus_flat_dict):

                    # print('flat_record_does_not_exist')

                    flat_locus_produced = False

                    seq_rep_id_list = list()
                    # print(org_id, org_flat, org_records_count)

                    new_seq_str = ''
                    sequence_alphabet_str = DB.bio_alphabet_to_string(
                        bio_alphabet=org_records[0].seq.alphabet)

                    ############################################################

                    if org_records_count > 1:

                        aln = wf.flatten_locus(
                            records=org_records,
                            reference_records=reference_records,
                            locus_dict=locus_dict,
                            log_file_path=LFP,
                            aln_program=FLAT_ALN_PROG,
                            aln_program_executable=FLAT_ALN_PROG_EXE,
                            aln_options=FLAT_ALN_PROG_OPTIONS,
                            min_locus_sequence_identity=FLAT_ID)

                        if aln:
                            flat_locus_produced = True

                            krbioio.write_alignment_file(
                                alignment=aln,
                                file_path=aln_file_path,
                                file_format='phylip-relaxed')

                            for ar in aln:

                                seq = DB.get_sequence_for_record(
                                    record_reference=int(ar.id),
                                    record_reference_type='gi'
                                    )

                                seq_rep_list = DB.produce_seq_edits(
                                    s1=str(seq).upper(),
                                    s2=str(ar.seq).upper())

                                # print(seq[0:100])
                                # print(str(ar.seq).upper()[0:100])
                                # print(seq_rep)

                                rec_id = DB.db_get_row_ids(
                                    table_name='records',
                                    where_dict={'ncbi_gi': int(ar.id)})[0]

                                seq_id = DB.db_get_row_ids(
                                    table_name='sequences',
                                    where_dict={'rec_id': rec_id})[0]

                                seq_rep_id = DB.add_sequence_representation(
                                    seq_id=seq_id,
                                    repr_list=seq_rep_list,
                                    rec_id=None,
                                    aln_id=None)[0]

                                seq_rep_id_list.append(seq_rep_id)

                            consensus = kralign.consensus(
                                alignment=aln,
                                threshold=0.4,
                                unknown='N',
                                resolve_ambiguities=FLAT_RESOLVE_AMBIGUITIES)

                            new_seq_str = str(consensus[0])

                            # print(consensus[0])

                            # print('--- --- --- --- --- --- --- ---')

                    else:

                        trimmed_rec = wf.trim_record_to_locus(
                            record=org_records[0],
                            locus_name=locus_name)

                        if trimmed_rec:

                            flat_locus_produced = True

                            new_seq_str = str(trimmed_rec.seq)
                            krbioio.write_sequence_file(
                                records=[trimmed_rec],
                                file_path=seq_file_path,
                                file_format='fasta')

                    ############################################################

                    if flat_locus_produced:

                        internal_reference = locus_name + '_' + krother.random_id(6)
                        parent_rec_id_list = [int(x.annotations['kr_seq_db_id']) for x in org_records]

                        rec_id = DB.add_record(
                            org_id=org_id,
                            ncbi_gi=None,
                            ncbi_version=None,
                            internal_reference=internal_reference,
                            description=None,
                            sequence_str=new_seq_str,
                            sequence_alphabet_str=sequence_alphabet_str,
                            parent_rec_id_list=parent_rec_id_list,
                            active=True,
                            action_str='Record for locus')

                        DB.add_record_annotation(
                            record_reference=rec_id,
                            type_str='locus_flat',
                            annotation_str=locus_name,
                            record_reference_type='raw')

                        if len(seq_rep_id_list) > 1:
                            DB.add_alignment(
                                name=aln_name,
                                seq_rep_id_list=seq_rep_id_list,
                                description=None,
                                rec_id=rec_id)

                # print('--- --- --- --- --- --- --- --- --- --- --- --- --- ---')

            ####################################################################

            DB.save()

            print()

    ############################################################################

    DB.close()

    ############################################################################
