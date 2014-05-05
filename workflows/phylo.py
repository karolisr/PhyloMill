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

    ARGS = PARSER.parse_args()

    ############################################################################

    COMMANDS = None
    if ARGS.commands:
        COMMANDS = set([x.strip() for x in ARGS.commands.split(',')])

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

    # Get configuration information
    CFG = ConfigParser.SafeConfigParser(allow_no_value=True)
    CFG.optionxform=str
    CFG.read(CFG_FILE_PATH)

    EMAIL = CFG.get('General', 'email')
    MAX_SEQ_LENGTH = CFG.getint('General', 'max_seq_length')

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
        else:
            TAX_IDS.append(tax)

    # Organism name resolution
    SYN = list()
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

            for l_sec in locus_cfg.sections():
                if l_sec == l:
                    l_db = locus_cfg.get(l, 'database')
                    l_query = locus_cfg.get(l, 'query')
                    l_dict['database'] = l_db
                    l_dict['query'] = l_query
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

            ncbi_db = LOCI[locus_name]['database']
            query_term_str = LOCI[locus_name]['query']

            msg = 'Searching NCBI ' + ncbi_db + ' database for ' + \
                  locus_name + '.'
            write_log(msg, LFP, newlines_before=1, newlines_after=0)

            gis = wf.search_genbank(
                ncbi_db=ncbi_db,
                query_term_str=query_term_str,
                ncbi_tax_ids=TAX_IDS,
                max_seq_length=MAX_SEQ_LENGTH,
                email=EMAIL)

            msg = 'Found ' + str(len(gis)) + ' records.'
            write_log(msg, LFP)

            gis_in_blacklist = list()
            gis_good = list()

            for gi in gis:

                in_blacklist = DB.in_blacklist(
                    record_reference = gi,
                    record_reference_type='gi')

                if in_blacklist:
                    gis_in_blacklist.append(gi)
                else:
                    gis_good.append(gi)

            msg = 'There are ' + str(len(gis_in_blacklist)) + \
                  ' blacklisted records.'
            write_log(msg, LFP)

            gis_in_db = list()
            gis_new = list()

            for gi in gis_good:

                in_db = DB.in_db(
                    record_reference = gi,
                    record_reference_type='gi')

                if in_db:
                    gis_in_db.append(gi)
                else:
                    gis_new.append(gi)

            msg = 'There are ' + str(len(gis_new)) + \
                  ' new records.'
            write_log(msg, LFP)

            if len(gis_new) > 0:

                msg = 'Downloading new records.'
                write_log(msg, LFP)

                print('')

                timestamp = krother.timestamp()
                timestamp = timestamp.replace('-', '_')
                timestamp = timestamp.replace(':', '_')
                timestamp = timestamp.replace(' ', '_')

                gb_file_name = locus_name + '_' + timestamp + '.gb'
                gb_file_path = DNLD_DIR_PATH + gb_file_name

                krncbi.download_sequence_records(
                    file_path=gb_file_path,
                    uids=gis_new,
                    db=ncbi_db,
                    entrez_email=EMAIL)

                print('')

                records_new = krbioio.read_sequence_file(
                    file_path=gb_file_path,
                    file_format='gb',
                    ret_type='list',
                    key='gi')

                records_to_add = records_new

                if len(records_to_add) > 0:

                    msg = 'Adding downloaded records to database.'
                    write_log(msg, LFP)

                    for record in records_to_add:
                        DB.add_genbank_record(
                            record=record,
                            action_str='Genbank search result.')

                        gi = int(record.annotations['gi'])

                        DB.add_record_annotation(
                            record_reference=gi,
                            type_str='locus',
                            annotation_str=locus_name,
                            record_reference_type='gi')

                        DB.add_record_annotation(
                            record_reference=gi,
                            type_str='source_file',
                            annotation_str=gb_file_name,
                            record_reference_type='gi')

            DB.save()

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

    DB.close()

    ############################################################################
