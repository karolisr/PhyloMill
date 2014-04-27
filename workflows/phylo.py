#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
# from __future__ import unicode_literals

def write_log(msg, log_file_path, append=True, newlines=0):
    from krpy import krother
    log_handle = None
    if append:
        log_handle = open(log_file_path, 'a')
    else:
        log_handle = open(log_file_path, 'wb')
    timestamp = krother.timestamp()
    msg = newlines*'\n' + timestamp + ' - ' + msg
    log_handle.write(msg)
    log_handle.write('\n')
    log_handle.close()
    print(msg)
    return

if __name__ == '__main__':

    import os
    import sys
    import argparse
    import inspect
    import shutil
    import ConfigParser

    from krpy import KRSequenceDatabase
    from krpy.workflow_functions import krphylowf as wf
    from krpy import krncbi
    from krpy import krbioio

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
        print('Project directory is required.')
        sys.exit(0)

    # Log file path
    LFP = PRJ_DIR_PATH + 'log.txt'

    if ARGS.project_dir:

        if os.path.exists(PRJ_DIR_PATH):

            msg = 'Using project directory at ' + PRJ_DIR_PATH.rstrip(PS)
            write_log(msg=msg, log_file_path=LFP, append=True, newlines=1)
            DB = KRSequenceDatabase.KRSequenceDatabase(
                PRJ_DIR_PATH + 'db.sqlite3')

        else:

            # Script filename
            script_file_path = inspect.getfile(inspect.currentframe())

            # Script directory path
            script_dir_path = os.path.dirname(os.path.abspath(script_file_path))

            prj_template_dir_path = script_dir_path.strip('workflows') + \
                'data' + PS + 'phylo-prj-tempate'

            shutil.copytree(prj_template_dir_path, PRJ_DIR_PATH, symlinks=False,
                            ignore=None)

            # wd = os.getcwd()
            # os.chdir(PRJ_DIR_PATH+'03-name_resolution')
            # call(['./get_ncbi_data.sh'], stdout=open(os.devnull, 'wb'))
            # os.chdir(wd)

            KRSequenceDatabase.KRSequenceDatabase(PRJ_DIR_PATH + 'db.sqlite3')

            print('')
            msg = 'Creating project directory at ' + PRJ_DIR_PATH.rstrip(PS)
            write_log(msg=msg, log_file_path=LFP, append=False)

            sys.exit(0)

    ############################################################################

    TEMP_DIR_PATH = PRJ_DIR_PATH + 'temporary' + PS
    OUT_DIR_PATH = PRJ_DIR_PATH + 'output' + PS

    ############################################################################

    if not COMMANDS:
        write_log('No commands given.', LFP)
        sys.exit(0)

    ############################################################################

    # Get configuration information
    config_file_path = PRJ_DIR_PATH + 'config'
    config = ConfigParser.SafeConfigParser(allow_no_value=True)
    config.optionxform=str
    config.read(config_file_path)

    email = config.get('General', 'email')
    max_seq_length = config.getint('General', 'max_seq_length')
    seq_type = config.get('General', 'seq_type')

    # Parse taxa
    tax_temp = config.items('Taxa')
    tax_temp = [x[0] for x in tax_temp]
    tax_ids = list()
    for tax in tax_temp:
        if tax.isalpha():
            tax_id = list(krncbi.esearch(tax, 'taxonomy', email))[0]
            tax_ids.append(str(tax_id))
        else:
            tax_ids.append(tax)

    # Parse loci
    loci_temp = config.items('Loci')
    loci = dict()
    for l in loci_temp:
        loci[l[0]] = [x.strip() for x in l[1].split(',')]

    ############################################################################

    # Search genbank
    if 'search' in COMMANDS:

        ncbi_db = None

        if seq_type == 'aa':
            ncbi_db = 'protein'
        elif seq_type == 'nt':
            ncbi_db = 'nuccore'

        for locus_name in loci.keys():

            query_terms = loci[locus_name]

            msg = 'Searching NCBI ' + ncbi_db + ' database for ' + \
                  locus_name + '.'
            write_log(msg, LFP, newlines=1)

            query_terms = loci[locus_name]

            gis = wf.search_genbank(
                ncbi_db=ncbi_db,
                query_terms=query_terms,
                ncbi_tax_ids=tax_ids,
                max_seq_length=max_seq_length,
                email=email)

            msg = 'Found ' + str(len(gis)) + ' records.'
            write_log(msg, LFP)

            new_gis = list()
            in_db_gis = list()

            for gi in gis:

                in_blacklist = DB.in_blacklist(
                    record_reference = gi,
                    record_reference_type='gi')

                in_db = DB.in_db(
                    record_reference = gi,
                    record_reference_type='gi')

                if in_blacklist:
                    pass
                elif in_db:
                    in_db_gis.append(gi)
                else:
                    new_gis.append(gi)

            msg = 'There are ' + str(len(new_gis)) + ' new records.'
            write_log(msg, LFP)

            gb_file_path = TEMP_DIR_PATH + locus_name + '.gb'

            if len(new_gis) > 0:

                msg = 'Downloading.'
                write_log(msg, LFP)

                print('')

                krncbi.download_sequence_records(
                    file_path=gb_file_path,
                    uids=new_gis,
                    db=ncbi_db,
                    entrez_email=email)

                print('')

                msg = 'Filtering records.'
                write_log(msg, LFP)

                records_new = krbioio.read_sequence_file(
                    file_path=gb_file_path,
                    file_format='gb',
                    ret_type='list',
                    key='gi')

                records_in_db = DB.get_records(
                    record_reference_list=in_db_gis,
                    record_reference_type='gi'
                    )

                # for r in records_in_db:
                #     print(r)
                #     print('-------------------------------------------------')

                records = records_new + records_in_db

                good_bad = wf.accept_records(
                    records=records,
                    temp_dir=TEMP_DIR_PATH,
                    min_clust_size=10)

                good = good_bad[0]
                bad = good_bad[1]

                good_new = list()
                bad_new = list()

                for g in good:
                    gi = int(g[1])
                    if gi in new_gis:
                        good_new.append(gi)

                for b in bad:
                    gi = int(b[1])
                    if gi in new_gis:
                        bad_new.append(gi)

                records_to_add = list()

                for record in records_new:
                    gi = int(record.annotations['gi'])
                    if gi in good_new:
                        records_to_add.append(record)


                msg = 'Accepted ' + str(len(records_to_add)) + ' out of '  + \
                      str(len(records_new)) + ' records.'
                write_log(msg, LFP)

                if len(records_to_add) > 0:

                    msg = 'Adding downloaded records to database.'
                    write_log(msg, LFP)

                    for record in records_to_add:
                        DB.add_genbank_record(
                            record=record,
                            action_str='initial search')

                        gi = int(record.annotations['gi'])

                        DB.add_record_annotation(
                            record_reference=gi,
                            type_str='locus',
                            annotation_str=locus_name,
                            record_reference_type='gi')

        DB.save()

    ############################################################################

    DB.close()

    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
