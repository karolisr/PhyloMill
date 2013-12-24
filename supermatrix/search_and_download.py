#!/usr/bin/env python

'''
    DOCSTRING
'''

from __future__ import print_function
from __future__ import unicode_literals

def search_and_download(
    queries,
    output_dir,
    file_name_sep,
    email,
    log_dir):

    '''
    This will search NCBI and download sequences.
    '''

    import os
    from krpy import krio
    from krpy import krncbi
    from krpy import krbioio

    ps = os.path.sep

    print('\nSearching NCBI.')
    print('\n\tPreparing output directory "', output_dir, '"', sep='')

    krio.prepare_directory(output_dir)
    krio.prepare_directory(log_dir)

    all_uids = set()

    for query_dict in queries:

        name1 = query_dict['name1']
        name2 = query_dict['name2']
        db = query_dict['database']
        query = query_dict['query']

        '''
        File name is generated based on the search query (periods will be
        replaced with file_name_sep):

            name1.name2.locus.minlen.feature_type.qualifier_label
            .database.[other].extension

        Short explanation:

            Full name of a query consists of a name1 and a name2. This
            is to allow for multiple query strings for the same locus.
            An example is our good friend trnK/matK. trnK contains matK,
            so we want to be sure we get all trnK and matK loci because it
            may happen that matK will not be annotated within trnK. While
            it is possible to introduce "OR" statements into NCBI queries
            and collapse the search query into one string, it does not help
            us later if we wish to treat the results independently.
        '''

        file_path = (output_dir.rstrip(ps) + ps +
                     name1 + file_name_sep +
                     name2 + file_name_sep +
                     '.gb')

        if query.lower().startswith('donotdownload'):
            all_records = dict()
            file_list = krio.parse_directory(output_dir, file_name_sep)
            for f in file_list:
                if not f['ext'].startswith('gb'):
                    continue
                if f['split'][0] == name1:
                    records = krbioio.read_sequence_file(
                        file_path=f['path'],
                        file_format='genbank',
                        ret_type='dict')
                    for key in records.keys():
                        if records[key] not in all_records.keys():
                            all_records[key] = records[key]

            all_records = all_records.values()
            print('\n\tWill not download:', name1, name2, '\n', sep=' ')
            krbioio.write_sequence_file(records, file_path, 'genbank')

        else:

            print('\n\tSearching for:', name1, name2, 'in', db, 'database.\n',
                  sep=' ')

            # Search NCBI.
            result_uids = krncbi.esearch(query, db, email)
            all_uids |= set(result_uids)
            # Download records.
            krncbi.download_sequence_records(file_path, result_uids, db, email)

    # downloaded_gis_file = log_dir.rstrip(ps) + ps + '01-downloaded-gis.csv'
    # novel_gis_file = log_dir.rstrip(ps) + ps + '01-novel-gis.csv'
    # missing_gis_file = log_dir.rstrip(ps) + ps + '01-missing-gis.csv'

    # # First find previously downloaded GIs and see how many of them are the same
    # previous_uids = krio.read_table_file(
    #     path=downloaded_gis_file,
    #     has_headers=False,
    #     headers=None,
    #     delimiter=',',
    #     quotechar=None,
    #     stripchar='"',
    #     rettype='set')

    # novel_uids = all_uids - previous_uids
    # missing_uids = previous_uids - all_uids

    # novel_uids = list(novel_uids)
    # handle = open(novel_gis_file, 'w')
    # for uid in novel_uids:
    #     handle.write(uid+'\n')
    # handle.close()

    # missing_uids = list(missing_uids)
    # handle = open(missing_gis_file, 'w')
    # for uid in missing_uids:
    #     handle.write(uid+'\n')
    # handle.close()

    # all_uids = list(all_uids)
    # handle = open(downloaded_gis_file, 'w')
    # for uid in all_uids:
    #     handle.write(uid+'\n')
    # handle.close()


# if __name__ == '__main__':

#     from krpy import krio

#     QUERY_FILE = '../testdata/search_and_align/search_queries'
#     QUERIES = krio.read_table_file(QUERY_FILE, has_headers=True,
#                                    headers=None, delimiter='\t',
#                                    quotechar="'")
#     OUTPUT_DIR = '../testdata/testoutputdir'
#     LOG_DIR = '../testdata/testlogdir'
#     FILE_NAME_SEP = '$'
#     EMAIL = 'test@test.com'

#     search_and_download(
#         QUERIES,
#         OUTPUT_DIR,
#         FILE_NAME_SEP,
#         EMAIL,
#         LOG_DIR)
