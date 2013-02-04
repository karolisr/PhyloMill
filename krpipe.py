#!/usr/bin/env python

from __future__ import print_function
from __future__ import unicode_literals

# Utility functions used by the pipeline functions.

def parse_directory(path, file_name_sep):

    '''
    Will parse a directory at a given path and return a list of dictionary
    objects with keys:
        name: name of a file without extension
        ext: file extension
        full: file name with extension
        path: full file path, relative to the input path
        split: file name split using file_name_sep input variable
    '''
    
    import os
    
    ps = os.path.sep
    
    file_list = os.listdir(path)
    
    return_list = list()
    
    for f in file_list:
    
        file_name = os.path.splitext(f)[0]
        file_ext = os.path.splitext(f)[1].split('.')[1]
        file_name_split = file_name.split(file_name_sep)
    
        file_dict = dict()
    
        file_dict['name'] = file_name
        file_dict['ext'] = file_ext
        file_dict['full'] = f
        file_dict['path'] = path + ps + f
        file_dict['split'] = file_name_split
    
        return_list.append(file_dict)
    
    return return_list

# Pipeline functions.

def search_and_download(queries, output_dir, file_name_sep, email):

    '''
    This will search NCBI and download sequences.
    '''

    import krio
    import krncbi

    print('\nSearching NCBI.')
    print('\tPreparing output directory "', output_dir, '"', sep='')

    krio.prepare_directory(output_dir)

    for query_dict in queries:

        name1 = query_dict['name1']
        name2 = query_dict['name2']
        locus = query_dict['locus']
        minlen = query_dict['minlen']
        feature_type = query_dict['ncbi_feature_type']
        qualifier_label = query_dict['ncbi_qualifier_label']
        db = query_dict['database']
        query = query_dict['query']

        print('\n\tSearching for:', name1, name2, 'Locus:', locus, 'in', db,
            'database.\n', sep=' ')

        '''
        File name is genrated based on the search query (periods will be
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

        file_name = (output_dir +
                     name1 + file_name_sep +
                     name2 + file_name_sep +
                     locus + file_name_sep +
                     minlen + file_name_sep +
                     feature_type + file_name_sep +
                     qualifier_label + file_name_sep +
                     db + '.gb')
        # Search NCBI.
        result_uids = krncbi.esearch(query, db, email)
        # Download records.
        krncbi.download_sequence_records(file_name, result_uids, db, email)

# End pipeline functions.

if __name__ == '__main__':
    
    import sys
    import argparse
    import os
    import krio
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--test', action='store_true', help='Run tests.')
    parser.add_argument('-c', '--command', type=unicode,
        choices=['search_and_download'], help='Run a command.')
    parser.add_argument('-q', '--query', type=unicode, help='Query file path.')
    parser.add_argument('-o', '--output', type=unicode,
        help='Output directory path.')
    parser.add_argument('-s', '--sep', type=unicode,
        help='Output file name separator.')
    parser.add_argument('-e', '--email', type=unicode, help='Email.')
    
    args = parser.parse_args()

    PS = os.path.sep

    if args.test:

        print('Running tests.')

        # Tests

        # parse_directory
        t_parse_directory = parse_directory('testdata'+PS+'parsedirectory',
            '$')
        for d in t_parse_directory:
            print(d)

        sys.exit(0)

    else:

        if args.command == 'search_and_download':

            queries = krio.read_table_file(args.query, has_headers=True,
                          headers=None, delimiter=b',', iterator=False)

            search_and_download(queries, args.output+PS, args.sep, args.email)
