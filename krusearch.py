from __future__ import print_function
from __future__ import unicode_literals

'''
    This module is designed for USEARCH 6. There are many differences between
    USEARCH versions 6 and 5. This module is incompatible with version 5.
'''

def cluster_file(

    input_file_path,
    output_file_path,
    identity_threshold,
    sorted_input=False,
    algorithm='fast', # fast smallmem
    strand='plus', # plus both
    threads=1,
    quiet=True,
    program='usearch6'

    ):

    import os
    import subprocess

    if quiet:
        quiet = ' -quiet'
    else:
        quiet = ''

    if algorithm == 'smallmem' and not sorted_input:
        subprocess.call(
            (program + quiet +
             ' -sortbylength ' + input_file_path +
             ' -output ' + input_file_path + '_sorted'
             ),
            shell=True)
        input_file_path = input_file_path + '_sorted'

    if algorithm == 'fast':
        threads = ' -threads ' + str(threads)
    else:
        threads = ''

    subprocess.call(
        (program + quiet +
         ' -cluster_' + algorithm + ' ' + input_file_path +
         ' -strand ' + strand +
         ' -id ' + str(identity_threshold) +
         threads +
         ' -uc ' + output_file_path
         ),
        shell=True)

    if algorithm == 'smallmem' and not sorted_input:
        os.remove(input_file_path)

    return output_file_path

def parse_uc_file(uc_file_path):

    import csv
    uc = csv.reader(open(uc_file_path, 'rb'), delimiter='\t')

    cluster_dict = {}

    for row in uc:
        if row[0].startswith('S'):
            cluster_dict[row[1]] = [row[8]]
        elif row[0].startswith('H'):
            cluster_dict[row[1]].append(row[8])

    return cluster_dict

if __name__ == '__main__':
    
    # Tests
    
    import os

    PS = os.path.sep

    to_cluster_file='testdata/tocluster.fasta'
    output_file='testdata/clustered.uc'
    cluster_file(to_cluster_file,
                 output_file,
                 identity_threshold=0.99,
                 sorted_input=False,
                 algorithm='smallmem',
                 strand='both',
                 threads=4,
                 quiet=False,
                 program='usearch6')
