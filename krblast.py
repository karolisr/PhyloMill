from __future__ import print_function
# from __future__ import unicode_literals


def search(program, db, query_fasta_file, output_directory, evalue, threads, output_prefix):
    import subprocess
    import os
    output_file = output_directory + os.path.sep + output_prefix + program + '_' + db + '.xml'
    print('Running ' + program + ' on ' + db + ' database using query file "'+os.path.split(query_fasta_file)[1]+'" with E-Value of ' + evalue + '.')
    subprocess.call(program+' -db '+db+' -query '+query_fasta_file+' -out '+output_file+' -evalue '+str(evalue)+' -outfmt 5 -num_threads '+str(threads), shell=True)
    return output_file

def gis_from_blast_results(file_name, min_sequence_length=-1, max_sequence_length=-1):

    from Bio.Blast import NCBIXML

    results_handle = open(file_name)
    blast_records = NCBIXML.parse(results_handle)

    gis = []

    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            if (min_sequence_length == -1 or alignment.length >= min_sequence_length) and \
               (max_sequence_length == -1 or alignment.length <= max_sequence_length):

                gi = alignment.title.split('|')[1]
                gis.append(str(gi))

    results_handle.close()
    print("Found " + str(len(gis)) + " gis.")
    gis = list(set(gis))
    print("There are " + str(len(gis)) + " unique gis.")
    return gis

if __name__ == '__main__':

    # Tests

    import os

    PS = os.path.sep

    pass
