from __future__ import print_function
# from __future__ import unicode_literals


def search(program, db, query_fasta_file, output_directory, evalue, threads, output_prefix):
    import subprocess
    import os
    output_file = output_directory + os.path.sep + output_prefix + program + '_' + db + '.xml'
    print('Running ' + program + ' on ' + db + ' database using query file "'+os.path.split(query_fasta_file)[1]+'" with E-Value of ' + evalue + '.')
    subprocess.call(program+' -db '+db+' -query '+query_fasta_file+' -out '+output_file+' -evalue '+str(evalue)+' -outfmt 5 -num_threads '+str(threads), shell=True)
    return(output_file)


if __name__ == '__main__':

    # Tests

    import os

    PS = os.path.sep

    pass
