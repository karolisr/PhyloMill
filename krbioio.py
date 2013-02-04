from __future__ import print_function
from __future__ import unicode_literals

def read_sequence_data(handle, data_format, ret_type='list'):

    '''
    Read sequence data from Entrez.efetch or file handle and return a list of
    Biopython SeqRecord objects.
    '''

    from Bio import SeqIO
    records_generator = SeqIO.parse(handle, data_format)

    records = None

    if ret_type == 'dict':
        records = dict()
        for record in records_generator:
            records[record.id] = record
    else:
        records = list()
        for record in records_generator:
            records.append(record)

    handle.close()
    return records

if __name__ == '__main__':
    
    # Tests
    
    import os

    PS = os.path.sep

    pass
