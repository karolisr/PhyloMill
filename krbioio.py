from __future__ import print_function
#from __future__ import unicode_literals

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

def read_sequence_file(file_path, file_format, ret_type='list'):
    handle = open(file_path, 'rU')
    return read_sequence_data(handle, file_format, ret_type)

def write_sequence_file(records, file_path, file_format):
    from Bio import SeqIO
    from Bio import SeqRecord
    if isinstance(records, dict):
        records = records.values()
    handle = open(file_path, 'w')
    count_written = SeqIO.write(records, handle, file_format)
    handle.close()
    return count_written

# ToDo: File formats should be listed for convenience.
# https://github.com/biopython/biopython/blob/master/Bio/SeqIO/QualityIO.py

def read_alignment_file(file_path, file_format):

    from Bio import AlignIO
    handle = open(file_path, 'rU')
    alignment = AlignIO.read(handle, file_format)
    handle.close()
    return alignment

def write_alignment_file(alignment, file_path, file_format):

    from Bio import AlignIO
    handle = open(file_path, 'w')
    count_written = AlignIO.write(alignment, handle, file_format)
    handle.close()
    return count_written

if __name__ == '__main__':
    
    # Tests
    
    import os

    PS = os.path.sep

    pass
