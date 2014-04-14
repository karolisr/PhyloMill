from __future__ import print_function
# from __future__ import unicode_literals


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
    if isinstance(records, dict):
        records = records.values()
    handle = open(file_path, 'w')
    count_written = SeqIO.write(records, handle, file_format)
    handle.close()
    return count_written

def export_records(records, file_format, file_path, seq_id=None):
    import copy
    if file_format == 'fasta':
        records_copy = []
        for rec in records:
            records_copy.append(copy.copy(rec))
        for record in records_copy:
            if seq_id is None or seq_id == 'gi':
                record.id = record.annotations['gi']
            if seq_id == 'accession':
                pass
            record.description = ''
        write_sequence_file(records_copy, file_path, file_format)
    elif file_format == 'genbank':
        write_sequence_file(records_copy, file_path, file_format)

# ToDo: File formats should be listed for convenience.
# https://github.com/biopython/biopython/blob/master/Bio/SeqIO/QualityIO.py


def read_alignment_file(file_path, file_format):

    from Bio import AlignIO
    handle = open(file_path, 'rU')
    alignment = None
    try:
        alignment = AlignIO.read(handle, file_format)
    except:
        pass
    handle.close()
    return alignment


def write_alignment_file(alignment, file_path, file_format):

    from Bio import AlignIO
    handle = open(file_path, 'w')
    count_written = None
    if alignment:
        count_written = AlignIO.write(alignment, handle, file_format)
    handle.close()
    return count_written


def split_fastq_file(pieces, output_dir, forward_reads_file_path,
                     reverse_reads_file_path=None, log_func=None,
                     log_file_path=None):

    import os
    import krio

    msg = 'Splitting FASTQ file into ' + str(pieces) + ' pieces.'
    print(msg)
    if log_func and log_file_path:
        log_func(msg, log_file_path)

    krio.prepare_directory(output_dir)
    print('Counting reads, this may take some time...')
    num_lines = krio.num_lines_in_file(forward_reads_file_path, print_every=400000)
    msg = 'There are ' + str(num_lines / 4) + ' records.'
    print(msg)
    if log_func and log_file_path:
        log_func(msg, log_file_path)
    records_per_file = num_lines / 4 / pieces

    forward_file_handles = list()
    reverse_file_handles = list()

    for piece in range(0, pieces):
        handle = open(output_dir + os.path.sep + 'f_' + str(piece + 1) +
                      '.fastq', 'wa')
        forward_file_handles.append(handle)
        if reverse_reads_file_path:
            handle = open(output_dir + os.path.sep + 'r_' + str(piece + 1) +
                          '.fastq', 'wa')
            reverse_file_handles.append(handle)

    forward_file_handles.reverse()
    reverse_file_handles.reverse()

    msg = '\nSplitting forward reads.\n'
    print(msg)
    if log_func and log_file_path:
        log_func(msg, log_file_path)
    with open(forward_reads_file_path) as f:
        write_handle = None
        lines_written = 0
        for i, l in enumerate(f):
            if (len(forward_file_handles) and
                    ((float(i) / 4) % records_per_file == 0)):
                if lines_written != 0:
                    msg = '\tWritten ' + str(lines_written / 4) + ' records.'
                    print(msg)
                    if log_func and log_file_path:
                        log_func(msg, log_file_path)
                    lines_written = 0
                msg = ('\t' + str(len(forward_file_handles)) +
                       ' files remaining.')
                print(msg)
                if log_func and log_file_path:
                    log_func(msg, log_file_path)
                write_handle = forward_file_handles.pop()
            write_handle.write(l)
            lines_written = lines_written + 1
            if num_lines == i + 1:
                msg = '\tWritten ' + str(lines_written / 4) + ' records.'
                print(msg)
                if log_func and log_file_path:
                    log_func(msg, log_file_path)

    if reverse_reads_file_path:
        msg = '\nSplitting reverse reads.\n'
        print(msg)
        if log_func and log_file_path:
            log_func(msg, log_file_path)
        with open(reverse_reads_file_path) as f:
            write_handle = None
            lines_written = 0
            for i, l in enumerate(f):
                if (len(reverse_file_handles) and
                        ((float(i) / 4) % records_per_file == 0)):
                    if lines_written != 0:
                        msg = '\tWritten ' + str(lines_written / 4) + ' records.'
                        print(msg)
                        if log_func and log_file_path:
                            log_func(msg, log_file_path)
                        lines_written = 0
                    msg = ('\t' + str(len(reverse_file_handles)) +
                           ' files remaining.')
                    print(msg)
                    if log_func and log_file_path:
                        log_func(msg, log_file_path)
                    write_handle = reverse_file_handles.pop()
                write_handle.write(l)
                lines_written = lines_written + 1
                if num_lines == i + 1:
                    msg = '\tWritten ' + str(lines_written / 4) + ' records.'
                    print(msg)
                    if log_func and log_file_path:
                        log_func(msg, log_file_path)


if __name__ == '__main__':

    # Tests

    import os

    PS = os.path.sep

    pass
