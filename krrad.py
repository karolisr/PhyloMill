from __future__ import print_function
#from __future__ import unicode_literals


def read_barcodes(file_path, delimiter, id_header, barcode_header):

    '''
    Read RAD barcodes and return a list of dictionaries with keys id and
    barcode.
    '''

    import krio

    barcodes = krio.read_table_file(
        path=file_path, has_headers=True, headers=None, delimiter=delimiter,
        quotechar='"', stripchar='', rettype='dict')

    id_header_new = 'id'
    barcode_header_new = 'barcode'

    for r in barcodes:
        if id_header_new != id_header:
            r[id_header_new] = r[id_header]
            del r[id_header]
        if barcode_header_new != barcode_header:
            r[barcode_header_new] = r[barcode_header]
            del r[barcode_header]

    return(barcodes)


def split_rad_fastq_file(pieces, output_dir, forward_reads_file_path,
                         reverse_reads_file_path=None):

    import os
    import krio

    krio.prepare_directory(output_dir)
    num_lines = krio.num_lines_in_file(forward_reads_file_path)
    print('There are ' + str(num_lines / 4) + ' records.')
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

    print('Splitting forward reads.')
    with open(forward_reads_file_path) as f:
        write_handle = None
        lines_written = 0
        for i, l in enumerate(f):
            if (len(forward_file_handles) and
                    ((float(i) / 4) % records_per_file == 0)):
                if lines_written != 0:
                    print('\tWritten', str(lines_written / 4), 'records.')
                    lines_written = 0
                print('\t' + str(len(forward_file_handles)) +
                      ' files remaining.')
                write_handle = forward_file_handles.pop()
            write_handle.write(l)
            lines_written = lines_written + 1
            if num_lines == i + 1:
                print('\tWritten', str(lines_written / 4), 'records.')

    if reverse_reads_file_path:
        print('Splitting reverse reads...')
        with open(reverse_reads_file_path) as f:
            write_handle = None
            lines_written = 0
            for i, l in enumerate(f):
                if (len(reverse_file_handles) and
                        ((float(i) / 4) % records_per_file == 0)):
                    if lines_written != 0:
                        print('\tWritten', str(lines_written / 4), 'records.')
                        lines_written = 0
                    print('\t' + str(len(reverse_file_handles)) +
                          ' files remaining.')
                    write_handle = reverse_file_handles.pop()
                write_handle.write(l)
                lines_written = lines_written + 1
                if num_lines == i + 1:
                    print('\tWritten', str(lines_written / 4), 'records.')

if __name__ == '__main__':
    # Tests
    import os
    ps = os.path.sep

    # read_barcodes
    # barcodes = read_barcodes(
    #     file_path='testdata' + ps + 'rad_barcodes.tsv',
    #     delimiter='\t',
    #     id_header='id',
    #     barcode_header='barcode')
    # for r in barcodes:
    #     print(r)

    # split_rad_fastq_file
    # split_rad_fastq_file(
    #     pieces=4,
    #     output_dir='../test/fastq-split',
    #     forward_reads_file_path='testdata/rad_forward.fastq',
    #     reverse_reads_file_path='testdata/rad_reverse.fastq')
