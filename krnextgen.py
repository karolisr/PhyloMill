from __future__ import print_function
#from __future__ import unicode_literals


def read_barcodes(file_path, delimiter, id_header, barcode_header):

    '''
    Read barcodes from a CSV table and return a list of dictionaries with keys
    id and barcode.
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


def _write_demultiplex_results_(barcodes,
                                reverse_reads_file_path,
                                result_batch_forward_other,
                                write_handle_forward_other,
                                result_batch_reverse_other,
                                write_handle_reverse_other,
                                output_file_format):
    from Bio import SeqIO

    for barcode in barcodes:
        SeqIO.write(barcode['result_batch_forward'],
                    barcode['write_handle_forward'], output_file_format)
        del barcode['result_batch_forward'][:]
        if reverse_reads_file_path is not None:
            SeqIO.write(barcode['result_batch_reverse'],
                        barcode['write_handle_reverse'], output_file_format)
            del barcode['result_batch_reverse'][:]

    SeqIO.write(result_batch_forward_other,
                write_handle_forward_other, output_file_format)
    del result_batch_forward_other[:]

    if reverse_reads_file_path is not None:
        SeqIO.write(result_batch_reverse_other,
                    write_handle_reverse_other, output_file_format)
        del result_batch_reverse_other[:]


def demultiplex(barcodes,
                forward_reads_file_path,
                reverse_reads_file_path=None,
                input_file_format='fastq',
                output_file_format='fastq',
                max_barcode_mismatch_count=1,
                output_dir='.',
                trim_barcode=True,
                trim_extra=0,
                write_every=1000
                ):

    import os
    import Levenshtein
    from Bio import SeqIO
    from Bio import SeqRecord
    import krio
    import krseq

    ps = os.path.sep
    output_dir = output_dir.rstrip(ps) + ps
    krio.prepare_directory(output_dir)

    forward_reads = SeqIO.parse(forward_reads_file_path, input_file_format)
    reverse_reads = None
    if reverse_reads_file_path is not None:
        reverse_reads = SeqIO.parse(reverse_reads_file_path, input_file_format)

    for barcode in barcodes:

        barcode['length'] = len(barcode['barcode'])
        base_file_name = output_dir + barcode['id'] + '_' + barcode['barcode']

        barcode['result_batch_forward'] = list()
        barcode['file_path_forward'] = (base_file_name + '_f.' +
                                        output_file_format)
        barcode['write_handle_forward'] = open(barcode['file_path_forward'],
                                               'wa')

        if reverse_reads_file_path is not None:

            barcode['result_batch_reverse'] = list()
            barcode['file_path_reverse'] = (base_file_name + '_r.' +
                                            output_file_format)
            barcode['write_handle_reverse'] = open(
                barcode['file_path_reverse'], 'wa')

    write_handle_forward_other = open(output_dir + 'mismatch_f.' +
                                      output_file_format, 'wa')
    result_batch_forward_other = list()
    write_handle_reverse_other = None
    result_batch_reverse_other = None
    if reverse_reads_file_path is not None:
        write_handle_reverse_other = open(output_dir + 'mismatch_r.' +
                                          output_file_format, 'wa')
        result_batch_reverse_other = list()

    # Loop over all reads
    for i, f_record in enumerate(forward_reads):
        # print(i+1, end='\r')
        f_record = SeqRecord.SeqRecord(
            f_record.seq,
            f_record.description.replace(' ', '|'),
            '',
            '',
            f_record.dbxrefs,
            f_record.features,
            f_record.annotations,
            f_record.letter_annotations)

        # Try all barcodes
        barcode_match_found = False
        for barcode in barcodes:

            l = barcode['length']
            b = barcode['barcode'].lower()
            r = str(f_record.seq[0:l]).lower()
            ld = Levenshtein.distance(b, r)

            # Barcode match found
            if ld <= max_barcode_mismatch_count:
                barcode_match_found = True
                if trim_barcode:
                    f_record = krseq.trim_residues(f_record, l, False)
                    if trim_extra > 0:
                        f_record = krseq.trim_residues(f_record, trim_extra,
                                                       False)
                barcode['result_batch_forward'].append(f_record)

                # If there are reverse reads to consider
                if reverse_reads_file_path is not None:
                    r_record = reverse_reads.next()
                    r_record = krseq.reverse_complement(r_record)
                    r_record = SeqRecord.SeqRecord(
                        r_record.seq,
                        r_record.description.replace(' ', '|'),
                        '',
                        '',
                        r_record.dbxrefs,
                        r_record.features,
                        r_record.annotations,
                        r_record.letter_annotations)

                    if trim_barcode and trim_extra > 0:
                        r_record = krseq.trim_residues(r_record, trim_extra,
                                                       True)

                    barcode['result_batch_reverse'].append(r_record)

                # Barcode match found, breakout of the loop
                break

        if not barcode_match_found:
            result_batch_forward_other.append(f_record)
            if reverse_reads_file_path is not None:
                r_record = reverse_reads.next()
                result_batch_reverse_other.append(r_record)

        if i % write_every == write_every - 1 and i > 0:
            _write_demultiplex_results_(barcodes,
                                        reverse_reads_file_path,
                                        result_batch_forward_other,
                                        write_handle_forward_other,
                                        result_batch_reverse_other,
                                        write_handle_reverse_other,
                                        output_file_format)

    _write_demultiplex_results_(barcodes,
                                reverse_reads_file_path,
                                result_batch_forward_other,
                                write_handle_forward_other,
                                result_batch_reverse_other,
                                write_handle_reverse_other,
                                output_file_format)

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

    # demultiplex
    # demultiplex(barcodes,
    #             forward_reads_file_path='testdata/rad_forward.fastq',
    #             reverse_reads_file_path='testdata/rad_reverse.fastq',
    #             input_file_format='fastq',
    #             output_file_format='fastq',
    #             max_barcode_mismatch_count=1,
    #             output_dir='../test/demultiplex',
    #             trim_barcode=True,
    #             trim_extra=5,
    #             write_every=1000
    #             )
