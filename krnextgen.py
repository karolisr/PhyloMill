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


def mask_low_quality_sites(bio_seq_record, quality_score_treshold,
                           low_quality_residue='N'):
    from Bio.SeqRecord import SeqRecord
    r = bio_seq_record
    quality_scores = r.letter_annotations['phred_quality']
    sequence = r.seq.tomutable()
    for index, score in enumerate(quality_scores):
        if score < quality_score_treshold:
            sequence[index] = low_quality_residue
    new_r = SeqRecord(seq=sequence.toseq(), id=r.id, name=r.name,
                      description=r.description, dbxrefs=r.dbxrefs,
                      features=r.features, annotations=r.annotations,
                      letter_annotations=r.letter_annotations)
    return(new_r)


def proportion_low_quality_sites(bio_seq_record, low_quality_residue='N'):
    sequence = str(bio_seq_record.seq).lower()
    low_quality_sites = sequence.count(low_quality_residue.lower())
    sequence_length = len(bio_seq_record.seq)
    prop_lq_sites = float(low_quality_sites) / float(sequence_length)
    return prop_lq_sites


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

    for barcode in barcodes:
        barcode['write_handle_forward'].close()
        if reverse_reads_file_path is not None:
            barcode['write_handle_reverse'].close()

    write_handle_forward_other.close()
    if reverse_reads_file_path is not None:
        write_handle_reverse_other.close()


def combine_demultiplexed_results(input_dir, output_dir):

    import os
    import shutil
    import krio
    import krpipe

    ps = os.path.sep
    input_dir = input_dir.rstrip(ps) + ps
    output_dir = output_dir.rstrip(ps) + ps
    krio.prepare_directory(output_dir)
    directory_list = krpipe.parse_directory(
        path=input_dir,
        file_name_sep=' ',
        sort='forward'
    )

    file_dir_path = directory_list[0]['path'].rstrip(ps) + ps
    file_list = krpipe.parse_directory(
        path=file_dir_path,
        file_name_sep='_',
        sort='forward'
    )

    for f in file_list:
        output_file_path = output_dir + f['full']
        output_file_handle = open(output_file_path, 'wa')
        for part in range(1, len(directory_list)+1):
            part_directory_path = input_dir + str(part) + ps
            part_file_path = part_directory_path + f['full']
            shutil.copyfileobj(open(part_file_path, 'rb'),
                               output_file_handle)
        output_file_handle.close()


# -----------------------------------------------------------------------------
# Functions that follow are used to estimate various statistics used in
# population genetics.
# -----------------------------------------------------------------------------


def like_homo(s, p, e):

    '''
        Calculates the likelihood (for a given site) of observed data
        conditional on the site being homozygous. Diploid individual is
        assumed. Lynch 2008.

        parameters:
            s - a list of counts of nucleotides (for a given site) in a cluster
                 0  1  2  3
                [A, C, G, T]

            p - average nucleotide frequencies in the region of analysis
                 0  1  2  3
                [A, C, G, T]

            e - error rate
    '''

    from scipy.stats import binom
    total_s = sum(s)
    likelihood = 0
    for i in range(0, 4):
        l = p[i] * binom.pmf(total_s-s[i], total_s, e)
        likelihood = likelihood + l
    return(likelihood)


def like_hetero(s, p, e):

    '''
        Calculates the likelihood (for a given site) of observed data
        conditional on the site being heterozygous. Diploid individual is
        assumed. Lynch 2008.

        parameters:
            s - a list of counts of nucleotides (for a given site) in a cluster
                 0  1  2  3
                [A, C, G, T]

            p - average nucleotide frequencies in the region of analysis
                 0  1  2  3
                [A, C, G, T]

            e - error rate
    '''

    from scipy.stats import binom
    total_s = sum(s)
    S = 1.0 - sum([x*x for x in p])
    likelihood = 0
    for i in range(0, 4):
        for j in range(i+1, 4):
            l = (2.0 * p[i] * p[j] *
                 binom.pmf(total_s-s[i]-s[j], total_s, (2*e)/3.0) *
                 binom.pmf(s[i], s[i]+s[j], 0.5) / S)
            likelihood = likelihood + l
    return(likelihood)


def like_homo_hetero(s, p, e, pi):

    '''
        Calculates the total likelihood (for a given site) of observed data.
        Diploid individual is assumed. Lynch 2008.

        parameters:
            s - a list of counts of nucleotides (for a given site) in a cluster
                 0  1  2  3
                [A, C, G, T]

            p - average nucleotide frequencies in the region of analysis
                 0  1  2  3
                [A, C, G, T]

            e - error rate

            pi - nucleotide diversity
    '''

    likelihood = (1-pi) * like_homo(s, p, e) + pi * like_hetero(s, p, e)
    return likelihood


def neg_ll_homo_hetero(ns, p, e, pi):

    '''
        Calculates the negative natural log of the total likelihood
        (for multiple sites) of observed data. Diploid individual is assumed.
        Lynch 2008.

        parameters:
            ns - a list of lists of counts of nucleotides (for a given site) in
            a cluster
                  0  1  2  3    0  1  2  3
                [[A, C, G, T], [A, C, G, T], ... ]
                    Site 1        Site 2     ...

            p - average nucleotide frequencies in the region of analysis
                 0  1  2  3
                [A, C, G, T]

            e - error rate

            pi - nucleotide diversity
    '''

    import numpy
    ll = 0
    for s in ns:
        l = like_homo_hetero(s, p, e, pi)
        ll = ll + numpy.log(l)
    ll = ll * (-1.0)
    return(ll)


def mle_e_and_pi(ns, p):

    '''
        Using neg_ll_homo_hetero, will produce a region-wide (could be whole
        genome) maximum likelihood estimate of e (error rate) and pi
        (nucleotide diversity).

        parameters:
            ns - a list of lists of counts of nucleotides (for a given site) in
            a cluster
                  0  1  2  3    0  1  2  3
                [[A, C, G, T], [A, C, G, T], ... ]
                    Site 1        Site 2     ...

            p - average nucleotide frequencies in the region of analysis
                 0  1  2  3
                [A, C, G, T]
    '''

    from scipy import optimize
    nll = lambda estimated, ns_l, p_l: (
        neg_ll_homo_hetero(ns_l, p_l, estimated[0], estimated[1])
    )
    return(optimize.fmin_l_bfgs_b(
        nll,
        x0=(0.0001, 0.0001),
        args=(ns, p),
        bounds=((1E-10, 0.99999), (1E-10, 0.99999)),
        approx_grad=True
    ))


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

    # combine_demultiplexed_results
    # combine_demultiplexed_results(
    #     input_dir='/home/karolis/Dropbox/code/test/rad/02-demultiplexed-fastq-parts',
    #     output_dir='/home/karolis/Dropbox/code/test/rad/03-demultiplexed-fastq-combined')

    # ns = [[99, 1, 0, 0],
    #      [50, 49, 1, 0],
    #      [98, 1, 1, 0],
    #      [50, 50, 0, 0],
    #      [100, 0, 0, 0],
    #      [50, 50, 0, 0],
    #      [100, 0, 0, 0],
    #      [100, 0, 0, 0],
    #      [50, 0, 0, 50],
    #      [98, 0, 2, 0],
    #      [48, 1, 50, 1],
    #      [0, 50, 50, 0]]

    # ns = [[100, 0, 0, 0],
    #      [100, 0, 0, 0],
    #      [100, 0, 0, 0],
    #      [100, 0, 0, 0],
    #      [100, 0, 0, 0],
    #      [100, 0, 0, 0]]

    # ns = [[100, 0, 100, 0],
    #      [100, 0, 0, 0],
    #      [100, 0, 100, 0],
    #      [100, 0, 0, 0],
    #      [100, 0, 100, 0],
    #      [100, 0, 0, 0]]

    # ns = [[100, 0, 100, 0],
    #      [100, 0, 100, 0],
    #      [100, 0, 100, 0],
    #      [100, 0, 100, 0],
    #      [100, 0, 100, 0],
    #      [100, 0, 100, 0]]

    # p = [0.25, 0.25, 0.25, 0.25]
    # e = 1E-10
    # pi = 0.5

    # import numpy

    # l1 = like_homo(ns[1], p, e)
    # print('Homozygote')
    # print(l1)
    # print(numpy.log(l1))
    # print('--- --- --- --- --- --- --- --- --- --- ---')

    # l2 = like_hetero(ns[1], p, e)
    # print('Heterozygote')
    # print(l2)
    # print(numpy.log(l2))
    # print('--- --- --- --- --- --- --- --- --- --- ---')

    # l3 = like_homo_hetero(ns[1], p, e, pi)
    # print('Total')
    # print(l3)
    # print(numpy.log(l3))
    # print('--- --- --- --- --- --- --- --- --- --- ---')

    # l4 = neg_ll_homo_hetero(ns, p, e, pi)
    # print('Negative LN total')
    # print(l4)
    # print('--- --- --- --- --- --- --- --- --- --- ---')

    # ml_est = mle_e_and_pi(ns, p)
    # print('                      e :', ml_est[0][0])
    # print('                     pi :', ml_est[0][1])
    # print('negative log likelihood :', ml_est[1])
