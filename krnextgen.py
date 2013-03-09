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
    sequence = str(bio_seq_record.seq)
    low_quality_sites = sequence.count(low_quality_residue)
    sequence_length = len(bio_seq_record.seq)
    prop_lq_sites = float(low_quality_sites) / float(sequence_length)
    return prop_lq_sites


def compare_sequences(s1, s2, ignore='N'):

    '''
        Compare two sequences of the same length. Returns a list of integers,
        one per site.

        Values:
            1 - match.
            0 - mismatch.
            2 - one or both of the nucleotides is an ignored character,
                N by default.

        Example output:
            1101111111112121
            ACGTACGTACGTNCGT
            ACCTACGTACGTNCNT
    '''

    mask = []
    for pos in zip(s1, s2):
        if ignore in pos:
            mask.append(2)
        elif pos[0] is pos[1]:
            mask.append(1)
        else:
            mask.append(0)
    return(mask)


def align_fr_reads(r1, r2, mmmr_cutoff=0.85, ignore='N'):

    '''
        Aligns forward and reverse reads. Returns position of the best
        alignment that has the most matching pairs of nucleotides and the
        match / (match + miss) is above mmmr_cutoff. Ignored
        character does not count as a miss.

        Returns:
            best_alignment = (match, total, ratio, (a, b), (c, d))

            [a:b] - is the alignment range on the first sequence
            [c:d] - is the alignment range on the second sequence
    '''

    l1 = len(r1)
    l2 = len(r2)

    # r1 is expected to be shorter than r2, if this is not the case, we switch
    # r1 and r2
    switch = False
    if l1 > l2:
        switch = True
        tmp_str = r1
        r1 = r2
        r2 = tmp_str
        tmp_int = l1
        l1 = l2
        l2 = tmp_int

    # Slide the reads past each other and find the best alignment

    # Default setting
    best_alignment = (0, 0, 0, (None, None), (None, None))

    for i in range(1, l1 + l2):

        #print('i =', i)

        b = min(l1, i)
        c = max(0, l2 - i)
        d = l2
        if i > l1:
            d = l2 + l1 - i
        a = abs(d - b - c)

        #print('r1[', a, ':', b, ']', r1[a:b], sep='')
        #print('r2[', c, ':', d, ']', r2[c:d], sep='')

        mask = compare_sequences(r1[a:b], r2[c:d], ignore=ignore)

        match = mask.count(1)
        miss = mask.count(0)
        total = match + miss
        ratio = 0

        # Decide if this alignment is better than the current best alignment
        if float(total) > 0:
            ratio = float(match) / float(total)
        # Consider the alignment to be acceptable only if match/(match+miss)
        # meets given ratio
        if ratio >= mmmr_cutoff:
            # Furthermore, consider the alignment to be acceptable only if it
            # contains more matches than the current best alignment
            if match > best_alignment[0]:
                best_alignment = (match, total, ratio, (a, b), (c, d))

    # If the switch between r1 and r2 was made, we need to account for that
    # in our output
    if switch:
        best_alignment = (
            best_alignment[0],
            best_alignment[1],
            best_alignment[2],
            best_alignment[4],
            best_alignment[3])

    return(best_alignment)


def consensus_fr_read(r1, r2, min_overlap=5, mmmr_cutoff=0.85, ignore='N'):

    best_alignment = align_fr_reads(r1, r2, mmmr_cutoff=mmmr_cutoff,
                                    ignore=ignore)

    # best_alignment = (match, total, ratio, (a, b), (c, d))

    o = best_alignment[1]     # length of overlap excluding ignored characters
    r = best_alignment[2]     # match / (match + miss)

    a = best_alignment[3][0]  # [a:b] - alignment range on the first sequence
    b = best_alignment[3][1]

    c = best_alignment[4][0]  # [c:d] - alignment range on the second sequence
    d = best_alignment[4][1]

    str1_aln = r1
    str2_aln = r2

    cons = None

    # messages:
    #   0 - normal. forward read first, second read concatenated.
    #   1 - overlap. forward read first, second read overlaps part of the
    #       forward read.
    #   2 - complete overlap. reverse read first, reverse read overlaps the
    #       beginning of the forward read and may read into the
    #       adaptor + barcode.
    #   3 - overlap mismatch. forward and reverse reads do not agree on the
    #       identity of several bases.

    message = 0

    if r != 0 and r < 1.0:
        message = 3
    elif o >= min_overlap:
        cons = []
        if a < c:
            r1 = r1[a:b]
            str2_aln = r2[c:d]
            message = 2
        else:
            str1_aln = '{:-<{pad}}'.format(r1, pad=len(r2) - d - c + len(r1))
            str2_aln = '{:->{pad}}'.format(r2, pad=a + len(r2))
            message = 1
        for pos in zip(str1_aln, str2_aln):
            spos = set(pos)
            if '-' in spos:
                spos.remove('-')
            if len(spos) == 1:
                cons.append(spos.pop())
            else:
                if ignore in spos:
                    spos.remove(ignore)
                    cons.append(spos.pop())
                else:
                    cons.append(ignore)
        cons = ''.join([char for char in cons])
    else:
        cons = r1 + r2
        message = 0

    ret_value = [o, r, message, cons]

    return(ret_value)


def bin_reads(f_record, r_record=None, max_prop_low_quality_sites=0.10,
              min_overlap=5, mmmr_cutoff=0.85, low_quality_residue='N'):

    from Bio import Seq
    from Bio import SeqRecord

    f_hq = False
    r_hq = False

    f_lq_sites = proportion_low_quality_sites(
        bio_seq_record=f_record,
        low_quality_residue=low_quality_residue)
    if f_lq_sites <= max_prop_low_quality_sites:
        f_hq = True

    if r_record:
        r_lq_sites = proportion_low_quality_sites(
            bio_seq_record=r_record,
            low_quality_residue=low_quality_residue)
        if r_lq_sites <= max_prop_low_quality_sites:
            r_hq = True

    consensus = None
    if f_hq and r_hq:
        consensus = consensus_fr_read(
            r1=str(f_record.seq),
            r2=str(r_record.seq),
            min_overlap=min_overlap,
            mmmr_cutoff=mmmr_cutoff,
            ignore=low_quality_residue)

        cons_message = consensus[2]
        # If reads were aligned but disagreed on certain nucleotides, we treat
        # them as separate. Hopefully there were more reads at this locus and
        # those were aligned successfully. These separate reads should cluster
        # with the aligned reads later in the analysis and correct nucleotides
        # will be judged based on the whole cluster.
        if cons_message == 3:
            consensus = None
        else:
            cons_seq = Seq.Seq(consensus[3])
            cons_id = f_record.id.split('|')[0] + '|C0NS:' + str(cons_message)
            consensus = SeqRecord.SeqRecord(
                seq=cons_seq, id=cons_id, name='', description='')

    ret_value = [f_hq, r_hq, consensus]

    return(ret_value)


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


def nucleotides_at_site(site):
    c_a = site.count('A')
    c_c = site.count('C')
    c_g = site.count('G')
    c_t = site.count('T')
    return([c_a, c_c, c_g, c_t])


def align_clusters(min_seq_cluster, max_seq_cluster, uc_file_path,
                   fasta_file_path, aln_output_file_path,
                   counts_output_file_path, temp_dir_path, temp_file_id):

    import krusearch
    import krbioio
    import kralign
    import krseq

    records_dict = krbioio.read_sequence_file(fasta_file_path, 'fasta',
                                              ret_type='dict')

    cluster_dict = krusearch.parse_uc_file(uc_file_path)

    handle_aln = open(aln_output_file_path, 'w')
    handle_counts = open(counts_output_file_path, 'w')

    keys = cluster_dict.keys()
    keys.sort(key=lambda x: x, reverse=False)

    for key in keys:
        records = list()
        members = cluster_dict[key]
        spc = len(members)
        if spc >= min_seq_cluster and spc <= max_seq_cluster:
            handle_aln.write('>CLUSTER_' + str(key) + '\n')
            handle_counts.write('>CLUSTER_' + str(key) + '\n')
            if spc > 1:
                for m in members:
                    if m[0] == '+':
                        records.append(records_dict[m[1]])
                    else:
                        records.append(
                            krseq.reverse_complement(records_dict[m[1]]))
                aln = kralign.align(
                    records, 'muscle', 1, options='', temp_dir=temp_dir_path,
                    temp_file_id=temp_file_id)
                for l in range(0, aln.get_alignment_length()):
                    column = aln[:, l]
                    column = column.upper()
                    handle_aln.write(column + '\n')
                    counts = nucleotides_at_site(column)
                    counts_str = (
                        str(counts[0]) + '\t' +
                        str(counts[1]) + '\t' +
                        str(counts[2]) + '\t' +
                        str(counts[3]) + '\n'
                    )
                    handle_counts.write(counts_str)
            else:
                record = records_dict[members[0][1]]
                record = str(record.seq)
                for l in range(0, len(record)):
                    column = record[l]
                    column = column.upper()
                    handle_aln.write(column + '\n')
                    counts = nucleotides_at_site(column)
                    counts_str = (
                        str(counts[0]) + '\t' +
                        str(counts[1]) + '\t' +
                        str(counts[2]) + '\t' +
                        str(counts[3]) + '\n'
                    )
                    handle_counts.write(counts_str)

    handle_aln.close()
    handle_counts.close()


def nt_freq(nt_counts_file):
    import krio
    nt_counts = krio.read_table_file(
        path=nt_counts_file,
        has_headers=False,
        headers=['A', 'C', 'G', 'T'],
        delimiter='\t',
        quotechar='"',
        stripchar='',
        commentchar=">",
        rettype='dict'
    )

    total = 0

    t_a = 0
    t_c = 0
    t_g = 0
    t_t = 0

    for r in nt_counts:

        c_a = int(r['A'])
        c_c = int(r['C'])
        c_g = int(r['G'])
        c_t = int(r['T'])

        t_a = t_a + c_a
        t_c = t_c + c_c
        t_g = t_g + c_g
        t_t = t_t + c_t

        total = t_a + t_c + t_g + t_t

    # print(t_a, t_c, t_g, t_t)
    # print(float(t_a)/float(total),
    #       float(t_c)/float(total),
    #       float(t_g)/float(total),
    #       float(t_t)/float(total))

    return([float(t_a)/float(total),
            float(t_c)/float(total),
            float(t_g)/float(total),
            float(t_t)/float(total)])


def nt_site_counts(nt_counts_file):
    import krio
    nt_counts = krio.read_table_file(
        path=nt_counts_file,
        has_headers=False,
        headers=['A', 'C', 'G', 'T'],
        delimiter='\t',
        quotechar='"',
        stripchar='',
        commentchar=">",
        rettype='dict'
    )

    ret_value = list()

    for r in nt_counts:

        c_a = int(r['A'])
        c_c = int(r['C'])
        c_g = int(r['G'])
        c_t = int(r['T'])

        ret_value.append([c_a, c_c, c_g, c_t])

    return(ret_value)

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
    # Keep track of [A, C, G, T] configurations to not repeat unnecessary
    # calculations
    s_log = list()
    # repeats_found = 0
    for s in ns:
        # Check if this configuration has occured already
        repeat = False
        for sl in s_log:
            if (
                sl[0][0] == s[0] and
                sl[0][1] == s[1] and
                sl[0][2] == s[2] and
                sl[0][3] == s[3]
            ):
                ll = ll + sl[1]
                repeat = True
                # repeats_found = repeats_found + 1
                break

        if not repeat:
            l = like_homo_hetero(s, p, e, pi)
            nl = numpy.log(l)
            ll = ll + nl
            s_log.append([s, nl])

    ll = ll * (-1.0)

    print(ll)
    # print('Repeat configurations:', len(s_log))
    # print('Repeats:', repeats_found)
    return(ll)


def mle_e_and_pi(ns, p, e0, pi0):

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

    ml_est = optimize.fmin_l_bfgs_b(
        nll,
        x0=(e0, pi0),
        args=(ns, p),
        bounds=((1E-10, 0.99999), (1E-10, 0.99999)),
        approx_grad=True
    )

    return([ml_est[0][0], ml_est[0][1], ml_est[1]])


if __name__ == '__main__':
    # Tests
    import os
    ps = os.path.sep

    # p = nt_freq('/home/karolis/Dropbox/code/krpy/testdata/nt.counts')
    # print(p)
    # ns = nt_site_counts('/home/karolis/Dropbox/code/krpy/testdata/nt.counts')
    # mle = mle_e_and_pi(ns, p, e0=0.001, pi0=0.001)
    # print(mle)

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
