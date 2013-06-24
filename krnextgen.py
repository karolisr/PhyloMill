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


def mask_low_quality_sites(seq_str, qual_str, quality_score_treshold,
                           low_quality_residue='N'):

    # phred_dict = dict((chr(min(126, qp + 33)), qp) for qp in range(0, 93 + 1))

    phred_dict = {
        '$': 3, '(': 7, ',': 11, '0': 15, '4': 19, '8': 23, '<': 27, '@': 31,
        'D': 35, 'H': 39, 'L': 43, 'P': 47, 'T': 51, 'X': 55, '\\': 59,
        '`': 63, 'd': 67, 'h': 71, 'l': 75, 'p': 79, 't': 83, 'x': 87,
        '|': 91, '#': 2, "'": 6, '+': 10, '/': 14, '3': 18, '7': 22, ';': 26,
        '?': 30, 'C': 34, 'G': 38, 'K': 42, 'O': 46, 'S': 50, 'W': 54, '[': 58,
        '_': 62, 'c': 66, 'g': 70, 'k': 74, 'o': 78, 's': 82, 'w': 86, '{': 90,
        '"': 1, '&': 5, '*': 9, '.': 13, '2': 17, '6': 21, ':': 25, '>': 29,
        'B': 33, 'F': 37, 'J': 41, 'N': 45, 'R': 49, 'V': 53, 'Z': 57, '^': 61,
        'b': 65, 'f': 69, 'j': 73, 'n': 77, 'r': 81, 'v': 85, 'z': 89, '~': 93,
        '!': 0, '%': 4, ')': 8, '-': 12, '1': 16, '5': 20, '9': 24, '=': 28,
        'A': 32, 'E': 36, 'I': 40, 'M': 44, 'Q': 48, 'U': 52, 'Y': 56, ']': 60,
        'a': 64, 'e': 68, 'i': 72, 'm': 76, 'q': 80, 'u': 84, 'y': 88, '}': 92}

    quality_scores = [phred_dict[y] for y in qual_str]
    masked = ''
    for i, q in enumerate(quality_scores):
        bp = seq_str[i]
        if q < quality_score_treshold:
            bp = low_quality_residue
        masked = masked + bp
    return(masked)


# def mask_low_quality_sites(bio_seq_record, quality_score_treshold,
#                            low_quality_residue='N'):
#     from Bio.SeqRecord import SeqRecord
#     r = bio_seq_record
#     quality_scores = r.letter_annotations['phred_quality']
#     sequence = r.seq.tomutable()
#     for index, score in enumerate(quality_scores):
#         if score < quality_score_treshold:
#             sequence[index] = low_quality_residue
#     new_r = SeqRecord(seq=sequence.toseq(), id=r.id, name=r.name,
#                       description=r.description, dbxrefs=r.dbxrefs,
#                       features=r.features, annotations=r.annotations,
#                       letter_annotations=r.letter_annotations)
#     return(new_r)


def proportion_low_quality_sites(seq_str, low_quality_residue='N'):
    low_quality_sites = seq_str.count(low_quality_residue)
    sequence_length = len(seq_str)
    if sequence_length > 0:
        prop_lq_sites = float(low_quality_sites) / float(sequence_length)
    else:
        prop_lq_sites = 1.0
    return(prop_lq_sites)


# def proportion_low_quality_sites(bio_seq_record, low_quality_residue='N'):
#     sequence = str(bio_seq_record.seq)
#     low_quality_sites = sequence.count(low_quality_residue)
#     sequence_length = len(bio_seq_record.seq)
#     prop_lq_sites = float(low_quality_sites) / float(sequence_length)
#     return(prop_lq_sites)


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


def align_reads(r1, r2, mmmr_cutoff=0.85, ignore='N'):

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

    best_alignment = align_reads(r1, r2, mmmr_cutoff=mmmr_cutoff,
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

    ret_value = [o, r, message, cons, best_alignment]

    return(ret_value)


def bin_reads(title, f_seq_str, r_seq_str=None,
              max_prop_low_quality_sites=0.10, min_overlap=5, mmmr_cutoff=0.85,
              low_quality_residue='N', f_oligo=None, r_oligo=None,
              min_read_length=10):

    # from Bio import Seq
    # from Bio import SeqRecord

    f_hq = False
    r_hq = False

    # Look for sequencing oligonucleotides within the reads

    #   (match, total, ratio, (a, b), (c, d))
    #
    #   [a:b] - is the alignment range on the first sequence
    #   [c:d] - is the alignment range on the second sequence

    mmmr_cutoff_seq_oligo = 0.8 ###

    # Forward
    if r_oligo and f_seq_str:
        f_seq_oligo = align_reads(
            r_oligo, f_seq_str, mmmr_cutoff=mmmr_cutoff_seq_oligo,
            ignore=low_quality_residue)

        if ((f_seq_oligo[1] >= 10) or
           (f_seq_oligo[4][1] == len(f_seq_str)) and f_seq_oligo[1] >= 1):

            # print('F READ', f_seq_oligo)
            # print(f_seq_oligo[4][0] * ' ' + r_oligo)
            # print(f_seq_oligo[3][0] * ' ' + f_seq_str)

            f_seq_str = f_seq_str[0:f_seq_oligo[4][0]]

            # print(f_seq_oligo[3][0] * ' ' + f_seq_str)

        # print('F READ', f_seq_str)

    # Reverse
    if f_oligo and r_seq_str:
        r_seq_oligo = align_reads(
            f_oligo, r_seq_str, mmmr_cutoff=mmmr_cutoff_seq_oligo,
            ignore=low_quality_residue)

        if ((r_seq_oligo[1] >= 10) or
           (r_seq_oligo[4][0] == 0) and r_seq_oligo[1] >= 1):

            # print('R READ', r_seq_oligo)
            # print(r_seq_oligo[4][0] * ' ' + f_oligo)
            # print(r_seq_oligo[3][0] * ' ' + r_seq_str)

            r_seq_str = r_seq_str[r_seq_oligo[4][1]:len(r_seq_str)]

            # print((r_seq_oligo[3][0] + r_seq_oligo[4][1]) * ' ' + r_seq_str)

        # print('R READ', r_seq_str)

    # Determine the proportion of low quality sites
    # Check the length of the read
    if len(f_seq_str) < min_read_length:
        f_hq = False
    else:
        f_lq_sites = proportion_low_quality_sites(
            seq_str=f_seq_str,
            low_quality_residue=low_quality_residue)
        if f_lq_sites <= max_prop_low_quality_sites:
            f_hq = True
    if r_seq_str:
        # Check the length of the read
        if len(r_seq_str) < min_read_length:
            r_hq = False
        else:
            r_lq_sites = proportion_low_quality_sites(
                seq_str=r_seq_str,
                low_quality_residue=low_quality_residue)
            if r_lq_sites <= max_prop_low_quality_sites:
                r_hq = True

    # Produce F/R read consensus
    cons_message = ''
    consensus = None
    cons_title = None
    cons_alignment = None
    if f_hq and r_hq:
        consensus = consensus_fr_read(
            r1=f_seq_str,
            r2=r_seq_str,
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
        elif len(consensus[3]) < min_read_length:
            consensus = None
        else:
            cons_alignment = consensus[4]
            consensus = consensus[3]
            cons_title = title.split('|')[0] + '|C0NS:' + str(cons_message)
            # consensus = SeqRecord.SeqRecord(
            #     seq=cons_seq, id=cons_id, name='', description='')

    ret_value = [f_hq, r_hq, f_seq_str, r_seq_str, consensus, cons_title,
                 cons_message, cons_alignment]

    return(ret_value)


def _write_demultiplex_results_(barcodes,
                                reverse_reads_file_path,
                                result_batch_forward_other,
                                write_handle_forward_other,
                                result_batch_reverse_other,
                                write_handle_reverse_other
                                ):
    # from Bio import SeqIO

    for barcode in barcodes:
        # SeqIO.write(barcode['result_batch_forward'],
                    # barcode['write_handle_forward'], output_file_format)
        for r in barcode['result_batch_forward']:
            barcode['write_handle_forward'].write(r[0] + '\n')
            barcode['write_handle_forward'].write(r[1] + '\n')
            barcode['write_handle_forward'].write('+\n')
            barcode['write_handle_forward'].write(r[2] + '\n')
        del barcode['result_batch_forward'][:]
        if reverse_reads_file_path is not None:
            # SeqIO.write(barcode['result_batch_reverse'],
                        # barcode['write_handle_reverse'], output_file_format)
            for r in barcode['result_batch_reverse']:
                barcode['write_handle_reverse'].write(r[0] + '\n')
                barcode['write_handle_reverse'].write(r[1] + '\n')
                barcode['write_handle_reverse'].write('+\n')
                barcode['write_handle_reverse'].write(r[2] + '\n')
            del barcode['result_batch_reverse'][:]

    for r in result_batch_forward_other:
        write_handle_forward_other.write(r[0] + '\n')
        write_handle_forward_other.write(r[1] + '\n')
        write_handle_forward_other.write('+\n')
        write_handle_forward_other.write(r[2] + '\n')
    # SeqIO.write(result_batch_forward_other,
                # write_handle_forward_other, output_file_format)
    del result_batch_forward_other[:]

    if reverse_reads_file_path is not None:
        for r in result_batch_reverse_other:
            write_handle_reverse_other.write(r[0] + '\n')
            write_handle_reverse_other.write(r[1] + '\n')
            write_handle_reverse_other.write('+\n')
            write_handle_reverse_other.write(r[2] + '\n')
        # SeqIO.write(result_batch_reverse_other,
                    # write_handle_reverse_other, output_file_format)
        del result_batch_reverse_other[:]


def demultiplex(barcodes,
                forward_reads_file_path,
                reverse_reads_file_path=None,
                input_file_format='fastq',
                max_barcode_mismatch_count=1,
                output_dir='.',
                trim_barcode=True,
                trim_extra=0,
                write_every=1000
                ):

    import os
    import Levenshtein
    # from Bio import SeqIO
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    from Bio.Seq import Seq
    # from Bio import SeqRecord
    import krio
    # import krseq

    ps = os.path.sep
    output_dir = output_dir.rstrip(ps) + ps
    krio.prepare_directory(output_dir)

# >>> handle = open("Quality/tricky.fastq", "rU")
# >>> for (title, sequence, quality) in FastqGeneralIterator(handle):
# ...     print title
# ...     print sequence, quality

    # forward_reads = SeqIO.parse(forward_reads_file_path, input_file_format)
    forward_reads_handle = open(forward_reads_file_path, "rU")
    forward_reads = FastqGeneralIterator(forward_reads_handle)
    reverse_reads = None
    if reverse_reads_file_path is not None:
        # reverse_reads = SeqIO.parse(reverse_reads_file_path, input_file_format)
        reverse_reads_handle = open(reverse_reads_file_path, "rU")
        reverse_reads = FastqGeneralIterator(reverse_reads_handle)

    for barcode in barcodes:

        barcode['length'] = len(barcode['barcode'])
        base_file_name = output_dir + barcode['id'] + '_' + barcode['barcode']

        barcode['result_batch_forward'] = list()
        barcode['file_path_forward'] = (base_file_name + '_f.' + 'fastq')
        barcode['write_handle_forward'] = open(barcode['file_path_forward'],
                                               'wa')

        if reverse_reads_file_path is not None:
            barcode['result_batch_reverse'] = list()
            barcode['file_path_reverse'] = (base_file_name + '_r.' + 'fastq')
            barcode['write_handle_reverse'] = open(
                barcode['file_path_reverse'], 'wa')

    write_handle_forward_other = open(output_dir + 'Mismatch_f.' + 'fastq',
                                      'wa')
    result_batch_forward_other = list()
    write_handle_reverse_other = None
    result_batch_reverse_other = None
    if reverse_reads_file_path is not None:
        write_handle_reverse_other = open(output_dir + 'Mismatch_r.' + 'fastq',
                                          'wa')
        result_batch_reverse_other = list()

    # Loop over all reads
    for i, (f_title, f_seq, f_qual) in enumerate(forward_reads):
        # print(i+1, end='\r')
        # f_record = SeqRecord.SeqRecord(
        #     f_record.seq,
        #     f_record.description.replace(' ', '|'),
        #     '',
        #     '',
        #     f_record.dbxrefs,
        #     f_record.features,
        #     f_record.annotations,
        #     f_record.letter_annotations)

        # Try all barcodes
        barcode_match_found = False
        for barcode in barcodes:

            l = barcode['length']
            b = barcode['barcode'].lower()
            r = str(f_seq[0:l]).lower()
            ld = Levenshtein.distance(b, r)

            # Barcode match found
            if ld <= max_barcode_mismatch_count:
                barcode_match_found = True
                f_title = f_title.replace(' ', '|')
                f_title = '@' + f_title
                if trim_barcode:
                    # f_record = krseq.trim_residues(f_record, l, False)
                    f_seq = f_seq[l:]
                    f_qual = f_qual[l:]
                    if trim_extra > 0:
                        # f_record = krseq.trim_residues(f_record, trim_extra,
                        #                                False)
                        f_seq = f_seq[trim_extra:]
                        f_qual = f_qual[trim_extra:]
                # barcode['result_batch_forward'].append(f_record)
                f_read = (f_title, f_seq, f_qual)
                barcode['result_batch_forward'].append(f_read)

                # If there are reverse reads to consider
                if reverse_reads_file_path is not None:
                    r_title, r_seq, r_qual = reverse_reads.next()
                    # r_record = reverse_reads.next()
                    # r_record = krseq.reverse_complement(r_record)
                    # r_record = SeqRecord.SeqRecord(
                    #     r_record.seq,
                    #     r_record.description.replace(' ', '|'),
                    #     '',
                    #     '',
                    #     r_record.dbxrefs,
                    #     r_record.features,
                    #     r_record.annotations,
                    #     r_record.letter_annotations)

                    r_seq = Seq(r_seq)
                    r_seq = str(r_seq.reverse_complement())
                    r_qual = r_qual[::-1]

                    r_title = r_title.replace(' ', '|')
                    r_title = '@' + r_title

                    if trim_barcode and trim_extra > 0:
                        # r_record = krseq.trim_residues(r_record, trim_extra,
                        #                                True)

                        r_seq = r_seq[:-trim_extra]
                        r_qual = r_qual[:-trim_extra]

                    # barcode['result_batch_reverse'].append(r_record)
                    r_read = (r_title, r_seq, r_qual)
                    barcode['result_batch_reverse'].append(r_read)

                # Barcode match found, breakout of the loop
                break

        if not barcode_match_found:
            f_title = f_title.replace(' ', '|')
            f_title = '@' + f_title
            f_read = (f_title, f_seq, f_qual)
            result_batch_forward_other.append(f_read)
            if reverse_reads_file_path is not None:
                # r_record = reverse_reads.next()
                r_title, r_seq, r_qual = reverse_reads.next()
                r_title = r_title.replace(' ', '|')
                r_title = '@' + r_title
                r_read = (r_title, r_seq, r_qual)
                result_batch_reverse_other.append(r_read)

        if i % write_every == write_every - 1 and i > 0:
            _write_demultiplex_results_(barcodes,
                                        reverse_reads_file_path,
                                        result_batch_forward_other,
                                        write_handle_forward_other,
                                        result_batch_reverse_other,
                                        write_handle_reverse_other)

    _write_demultiplex_results_(barcodes,
                                reverse_reads_file_path,
                                result_batch_forward_other,
                                write_handle_forward_other,
                                result_batch_reverse_other,
                                write_handle_reverse_other)

    for barcode in barcodes:
        barcode['write_handle_forward'].close()
        if reverse_reads_file_path is not None:
            barcode['write_handle_reverse'].close()

    write_handle_forward_other.close()
    if reverse_reads_file_path is not None:
        write_handle_reverse_other.close()

    return()


def combine_demultiplexed_results(input_dir, output_dir):

    import os
    import shutil
    import krio

    ps = os.path.sep
    input_dir = input_dir.rstrip(ps) + ps
    output_dir = output_dir.rstrip(ps) + ps
    krio.prepare_directory(output_dir)
    directory_list = krio.parse_directory(
        path=input_dir,
        file_name_sep=' ',
        sort='forward'
    )

    file_dir_path = directory_list[0]['path'].rstrip(ps) + ps
    file_list = krio.parse_directory(
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
                   fasta_file_path, aln_clustal_phylip_file_path=None,
                   aln_output_file_path=None, counts_output_file_path=None,
                   program='mafft', options='', program_executable='mafft'):

    from Bio import SeqIO
    from Bio import AlignIO

    import krusearch
    # import krbioio
    import kralign
    import krseq
    # import krcl

    # Uses too much RAM
    # records_dict = krbioio.read_sequence_file(fasta_file_path, 'fasta',
    #                                           ret_type='dict')

    # print(program)

    records_dict = SeqIO.index(fasta_file_path, 'fasta')
    cluster_dict = krusearch.parse_uc_file(uc_file_path)

    handle_aln = None
    if aln_output_file_path:
        handle_aln = open(aln_output_file_path, 'w')

    handle_counts = None
    if counts_output_file_path:
        handle_counts = open(counts_output_file_path, 'w')

    # f_id = uc_file_path.split('.')[0].split('/')[-1].split('_')[0]

    alignments = list()

    keys = cluster_dict.keys()
    keys.sort(key=lambda x: x, reverse=False)
    # cluster_count = len(keys)

    cluster_depths = list()
    # krcl.hide_cursor()
    for i, key in enumerate(keys):
        records = list()
        members = cluster_dict[key]
        rpc = len(members)
        cluster_depths.append(rpc)
        if rpc >= min_seq_cluster and (rpc <= max_seq_cluster or
                                       max_seq_cluster == 0):
            # krcl.print_progress(i, cluster_count, 50, '')
            # print(f_id, i, '/', cluster_count)
            if handle_aln:
                handle_aln.write('>CLUSTER_' + str(key) + '\n')
            if handle_counts:
                handle_counts.write('>CLUSTER_' + str(key) + '\n')
            # X = False ####
            if rpc > 1:
                for m in members:
                    # if m[1] == "HWI-ST1155:95:D0KL1ACXX:4:1101:2261:6701|C0NS:0": ####
                    #     X = True ####
                    # if X: ####
                    #     print(m) ####
                    if m[0] == '+':
                        records.append(records_dict[m[1]])
                    else:
                        records.append(
                            krseq.reverse_complement(records_dict[m[1]]))
                # if X: ####
                #     print('---') ####
                aln = kralign.align(
                    records, program=program, options=options, program_executable=program_executable)
                    # options='--retree 1 --thread '+str(threads))
                # aln = kralign.align(records, 'muscle', options='')
                # if X: ####
                #     AlignIO.write(aln, '/home/karolis/Dropbox/Code/aln.fasta', "phylip-relaxed") ####
                if aln_clustal_phylip_file_path:
                    # import krother
                    # print(krother.attr(aln))
                    ids = list()
                    for seq in aln:
                        # print(seq.id)
                        seq.id = seq.id.split('_')[0] ###
                        ids.append(seq.id)
                    if len(set(ids)) == len(ids):
                        alignments.append(aln)
                    else:
                        print('Warning: multiple sequences from the same sample.')
                if handle_aln or handle_counts:
                    for l in range(0, aln.get_alignment_length()):
                        column = aln[:, l]
                        column = column.upper()
                        if handle_aln:
                            handle_aln.write(column + '\n')
                        counts = nucleotides_at_site(column)
                        counts_str = (
                            str(counts[0]) + '\t' +
                            str(counts[1]) + '\t' +
                            str(counts[2]) + '\t' +
                            str(counts[3]) + '\n'
                        )
                        if handle_counts:
                            handle_counts.write(counts_str)
            else:
                if handle_aln or handle_counts:
                    record = records_dict[members[0][1]]
                    record = str(record.seq)
                    for l in range(0, len(record)):
                        column = record[l]
                        column = column.upper()
                        if handle_aln:
                            handle_aln.write(column + '\n')
                        counts = nucleotides_at_site(column)
                        counts_str = (
                            str(counts[0]) + '\t' +
                            str(counts[1]) + '\t' +
                            str(counts[2]) + '\t' +
                            str(counts[3]) + '\n'
                        )
                        if handle_counts:
                            handle_counts.write(counts_str)
    # krcl.show_cursor()
    if handle_aln:
        handle_aln.close()
    if handle_counts:
        handle_counts.close()

    if aln_clustal_phylip_file_path:
        AlignIO.write(alignments, aln_clustal_phylip_file_path, "phylip-relaxed")

    return(cluster_depths)


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


def nt_site_counts(nt_counts_file, min_total_per_site=1, max_total_per_site=0,
                   rettype='list'):

    # from Bio import trie  # ##
    import string  # ##
    import datrie  # ##
    # import krio

    # commentchar = '>'
    # if rettype == 'dict':
    #     commentchar = '#'

    # nt_counts = krio.read_table_file(
    #     path=nt_counts_file,
    #     has_headers=False,
    #     # headers=['A', 'C', 'G', 'T'],
    #     headers=None,
    #     delimiter='\t',
    #     quotechar='"',
    #     stripchar='',
    #     commentchar=commentchar,
    #     rettype='list'  # This rettype is always a list
    # )

    input_handle = open(nt_counts_file, 'rb')
    nt_counts = input_handle.readlines()

    # print('Done ' + nt_counts_file)

    ret_value = list()

    if rettype == 'dict':
        # ret_value = dict()
        # ret_value = trie.trie()  # ##
        ret_value = datrie.Trie(string.printable)  # ##

    cluster_name = None

    for r in nt_counts:

        # print(r)

        r = r.strip()

        # if r[0].startswith('>'):
        if rettype == 'dict' and r.startswith('>'):
            cluster_name = unicode(r.split('>')[1])  # ## datrie wants unicode keys
            # print(cluster_name)
            ret_value[cluster_name] = list()
            # print(len(ret_value.keys()))
            # print(ret_value[cluster_name])
            continue

        if rettype == 'list' and r.startswith('>'):
            continue

        # c_a = int(r['A'])
        # c_c = int(r['C'])
        # c_g = int(r['G'])
        # c_t = int(r['T'])

        r = r.split('\t')

        c_a = int(r[0])
        c_c = int(r[1])
        c_g = int(r[2])
        c_t = int(r[3])

        t = c_a + c_c + c_g + c_t

        if (t >= min_total_per_site and
           (t <= max_total_per_site or max_total_per_site == 0)):
            if rettype == 'dict':
                ret_value[cluster_name].append((c_a, c_c, c_g, c_t))
                # print(cluster_name)
                # print(ret_value[cluster_name])
            else:
                ret_value.append((c_a, c_c, c_g, c_t))

    # print(ret_value)

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

    likelihood = (1.0-pi) * like_homo(s, p, e) + pi * like_hetero(s, p, e)
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

    if (e > 1.0 or e < 0) or (pi > 1.0 or pi < 0):
        return(float("inf"))

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
            if l > 0:
                nl = numpy.log(l)
                ll = ll + nl
                # print(ll)
                s_log.append([s, nl])

    ll = ll * (-1.0)
    # print(ll)
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

    ml_est = optimize.fmin(
        nll,
        x0=(e0, pi0),
        args=(ns, p),
        # xtol=0.0001,
        # ftol=0.0001,
        # maxiter=None,
        # maxfun=None,
        full_output=1,
        disp=0
        # retall=0,
        # callback=None
    )

    # ml_est = optimize.fmin_l_bfgs_b(
    #     nll,
    #     x0=(e0, pi0),
    #     args=(ns, p),
    #     bounds=((1E-10, 0.99999), (1E-10, 0.99999)),
    #     approx_grad=True
    # )

    ret_value = [ml_est[0][0], ml_est[0][1], ml_est[1]]

    # print(ret_value)

    return(ret_value)


def consensus_base(s, e, pi, p=0.95, low_quality_residue='N', min_total_per_site=4):

    '''
        Given nucleotide counts (from multiple NextGen reads) at a site,
        determine if the site is heterozygous and return nucleotides present.

        Parameters:
            s - A list of counts of nucleotides for a given site
                 0  1  2  3
                [A, C, G, T]

            e - Error rate

            pi - Nucleotide diversity

            p - Threshold probability value. When the relative probability of a
                base at site is below this value, the site will be called as
                low quality.

        Returns:
            (rel_prob, het, bases_at_site, consensus)

            rel_prob - relative probability of the base at site
            het - heterozygous (True/False)
            bases_at_site - list of bases at site
            consensus - consensus base or ambiguity
    '''

    from scipy import special
    from scipy import stats
    from heapq import nlargest

    import kriupac

    # Indexes of the two most common bases at the site
    # TODO What happens if there are more than two common bases?
    indexes = [0, 1, 2, 3]
    common = nlargest(2, indexes, key=lambda i: s[i])

    # print(common[0], common[1])

    bases = ['A', 'C', 'G', 'T']
    b1 = bases[common[0]]
    b2 = bases[common[1]]

    # print(b1, b2)

    k1 = s[common[0]]
    k2 = s[common[1]]
    n = s[common[0]] + s[common[1]]

    if n < min_total_per_site:
        ret_value = (0, False, (low_quality_residue, low_quality_residue), low_quality_residue)
        return(ret_value)

    # print(k1, k2, n)

    prob_het = special.binom(n, k1) / (2.0 ** n)
    prob_hom_1 = stats.binom.pmf(k1, n, e)
    prob_hom_2 = stats.binom.pmf(k2, n, e)

    prior_het = pi
    prior_hom = (1.0-pi) / 2.0

    prob_het = prob_het * prior_het
    prob_hom_1 = prob_hom_1 * prior_hom
    prob_hom_2 = prob_hom_2 * prior_hom

    probs = [prob_het, prob_hom_1, prob_hom_2]
    prob_max = max(probs)
    rel_prob = prob_max / sum(probs)

    bases_at_site = (b1, b1)
    consensus = b1

    het = False
    if probs.index(prob_max) == 0:
        het = True

    if rel_prob < p:
        bases_at_site = (low_quality_residue, low_quality_residue)
        consensus = low_quality_residue
    elif het:
        bases_at_site = [b1, b2]
        bases_at_site.sort()
        consensus = kriupac.IUPAC_DOUBLE_DNA_DICT[''.join(bases_at_site)]

    ret_value = (rel_prob, het, bases_at_site, consensus)

    return(ret_value)


if __name__ == '__main__':
    # Tests
    import os
    ps = os.path.sep

    # ns = nt_site_counts('/home/karolis/Dropbox/Code/krpy/testdata/nt.counts',
    #                     rettype='dict')

    # # consensus_base
    # s = [69, 3, 10, 600]
    # e = 0.001
    # pi = 0.01
    # cb = consensus_base(s, e, pi, p=0.95, low_quality_residue='N')
    # print(cb)

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
    #     input_dir='/home/karolis/Dropbox/code/test/rad/
    #02-demultiplexed-fastq-parts',
    #     output_dir='/home/karolis/Dropbox/code/test/rad/
    #03-demultiplexed-fastq-combined')

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
