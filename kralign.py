from __future__ import print_function
# from __future__ import unicode_literals


def concatenate(alignments, padding_length=0, partitions=None):

    '''
    Concatenate alignments based on the Seq ids; row order does not
    matter. If one alignment contains a Seq id that another one does
    not, gaps will be introduced in place of the missing Seq.

    Args:
        alignments: (tuple, list) Alignments to be concatenated.

        padding_length: Introduce this many gaps between concatenated
            alignments.
    '''

    from Bio import Alphabet
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import MultipleSeqAlignment
    if not isinstance(alignments, (list, tuple)):
        raise ValueError('Argument must be a list or a tuple.')
    elif len(alignments) == 1:
        return alignments[0]
    if isinstance(alignments, tuple):
        alignments = list(alignments)
    aln1 = None
    aln2 = None
    if len(alignments) > 2:
        aln2 = alignments.pop()
        result1 = concatenate(alignments=alignments,
                              padding_length=padding_length,
                              partitions=partitions)
        aln1 = result1[0]
        partitions = result1[1]
    elif len(alignments) == 2:
        aln1 = alignments[0]
        aln2 = alignments[1]
    if (not isinstance(aln1, MultipleSeqAlignment) or
            not isinstance(aln2, MultipleSeqAlignment)):
        raise ValueError(
            'Argument must inherit from Bio.Align.MultipleSeqAlignment.')
    alphabet = Alphabet._consensus_alphabet([aln1._alphabet, aln2._alphabet])
    aln1_dict = dict()
    aln2_dict = dict()
    for aln1_s in aln1:
        aln1_dict[aln1_s.id] = aln1_s
    for aln2_s in aln2:
        aln2_dict[aln2_s.id] = aln2_s
    aln1_length = aln1.get_alignment_length()
    aln2_length = aln2.get_alignment_length()
    aln1_gaps = SeqRecord(Seq('-' * aln1_length, alphabet))
    aln2_gaps = SeqRecord(Seq('-' * aln2_length, alphabet))
    padding = SeqRecord(Seq('N' * padding_length, alphabet))

    if not partitions:
        partitions = [(1, aln1_length)]
    partitions.append((1 + aln1_length, padding_length + aln1_length + aln2_length))

    result_seq_list = list()
    for aln1_key in aln1_dict.keys():
        merged_Seq = None
        if aln1_key in aln2_dict:
            merged_Seq = aln1_dict[aln1_key] + padding + aln2_dict[aln1_key]
            merged_Seq.id = aln1_dict[aln1_key].id
            merged_Seq.name = ''
            merged_Seq.description = ''
            aln2_dict.pop(aln1_key)
        else:
            aln1_seq_record = aln1_dict[aln1_key]
            merged_Seq = aln1_seq_record + padding + aln2_gaps
            merged_Seq.id = aln1_seq_record.id
            merged_Seq.name = ''
            merged_Seq.description = ''
        result_seq_list.append(merged_Seq)
    for aln2_seq_record in aln2_dict.values():
        merged_Seq = aln1_gaps + padding + aln2_seq_record
        merged_Seq.id = aln2_seq_record.id
        merged_Seq.name = ''
        merged_Seq.description = ''
        result_seq_list.append(merged_Seq)
    result_alignment = MultipleSeqAlignment(result_seq_list, alphabet)
    result_alignment.sort()
    return((result_alignment, partitions))


def align(records, program, options='', program_executable=''):

    import subprocess
    from StringIO import StringIO
    from Bio import AlignIO
    from Bio import SeqIO
    import shlex

    input_handle = StringIO()
    SeqIO.write(records, input_handle, 'fasta')

    args = None

    options = shlex.split(options)

    if program_executable == '':
        program_executable = program

    if program == 'muscle':
        args = [program_executable, '-quiet'] + options + ['-in', '-', '-out', '-']

    elif program == 'mafft':
        args = [program_executable, '--quiet'] + options + ['-']

    if program == 'clustalo':
        args = [program_executable] + options + ['-i', '-']

    alignment = None

    if args:
        # print(args)
        pipe = subprocess.Popen(
            args=args,
            bufsize=0,
            executable=None,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            preexec_fn=None,
            close_fds=False,
            shell=False,
            cwd=None,
            env=None,
            universal_newlines=True,
            startupinfo=None,
            creationflags=0)

        data = pipe.communicate(input=input_handle.getvalue())
        alignment = AlignIO.read(StringIO(data[0]), 'fasta')

    return alignment


def pairwise_identity(

    alignment,
    unknown_letters=set(['N']),
    unknown_id=0.0,
    free_unknowns=True,
    gap_id=0.0,
    free_gaps=True,
    end_gap_id=0.0,
    free_end_gaps=True):

    import sys

    from krpy import kriupac

    if len(alignment) != 2:
        print('Alignment must contain exactly two sequences.')
        sys.exit(1)

    end_gap_letter = '#'
    col_count = alignment.get_alignment_length()

    # Produce a list of string representations of the sequences in alignment.
    # Leading and trailing gaps will be replaced with term_gap_letter.
    aln_seq_str_list = list()
    for aln_seq in alignment:
        aln_str = str(aln_seq.seq)
        aln_str_l_strip = aln_str.lstrip(kriupac.IUPAC_DNA_GAPS_STRING)
        left_gap_count = len(aln_str) - len(aln_str_l_strip)
        aln_str_l_r_strip = aln_str_l_strip.rstrip(kriupac.IUPAC_DNA_GAPS_STRING)
        right_gap_count = len(aln_str_l_strip) - len(aln_str_l_r_strip)
        aln_str_term_gaps = left_gap_count * end_gap_letter + aln_str_l_r_strip + right_gap_count * end_gap_letter
        aln_seq_str_list.append(aln_str_term_gaps)

    # Produce a list of alignment column strings.
    aln_column_str_list = list()
    for col_idx in range(0, col_count):
        aln_column_str = ''
        for aln_seq_str in aln_seq_str_list:
            aln_column_str = aln_column_str + aln_seq_str[col_idx]
        aln_column_str_list.append(aln_column_str)

    # print('--- --- --- --- --- --- --- --- --- --- --- ---')

    score_list = list()
    weights_list = list()

    for col_idx in range(0, col_count):
        col_str = aln_column_str_list[col_idx]

        l1 = col_str[0]
        l2 = col_str[1]

        if l1 in kriupac.IUPAC_DNA_DICT_REVERSE.keys():
            l1 = kriupac.IUPAC_DNA_DICT_REVERSE[l1]

        if l2 in kriupac.IUPAC_DNA_DICT_REVERSE.keys():
            l2 = kriupac.IUPAC_DNA_DICT_REVERSE[l2]

        l1 = set(l1)
        l2 = set(l2)

        #

        end_gap_in_l1 = False
        end_gap_in_l2 = False
        end_gap_in_col = False

        if end_gap_letter in l1:
            end_gap_in_l1 = True
        if end_gap_letter in l2:
            end_gap_in_l2 = True

        if end_gap_in_l1 or end_gap_in_l2:
            end_gap_in_col = True

        #

        gap_in_l1 = False
        gap_in_l2 = False
        gap_in_col = False

        for g in list(kriupac.IUPAC_DNA_GAPS):

            if g in l1:
                gap_in_l1 = True
            if g in l2:
                gap_in_l2 = True

        if gap_in_l1 or gap_in_l2:
            gap_in_col = True

        #

        unknown_in_l1 = False
        unknown_in_l2 = False
        unknown_in_col = False

        for u in list(unknown_letters):

            if u in l1:
                unknown_in_l1 = True
            if u in l2:
                unknown_in_l2 = True

        if unknown_in_l1 or unknown_in_l2:
            unknown_in_col = True

        #

        score = 0.0
        weight = 0.0

        if end_gap_in_col and gap_in_col:
            weight = 0.0

        elif unknown_in_l1 and unknown_in_l2:
            weight = 0.0

        elif not free_end_gaps and end_gap_in_col:
            score = end_gap_id
            weight = 1.0

        elif not free_gaps and gap_in_col:
            score = gap_id
            weight = 1.0

        elif not free_unknowns and unknown_in_col:
            score = unknown_id
            weight = 1.0

        elif (not end_gap_in_col) and (not gap_in_col) and (not unknown_in_col):
            intersection = l1 & l2
            union = l1 | l2
            score = float(len(intersection)) / float(len(union))
            weight = 1.0

        score_list.append(score)
        weights_list.append(weight)

        # print(l1, l2, score, weight)

        # print('--- --- --- --- --- --- --- --- --- --- --- ---')

    pair_id = sum(score_list) / sum(weights_list)

    # print(pair_id)

    return pair_id


def identity(

    alignment,
    unknown_letters=set(['N']),
    unknown_id=0.0,
    free_unknowns=True,
    gap_id=0.0,
    free_gaps=True,
    end_gap_id=0.0,
    free_end_gaps=True):

    from Bio.Align import MultipleSeqAlignment

    row_count = len(alignment)

    pair_id_list = list()

    for i in range(0, row_count):
        for j in range(0, row_count):
            if i == j:
                continue

            aln = MultipleSeqAlignment(records=[alignment[i], alignment[j]])

            pair_id = pairwise_identity(
                alignment=aln,
                unknown_letters=unknown_letters,
                unknown_id=unknown_id,
                free_unknowns=free_unknowns,
                gap_id=gap_id,
                free_gaps=free_gaps,
                end_gap_id=end_gap_id,
                free_end_gaps=free_end_gaps)

            # print(alignment[i].id, alignment[j].id, pair_id)

            pair_id_list.append(pair_id)

    ident = sum(pair_id_list) / len(pair_id_list)

    return ident


def consensus(
    alignment,
    threshold=0.0,
    unknown='N',
    unknown_penalty=0.0,
    resolve_ambiguities=False,
    gap_penalty=0.0,
    end_gap_penalty=0.0
    ):

    from Bio import Seq
    from Bio.Alphabet import generic_dna
    from Bio.Alphabet import generic_rna
    from krpy import krseq
    from krpy import kriupac

    end_gap_letter = '#'

    uracil = False

    col_count = alignment.get_alignment_length()
    row_count = len(alignment)

    cons_str = ''
    pairwise_list = list()
    pairwise_list_weights = list()

    # Produce a list of string representations of the sequences in alignment.
    # Leading and trailing gaps will be replaced with term_gap_letter.
    aln_seq_str_list = list()
    for aln_seq in alignment:
        aln_str = str(aln_seq.seq)
        aln_str_l_strip = aln_str.lstrip(kriupac.IUPAC_DNA_GAPS_STRING)
        left_gap_count = len(aln_str) - len(aln_str_l_strip)
        aln_str_l_r_strip = aln_str_l_strip.rstrip(kriupac.IUPAC_DNA_GAPS_STRING)
        right_gap_count = len(aln_str_l_strip) - len(aln_str_l_r_strip)
        aln_str_term_gaps = left_gap_count * end_gap_letter + aln_str_l_r_strip + right_gap_count * end_gap_letter
        aln_seq_str_list.append(aln_str_term_gaps)

    # Produce a list of alignment column strings.
    aln_column_str_list = list()
    for col_idx in range(0, col_count):
        aln_column_str = ''
        for aln_seq_str in aln_seq_str_list:
            aln_column_str = aln_column_str + aln_seq_str[col_idx]
        aln_column_str_list.append(aln_column_str)

    for col_idx in range(0, col_count):

        col_str = aln_column_str_list[col_idx]
        col_str_clean = ''
        col_counts = dict()
        col_counts_expanded = dict()
        col_total = float()
        col_proportions = dict()
        col_cons_set = set()

        # Count bases in column.
        for letter in col_str:

            letter = letter.upper()

            if letter == 'U':
                uracil = True
                letter = 'T'

            # factor = 1.0

            # if letter == end_gap_letter:
            #     factor = end_gap_penalty

            # if letter == unknown:
            #     factor = unknown_penalty

            # if letter in kriupac.IUPAC_DNA_GAPS:
            #     factor = gap_penalty

            # col_count = 1.0 * factor

            # col_counts[letter] = col_counts.get(letter, 0) + col_count

            col_counts[letter] = col_counts.get(letter, 0) + 1.0
            col_str_clean = col_str_clean + letter

        for k in col_counts.keys():
            if k in kriupac.IUPAC_DNA_DICT_REVERSE:
                for letter in kriupac.IUPAC_DNA_DICT_REVERSE[k]:
                    col_counts_expanded[letter] = col_counts_expanded.get(letter, 0) + col_counts[k]
            else:
                col_counts_expanded[k] = col_counts_expanded.get(k, 0) + col_counts[k]

        for k in col_counts_expanded.keys():
            base_count = col_counts_expanded[k]
            col_total = col_total + base_count

        for k in col_counts_expanded.keys():
            base_count = col_counts_expanded[k]

            base_prop = 0.0
            if col_total > 0.0:
                base_prop = base_count / col_total

            col_proportions[k] = base_prop

        # Keep only the bases that occur at a high enough frequency
        if len(col_proportions) > 0 and threshold == 0:
            max_prop = max(col_proportions.values())
            if max_prop != 0.0:
                for k in col_proportions.keys():
                    if col_proportions[k] == max_prop:
                        col_cons_set.add(k)
        else:
            for k in col_proportions.keys():
                if col_proportions[k] >= threshold:
                    col_cons_set.add(k)

        for g in kriupac.IUPAC_DNA_GAPS:
            if g in col_cons_set:
                col_cons_set.remove(g)

        if end_gap_letter in col_cons_set:
            col_cons_set.remove(end_gap_letter)

        if len(col_cons_set) == 0:
            col_cons_set.add(unknown)

        col_cons_list = list(col_cons_set)
        col_cons_list.sort()
        col_str_new = ''.join(col_cons_list)

        if (unknown in col_str_new) and len(col_str_new) > 1:
            col_str_new = col_str_new.replace(unknown, '')

        site = unknown
        if col_str_new == unknown:
            site = unknown
        elif col_str_new == kriupac.IUPAC_DNA_STRING:
            site = unknown
        else:
            site = kriupac.IUPAC_DNA_DICT[col_str_new]

        cons_str = cons_str + site

        # Calculate pairwise identities
        do_not_compare = set([end_gap_letter, unknown]) | kriupac.IUPAC_DNA_GAPS
        same = 0
        diff = 0
        for i, l_1 in enumerate(col_str_clean):
            for j, l_2 in enumerate(col_str_clean):
                if i != j:

                    if (l_1 in do_not_compare) and (l_2 in do_not_compare):
                        continue

                    if l_1 == l_2:
                        same = same + 1.0

                    else:

                        if (l_1 == end_gap_letter) or (l_2 == end_gap_letter):
                            diff = diff + end_gap_penalty
                            # if end_gap_penalty > 0:
                            #     same = same + 1.0 - end_gap_penalty

                        elif (l_1 in kriupac.IUPAC_DNA_GAPS) or (l_2 in kriupac.IUPAC_DNA_GAPS):
                            diff = diff + gap_penalty
                            # if gap_penalty > 0:
                            #     same = same + 1.0 - gap_penalty

                        elif (l_1 == unknown) or (l_2 == unknown):
                            diff = diff + unknown_penalty
                            # if unknown_penalty > 0:
                            #     same = same + 1.0 - unknown_penalty
                        else:
                            diff = diff + 1

        pairwise = 0.0
        total = float(same + diff)
        if total > 0.0:
            pairwise = float(same) / total

        site_end_gap_count = col_str_clean.count(end_gap_letter)

        site_gap_count = 0
        for g in kriupac.IUPAC_DNA_GAPS:
            site_gap_count = site_gap_count + col_str_clean.count(g)

        site_unknown_count = col_str_clean.count(unknown)

        site_nt_count = len(col_str_clean) - site_end_gap_count - site_gap_count - site_unknown_count

        site_row_count = site_nt_count + (site_end_gap_count * (1-end_gap_penalty)) + (site_gap_count * (1-gap_penalty)) + (site_unknown_count * (1-unknown_penalty))
        pairwise = pairwise * site_row_count

        pairwise_list.append(pairwise)

        site_row_weight_count = site_nt_count + (site_end_gap_count * end_gap_penalty) + (site_gap_count * gap_penalty) + (site_unknown_count * unknown_penalty)
        pairwise_weight = 1
        if site_row_weight_count <= 1:
            pairwise_weight = 0

        pairwise_list_weights.append(pairwise_weight)

        # print(col_str_clean)
        # print('same:', same, 'diff:', diff, 'site_row_count:', site_row_count, 'pairwise:', pairwise, 'pairwise_weight:', pairwise_weight)

        # print('--- --- ---')

    pairwise_entire = sum(pairwise_list) / sum(pairwise_list_weights) / row_count

    # print(pairwise_entire)

    if resolve_ambiguities:
        cons_str = krseq.resolve_ambiguities(cons_str)

    alphabet = generic_dna
    if uracil:
        cons_str = cons_str.replace('T', 'U')
        alphabet = generic_rna

    cons_seq = Seq.Seq(cons_str, alphabet)

    ret_value = (cons_seq, pairwise_entire)

    return ret_value


def cluster(
    records,
    threshold=0.95,
    unknown='N',
    key='gi',
    aln_program='mafft',
    aln_executable='mafft',
    aln_options='--auto --reorder --adjustdirection'):

    results_dict = dict()
    consumed_ids = list()

    records = sorted(records, key=lambda x: len(x.seq), reverse=True)

    for a_rec in records:

        # print(a_rec.id, len(a_rec.seq))

        key_value = None
        if key == 'accession':
            key_value = a_rec.id
        elif key == 'gi':
            key_value = a_rec.annotations['gi']
        elif key == 'description':
            key_value = a_rec.description
        else:
            key_value = a_rec.id

        a_id = key_value

        if a_id in consumed_ids:
            continue

        results_dict[a_id] = list()
        results_dict[a_id].append(a_id)
        consumed_ids.append(a_id)

        for b_rec in records:

            key_value = None
            if key == 'accession':
                key_value = b_rec.id
            elif key == 'gi':
                key_value = b_rec.annotations['gi']
            elif key == 'description':
                key_value = b_rec.description
            else:
                key_value = b_rec.id

            b_id = key_value

            if a_id == b_id:
                continue

            if b_id in consumed_ids:
                continue

            aln = align(
                records=[a_rec, b_rec],
                program=aln_program,
                options=aln_options,
                program_executable=aln_executable)

            direction = '+'
            for a in aln:
                if a.id.startswith('_R_'):
                    direction = '-'
                    break

            cons = consensus(
                alignment=aln,
                threshold=0.000001,
                unknown=unknown,
                unknown_penalty=0.0,
                resolve_ambiguities=False,
                gap_penalty=0.0,
                end_gap_penalty=0.0
                )

            score = cons[1]

            if score >= threshold:
                results_dict[a_id].append([direction, b_id, score])
                consumed_ids.append(b_id)

            print(a_id, ':', b_id, '=', score)

        # for k in results_dict.keys():
        #     print(k, results_dict[k])
        # print('=== === === === === === === ===')

    return results_dict


def determine_conserved_regions(alignment_file, matrix, window, min_length, cutoff):

    import subprocess
    import csv
    import os

    directory = os.path.split(alignment_file)[0]
    cons_scores = directory + os.path.sep + 'conservation_scores.tsv'

    subprocess.call('score_conservation.py -o '+cons_scores+' -m /usr/local/conservation_code/matrix/'+matrix+'.bla -w '+str(window)+' '+alignment_file, shell=True)

    cons_csv = csv.reader(open(cons_scores, 'rb'), delimiter='\t')

    regions = []
    region = []

    for row in cons_csv:
        if row[0].startswith('#'):
            continue

        pos = int(row[0])+1
        con = float(row[1])

        if con >= float(cutoff):
            region.append([pos,con])
        else:
            if len(region) >= min_length:
                regions.append(region)
            region = []

    if len(region) >= min_length:
        regions.append(region)

    print('There are '+str(len(regions))+' conserved regions.')

    # for region in regions:
    #     print('----------------')
    #     print('There are '+str(len(region))+' residues in this region.')
    #     for position in region:
    #         print(position)

    return regions


def slice_out_conserved_regions(regions, alignment_file, name_prefix, output_dir_path):

    from Bio import AlignIO
    from Bio.Align import MultipleSeqAlignment
    import subprocess
    import os

    # directory = os.path.split(alignment_file)[0]
    directory = output_dir_path.strip(os.path.sep)
    alignment = AlignIO.read(open(alignment_file), "fasta")

    for i in range(0,len(regions)):

        region = regions[i]
        start = region[0][0]-1
        stop = region[-1][0]
        name = name_prefix + str(i+1)
        sliced_alignment = alignment[:,start:stop]
        sliced_alignment_edited = MultipleSeqAlignment(None)

        output_path = directory + os.path.sep + name + '.fasta'

        for record in sliced_alignment:
            if not "-" in str(record.seq):
                sliced_alignment_edited.append(record)

        AlignIO.write(sliced_alignment_edited, output_path, "fasta")
        subprocess.call('usearch -quiet -minseqlength 1 -derep_fulllength '+output_path+' -output '+output_path, shell=True)

        sliced_alignment_new = AlignIO.read(open(output_path), "fasta")

        j=1

        for record in sliced_alignment_new:
            record.id = name+'_'+str(j)
            record.description = ''
            record.name = ''
            j = j+1

        AlignIO.write(sliced_alignment_new, output_path, "fasta")

    return


# if __name__ == '__main__':

    # # Tests

    # import os

    # PS = os.path.sep

    # import krbioio

    # aln = krbioio.read_alignment_file('/Users/karolis/Desktop/aln_1.phy', 'phylip-relaxed')

    # ident = identity(

    #     alignment=aln,
    #     unknown_letters=set(['N']),
    #     unknown_id=0.0,
    #     free_unknowns=True,
    #     gap_id=0.0,
    #     free_gaps=True,
    #     end_gap_id=0.0,
    #     free_end_gaps=True)

    # print(ident)

    # pid = pairwise_identity(

    #     alignment=aln,
    #     unknown_letters=set(['N']),
    #     unknown_id=0.0,
    #     free_unknowns=True,
    #     gap_id=0.0,
    #     free_gaps=True,
    #     end_gap_id=0.0,
    #     free_end_gaps=True)

    # print(pid)

    # cons = consensus(
    #     alignment=aln,
    #     threshold=0.1,
    #     unknown='N',
    #     unknown_penalty=0.0,
    #     resolve_ambiguities=False,
    #     gap_penalty=0.0,
    #     end_gap_penalty=0.0)

    # print(cons)

    # recs = krbioio.read_sequence_file(
    #     file_path='/Users/karolis/Desktop/Actinidia_chinensis__mRNA.gb',
    #     file_format='genbank',
    #     ret_type='list'
    #     )

    # cluster(
    #     records=recs,
    #     threshold=0.95,
    #     unknown='N',
    #     key='gi',
    #     aln_program='mafft',
    #     aln_executable='mafft',
    #     aln_options='--auto --reorder --adjustdirection')
