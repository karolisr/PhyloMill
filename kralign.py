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
    from krpy import kriupac

    end_gap_letter = '#'

    uracil = False

    col_count = alignment.get_alignment_length()

    identities = list()
    cons_str = ''

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

    # Count bases in column.
    for col_idx in range(0, col_count):

        col_str = aln_column_str_list[col_idx]
        col_counts = dict()
        col_counts_expanded = dict()
        col_total = float()
        col_proportions = dict()
        col_cons_set = set()

        for letter in col_str:

            letter = letter.upper()

            if letter == 'U':
                uracil = True
                letter = 'T'

            factor = 1.0

            if letter == end_gap_letter:
                factor = end_gap_penalty

            if letter == unknown:
                factor = unknown_penalty

            if letter in kriupac.IUPAC_DNA_GAPS:
                factor = gap_penalty

            col_count = 1.0 * factor

            col_counts[letter] = col_counts.get(letter, 0) + col_count

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

        if len(col_cons_set) == 0:
            col_cons_set.add(unknown)

        if len(col_cons_set) > 1:
            identities.append(0)
        else:
            identities.append(1)

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

    if resolve_ambiguities:
        cons_str = krseq.resolve_ambiguities(cons_str)

    alphabet = generic_dna
    if uracil:
        cons_str = cons_str.replace('T', 'U')
        alphabet = generic_rna

    cons_seq = Seq.Seq(cons_str, alphabet)

    proportion_identical = float(sum(identities)) / float(len(identities))

    ret_value = (cons_seq, proportion_identical)

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
                results_dict[a_id].append(b_id)
                consumed_ids.append(b_id)

            print(a_id, ':', b_id, '=', score)

        # for k in results_dict.keys():
        #     print(k, results_dict[k])
        # print('=== === === === === === === ===')

    return results_dict


def consensus_DEPRECATED(alignment, threshold=0.0, unknown='N', resolve_ambiguities=False, gap_mismatch=False):
    from Bio import Seq
    import kriupac
    import krseq
    consensus = ''
    accepted_bases_at_sites = list()
    raw_counts_at_sites = list()
    identities = list()
    column_count = alignment.get_alignment_length()
    for column in range(0, column_count):
        # print(alignment[:, column])
        # Count individual characters in the column
        counts = dict()
        raw_count = 0
        for c in alignment[:, column]:
            c = c.upper()
            c.replace('U', 'T')
            if c in kriupac.IUPAC_DNA_CHARACTERS:
                counts[c] = counts.get(c, 0) + 1
                # if (c not in kriupac.IUPAC_DNA_GAPS) and (c not in kriupac.IUPAC_DNA_UNKNOWN):
                if (gap_mismatch or
                (((not gap_mismatch) and (c not in kriupac.IUPAC_DNA_GAPS)) and
                (c not in kriupac.IUPAC_DNA_UNKNOWN))):
                    raw_count = raw_count + 1
        raw_counts_at_sites.append(raw_count)
        # print(counts)
        counts_expanded = dict()
        for k in counts.keys():
            if k in kriupac.IUPAC_DNA_DICT_REVERSE:
                for char in kriupac.IUPAC_DNA_DICT_REVERSE[k]:
                    counts_expanded[char] = counts_expanded.get(char, 0) + counts[k]
            else:
                counts_expanded[k] = counts_expanded.get(k, 0) + counts[k]
        # print(counts_expanded)
        # Get the total characters in the column
        # TO DO: should N, and gaps be excluded?
        total = 0
        for k in counts_expanded.keys():
            # if (k not in kriupac.IUPAC_DNA_GAPS) and (k not in kriupac.IUPAC_DNA_UNKNOWN):
            if (gap_mismatch or
            (((not gap_mismatch) and (k not in kriupac.IUPAC_DNA_GAPS)) and
            (k not in kriupac.IUPAC_DNA_UNKNOWN))):
                total = total + counts_expanded[k]
        proportions = dict()
        for k in counts_expanded.keys():
            if (gap_mismatch or
            (((not gap_mismatch) and (k not in kriupac.IUPAC_DNA_GAPS)) and
            (k not in kriupac.IUPAC_DNA_UNKNOWN))):
                proportions[k] = float(counts_expanded[k]) / float(total)
        # print(proportions)
        # print('*** *** *** *** ***')

        site_set = set()
        if len(proportions) > 0 and threshold == 0:
            max_prop = max(proportions.values())
            for k in proportions.keys():
                if proportions[k] == max_prop:
                    site_set.add(k)
        else:
            for k in proportions.keys():
                if proportions[k] >= threshold:
                    site_set.add(k)

        if len(site_set) == 0:
            site_set.add(unknown)

        if len(site_set-set(unknown)) > 1:
            identities.append(0)
        else:
            identities.append(1)

        site_list = list(site_set)
        # print(site_list)
        site_list.sort()
        site_str = ''.join(site_list)
        site_str = site_str.replace('-', '')
        site_str = site_str.replace('.', '')

        if (unknown in site_str) and len(site_str) > 1:
            site_str = site_str.replace(unknown, '')

        site = unknown
        if site_str == unknown:
            site = unknown
        elif site_str == kriupac.IUPAC_DNA_STRING:
            site = unknown
        elif site_str == '':
            site = unknown
        else:
            site = kriupac.IUPAC_DNA_DICT[site_str]
        consensus = consensus + site
        accepted_bases_at_sites.append(site_set)

    if resolve_ambiguities:
        consensus = krseq.resolve_ambiguities(consensus)

    # print(consensus)
    # print(accepted_bases_at_sites)

    consensus = Seq.Seq(consensus)

    count_per_site = float(sum(raw_counts_at_sites)) / float(len(raw_counts_at_sites))

    len_identities = len(identities)

    if gap_mismatch:

        front_gap_boundary = None
        end_gap_boundary = None

        terminal_gap_window = 0
        number_of_unknown_in_range = 0

        number_of_internal_gaps = 0
        number_of_internal_gap_opens = 0

        terminal_gap_counter = 0
        accepted_bases_at_sites.reverse()
        for i, ab in enumerate(accepted_bases_at_sites):
            if end_gap_boundary is None:
                if '-' in ab:
                    # pass
                    if unknown in ab:
                        number_of_unknown_in_range = number_of_unknown_in_range + 1
                elif terminal_gap_counter <= terminal_gap_window:
                    terminal_gap_counter = terminal_gap_counter + 1
                else:
                    end_gap_boundary = i

        terminal_gap_counter = 0
        accepted_bases_at_sites.reverse()
        for i, ab in enumerate(accepted_bases_at_sites):
            if front_gap_boundary is None:
                if '-' in ab:
                    # pass
                    if unknown in ab:
                        number_of_unknown_in_range = number_of_unknown_in_range + 1
                elif terminal_gap_counter <= terminal_gap_window:
                    terminal_gap_counter = terminal_gap_counter + 1
                else:
                    front_gap_boundary = i

        gap_open_prev = False
        gap_open = False
        for i, ab in enumerate(accepted_bases_at_sites[front_gap_boundary:len(accepted_bases_at_sites)-end_gap_boundary]):
            if ('-' in ab) and (unknown not in ab):
                gap_open = True
                number_of_internal_gaps = number_of_internal_gaps + 1
            if unknown in ab:
                number_of_unknown_in_range = number_of_unknown_in_range + 1
            else:
                gap_open = False

            if gap_open != gap_open_prev:
                number_of_internal_gap_opens = number_of_internal_gap_opens + 1

            gap_open_prev = gap_open

        number_of_internal_gap_opens = number_of_internal_gap_opens / 2

        len_identities = len_identities - front_gap_boundary - end_gap_boundary + number_of_unknown_in_range - number_of_internal_gaps

        # print(sum(identities), len(identities), len_identities, number_of_internal_gap_opens)

    proportion_identical = float(sum(identities)) / float(len_identities)
    if proportion_identical > 1.0:
        proportion_identical = 1.0

    ret_value = (consensus, accepted_bases_at_sites, raw_counts_at_sites, count_per_site, identities, proportion_identical)

    return(ret_value)


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

#     # Tests

#     import os

#     PS = os.path.sep

#     import krbioio
#     # aln = krbioio.read_alignment_file('/Users/karolis/Desktop/aln_1.phy', 'phylip-relaxed')
#     # print()
#     # print(aln)
#     # print()
#     # print(consensus(alignment=aln))

#     recs = krbioio.read_sequence_file(
#         file_path='/Users/karolis/Desktop/Actinidia_chinensis__mRNA.gb',
#         file_format='genbank',
#         ret_type='list'
#         )

#     cluster(
#         records=recs,
#         threshold=0.95,
#         unknown='N',
#         key='gi',
#         aln_program='mafft',
#         aln_executable='mafft',
#         aln_options='--auto --reorder --adjustdirection')
