from __future__ import print_function
#from __future__ import unicode_literals


def concatenate(alignments, padding_length=0):

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
        aln1 = concatenate(alignments=alignments,
                           padding_length=padding_length)
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
    padding = SeqRecord(Seq('-' * padding_length, alphabet))
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
    return result_alignment


def align(records, program, threads, options='', temp_dir='.',
          temp_file_id='0'):

    import os
    import subprocess
    import krbioio

    working_dir = os.getcwd()

    temp_input_file = 'temp_to_be_aligned_' + temp_file_id + '.fasta'
    temp_output_file = 'temp_aligned_' + temp_file_id + '.fasta'

    os.chdir(temp_dir)

    krbioio.write_sequence_file(records, temp_input_file, 'fasta')

    if program == 'muscle':
        subprocess.call(program + ' -quiet -in ' + temp_input_file +
                        ' -out ' + temp_output_file, shell=True)

    elif (program == 'mafft') or (program == 'einsi') or (program == 'linsi'):
        prg_str = (program + ' --quiet --thread ' + str(threads) + ' ' +
                   options + ' ' + temp_input_file + ' > ' + temp_output_file)
        print(prg_str)
        subprocess.call(prg_str, shell=True)

    results = krbioio.read_alignment_file(temp_output_file, 'fasta')

    os.remove(temp_input_file)
    os.remove(temp_output_file)

    os.chdir(working_dir)

    return results
