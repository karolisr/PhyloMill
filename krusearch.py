from __future__ import print_function
#from __future__ import unicode_literals

'''
    This module is designed for USEARCH 6. There are many differences between
    USEARCH versions 6 and 5. This module is incompatible with version 5.
'''


def cluster_file(
    input_file_path,
    identity_threshold,
    output_file_path,

    consensus_file_path=False,

    sorted_input=False,
    algorithm='fast',  # fast smallmem
    strand='plus',  # plus both
    threads=1,
    quiet=True,
    program='usearch',
    heuristics=True,
    query_coverage=0.5,
    target_coverage=0.5,

    sizein=False,
    sizeout=False,
    usersort=False
):

    # import os
    import subprocess
    import os

    if quiet:
        quiet = ' -quiet'
    else:
        quiet = ''

    if algorithm == 'smallmem' and not sorted_input:
        ifp_split = os.path.splitext(input_file_path)
        subprocess.call(
            (program + quiet +
             ' -sortbylength ' + input_file_path +
             ' -output ' + ifp_split[0] + '_sorted' + ifp_split[1]
             ),
            shell=True)
        input_file_path = ifp_split[0] + '_sorted' + ifp_split[1]

    if algorithm == 'fast':
        threads = ' -threads ' + str(threads)
    else:
        threads = ''

    if heuristics:
        heuristics = ''
    else:
        heuristics = ' -fulldp'

    if sizein:
        sizein = ' -sizein'
    else:
        sizein = ''

    if sizeout:
        sizeout = ' -sizeout'
    else:
        sizeout = ''

    if usersort:
        usersort = ' -usersort'
    else:
        usersort = ''

    if consensus_file_path:
        consensus = ' -consout ' + consensus_file_path
    else:
        consensus = ''

    command = (
        program + quiet + heuristics + sizein + sizeout + usersort +
        ' -query_cov ' + str(query_coverage) +
        ' -target_cov ' + str(target_coverage) +
        ' -cluster_' + algorithm + ' ' + input_file_path +
        ' -strand ' + strand +
        ' -id ' + str(identity_threshold) +
        threads +
        ' -uc ' + output_file_path +
        consensus
    )

    # print(command)

    subprocess.call(command, shell=True)

    # if algorithm == 'smallmem' and not sorted_input:
    #     os.remove(input_file_path)

    return output_file_path


def parse_uc_file(uc_file_path, key='clust_number'):

    '''
        UC file

        Edited, from:
            http://www.drive5.com/usearch/manual/ucout.html

        USEARCH cluster format (UC) is a tab-separated text file. By
        convention, the .uc filename extension is used. Each line is either a
        comment (starts with #) or a record. Every input sequence generates one
        record (H, S or N); additional record types give information about
        clusters. If an input sequence matched a target sequence, then the
        alignment and the identity computed from that alignment are also
        provided. Fields that do not apply to a given record type are filled
        with an asterisk placeholder (*).

        Record  Description

        H       Hit. Represents a query-target alignment. For clustering,
                    indicates the cluster assignment for the query.

        S       Centroid (clustering only). There is one S record for each
                    cluster, this gives the centroid (representative) sequence
                    label in the 9th field. Redundant with the C record;
                    provided for backwards compatibility.

        C       Cluster record (clustering only). The 3rd field is set to the
                    cluster size (number of sequences in the cluster) and the
                    9th field is set to the label of the centroid sequence.

        Field   Description

        0       Record type S, H, C or N (see table below).
        1       Cluster number (0-based).
        2       Sequence length (S, N and H) or cluster size (C).
        3       For H records, percent identity with target.
        4       For H records, the strand: + or - for nucleotides,
                    . for proteins.
        5       Not used. Included for backwards compatibility.
        6       Not used. Included for backwards compatibility.
        7       Compressed alignment or the symbol '=' (equals sign).
                    The = indicates that the query is 100% identical to the
                    target sequence (field 10).
        8       Label of query sequence (always present).
        9       Label of target sequence (H records only).


        Compressed Alignments

        Edited, from:
            http://drive5.com/usearch/manual/compressed_alignments.html

        A compressed alignments represents an alignment in a compact format
        that does not include the sequence letters. The representation uses
        run-length encoding, as follows. Each column in the alignment is
        classified as M, D or I.

        Column  Description
        M       Match. A pair of letters.
        D       Delete. A gap in the target.
        I       Insert. A gap in the query.

        If there are n consecutive columns of type C, this is represented as
        nC. For example, 123M is 123 consecutive matches. As a special case, if
        n=1 then n is omitted. So for example, D5M2I3M represents an alignment
        of this form:

           Query    XXXXXX--XXX
           Target   -XXXXXXXXXX
           Column   DMMMMMIIMMM
    '''

    import csv
    import string
    import datrie

    uc = csv.reader(open(uc_file_path, 'rb'), delimiter=b'\t')

    # cluster_dict = {}
    cluster_dict = datrie.Trie(string.printable)

    for row in uc:

        rec_type = row[0]
        # clust_number = int(row[1])
        clust_number = unicode(row[1])
        frac_id = row[3]
        if frac_id != '*':
            frac_id = float(frac_id)
        strand = row[4]
        if strand == '*':
            strand = '+'
        #comp_aln = row[7]
        query = unicode(row[8])
        target = unicode(row[9])

        if rec_type == 'S':
            if key == 'clust_number':
                cluster_dict[clust_number] = [[strand, str(query), 100.0]]
            elif key == 'centroid':
                cluster_dict[query] = [[strand, str(query), 100.0]]

        elif rec_type == 'H':
            if key == 'clust_number':
                cluster_dict[clust_number].append([strand, str(query), frac_id])
            elif key == 'centroid':
                cluster_dict[target].append([strand, str(query), frac_id])

    return cluster_dict


def write_uc_file(cluster_dict, uc_file_path):

    # Writes a semi-compatible uc file

    handle = open(uc_file_path, 'wb')
    for c, k in enumerate(cluster_dict.keys()):
        seed_name = None
        for l, seq in enumerate(cluster_dict[k]):
            if l == 0:
                seed_name = str(seq[1])
            record = "\t*\t" + str(seq[2]) + "\t" + str(seq[0]) + "\t*\t*\t*\t" + str(seq[1]) + "\t"
            if l == 0:
                record = 'S\t' + str(k) + record + "*"
            else:
                record = 'H\t' + str(k) + record + seed_name
            handle.write(record + '\n')
    handle.close()


def cluster_records(records, similarity, temp_dir):

    '''
    Cluster records based on similarity treshold.
    '''

    import os
    import krbioio

    ps = os.path.sep

    to_cluster_path = temp_dir + ps + 'to_cluster_temp.fasta'
    clustered_path = temp_dir + ps + 'clustered_temp.uc'

    krbioio.write_sequence_file(records, to_cluster_path, 'fasta')
    cluster_file(to_cluster_path, clustered_path, similarity)
    cluster_dict = parse_uc_file(clustered_path)

    os.remove(to_cluster_path)
    os.remove(clustered_path)

    return cluster_dict


def decode_compressed_alignment(comp_aln):

    # >>> import re
    # >>> s = "5d4h2s"
    # >>> p = re.compile("([0-9])([a-z])")
    # >>> for m in p.findall(s):
    # ...   print m
    # ...
    # ('5', 'd')
    # ('4', 'h')
    # ('2', 's')

    idxs_in_str = lambda x: [i for i, ltr in enumerate(x[0]) if ltr in x[1]]
    category_idxs = idxs_in_str((comp_aln, 'MDI'))

    ret_value = list()

    for i, j in enumerate(category_idxs):
        count = 1
        if j > 0:
            if i == 0:
                count = comp_aln[0:j]
            else:
                count = comp_aln[category_idxs[i-1]+1:j]
        if count == '':
            count = 1
        count = int(count)
        ret_value.append([count, comp_aln[j:j+1]])

    return(ret_value)


def alignment_with_compressed_alignment(query, target, comp_aln):

    decoded = decode_compressed_alignment(comp_aln=comp_aln)

    q = bytearray(query)
    t = bytearray(target)

    seq_index = 0

    for i in decoded:
        count = i[0]
        if i[1] == 'M':
            seq_index = seq_index + count
        elif i[1] == 'D':
            for j in range(0, count):
                t.insert(seq_index, '-')
                seq_index = seq_index + 1
        elif i[1] == 'I':
            for j in range(0, count):
                q.insert(seq_index, '-')
                seq_index = seq_index + 1

    return([q, t])


# def multiple_alignment_with_compressed_alignment(queries, target):
#     import re
#     decoded = list()

#     for c in queries:
#         d = decode_compressed_alignment(c['a'])
#         l = 0
#         for x in d:
#             l = l + x[0]

#         new_d = ''
#         for x in d:
#             new_d = new_d + x[0] * x[1]
#         splitter = re.compile(r'.')
#         new_d = splitter.findall(new_d)

#         decoded.append({'q': c['q'], 'a': new_d, 'l': l})

#     decoded.sort(key=lambda x: x['l'], reverse=True)

#     for i in range(0, len(decoded[0]['a'])):
#         longest = decoded[0]['a'][i]
#         for q in decoded[1:len(decoded)]:
#             s = set([q['a'][i], longest])
#             if len(s) > 1 and 'D' in s:
#                 q['a'][i] = ''
#         print(i)

#     return(decoded)


if __name__ == '__main__':

    pass

    # queries = ['TTCGTACGT',
    #            'TTCAAAGTAAAACGTA',
    #            'TTTCAAAGCTAAAACGTA']

    # comp_alns = ['2D1M3I6M',
    #              '2D7M3D3M1D',
    #              '3D5M1D2M3D3M1D']

    # queries = [{'q': 'TTCGTACGT', 'a': '2D1M3I6M'},
    #            {'q': 'TTCAAAGTAAAACGTA', 'a': '2D7M3D3M1D'},
    #            {'q': 'TTTCAAAGCTAAAACGTA', 'a': '3D5M1D2M3D3M1D'}]

    # target = 'CAAAGTACGT'

    # x = multiple_alignment_with_compressed_alignment(queries, target)

    # Tests

    # import os

    # PS = os.path.sep

    # to_cluster_file = 'testdata/to_cluster.fasta'
    # output_file = 'testdata/clustered.uc'

    # # cluster_file
    # cluster_file(to_cluster_file,
    #              output_file,
    #              identity_threshold=0.99,
    #              sorted_input=False,
    #              algorithm='smallmem',
    #              strand='both',
    #              threads=4,
    #              quiet=False,
    #              program='usearch6')

    # # cluster_records
    # import krbioio
    # records = krbioio.read_sequence_file(to_cluster_file, 'fasta')
    # cluster_records(records, 0.99, 'testdata')

    # # comp_aln1 = '2DM3I6M'
    # comp_aln1 = 'IDDMIIIMIMMIIIMMMI'
    # query1 = 'TTCGTACGT'
    # target1 = 'TTCAAAGTAAAACGTA'
    # aln1 = alignment_with_compressed_alignment(query1, target1, comp_aln1)
    # print(aln1[0])
    # # print(aln1[1])

    # # comp_aln2 = '2D7M3D3M'
    # comp_aln2 = 'IDDMMMMMIMMDDDMMMD'
    # query2 = 'TTCAAAGTAAAACGTA'
    # target2 = 'TTTCAAAGCTAAAACGTA'
    # aln2 = alignment_with_compressed_alignment(query2, target2, comp_aln2)
    # print(aln2[0])
    # # print(aln2[1])

    # # comp_aln3 = '3D5M1D2M3D3M'
    # comp_aln3 = 'DDDMMMMMDMMDDDMMMD'
    # query3 = 'TTTCAAAGCTAAAACGTA'
    # target3 = 'CAAAGTACGT'
    # aln3 = alignment_with_compressed_alignment(query3, target3, comp_aln3)
    # print(aln3[0])
    # print(aln3[1])
    # print()
