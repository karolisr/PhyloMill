from __future__ import print_function
# from __future__ import unicode_literals


def search(program, db, query_fasta_file, output_directory, evalue, threads, output_prefix):
    import subprocess
    import os
    output_directory = output_directory.rstrip(os.path.sep)
    output_file = output_directory + os.path.sep + output_prefix + program + '_' + db + '.xml'
    print('Running ' + program + ' on ' + db + ' database using query file "'+os.path.split(query_fasta_file)[1]+'" with E-Value of ' + str(evalue) + '.')
    subprocess.call(program+' -db '+db+' -query '+query_fasta_file+' -out '+output_file+' -evalue '+str(evalue)+' -outfmt 5 -num_threads '+str(threads), shell=True)
    return output_file


def gis_from_blast_results(file_name, min_sequence_length=-1, max_sequence_length=-1):

    from Bio.Blast import NCBIXML

    results_handle = open(file_name)
    blast_records = NCBIXML.parse(results_handle)

    gis = []

    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            if (min_sequence_length == -1 or alignment.length >= min_sequence_length) and \
               (max_sequence_length == -1 or alignment.length <= max_sequence_length):

                gi = alignment.title.split('|')[1]
                gis.append(str(gi))

    results_handle.close()
    # print("Found " + str(len(gis)) + " gis.")
    gis = list(set(gis))
    # print("There are " + str(len(gis)) + " unique gis.")
    return gis


def annotate_blast_hits(blast_results_xml, gb_records, annotation_type, qualifiers_dict=None):

    from Bio.Blast import NCBIXML
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    from krpy import krbioio

    blast_results_handle = open(blast_results_xml)
    blast_records = NCBIXML.parse(blast_results_handle)

    gb_records_dict = gb_records
    if isinstance(gb_records_dict, basestring):
        gb_records_dict = krbioio.read_sequence_file(
            file_path=gb_records,
            file_format='gb',
            ret_type='dict',
            key='gi')

    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            alignment_id = (alignment.title.split("|")[1]).split(" ")[0]
            if gb_records_dict.has_key(str(alignment_id)):
                gb_record = gb_records_dict[str(alignment_id)]
                for hsp in alignment.hsps:
                    f_start = hsp.sbjct_start-1
                    f_end = hsp.sbjct_end
                    f_strand = 1
                    if hsp.frame[1] < 0:
                        f_strand = -1
                    default_qualifiers = {'query_start':hsp.query_start, 'query_end':hsp.query_end, 'label':blast_record.query}
                    if qualifiers_dict:
                        default_qualifiers = dict(default_qualifiers.items() + qualifiers_dict.items())
                    alignment_feature = SeqFeature(FeatureLocation(f_start, f_end), strand=f_strand, type=annotation_type, qualifiers=default_qualifiers)
                    gb_record.features.append(alignment_feature)

                    gb_records_dict[str(alignment_id)] = gb_record

    blast_results_handle.close()

    return gb_records_dict


def merge_blast_hit_annotations(gb_records, annotation_type_to_merge, annotation_type_merged, merged_label, qualifiers_dict=None):

    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from krpy import krbioio
    from krpy import krother

    gb_records_dict = gb_records
    if isinstance(gb_records_dict, basestring):
        gb_records_dict = krbioio.read_sequence_file(
            file_path=gb_records,
            file_format='gb',
            ret_type='dict',
            key='gi')

    for gb_record in gb_records_dict.values():
        merged_features = []
        features = []
        for feature in gb_record.features:
            if feature.type == annotation_type_to_merge:
                features.append(feature)
        features.sort(key=lambda x: x.location.start, reverse=False)
        prev_range = [-1,-1]
        start = -1
        q_start = -1
        q_end = -1
        strand = None
        for feature in features:
            if prev_range[0] == -1:
                start = int(feature.location.nofuzzy_start)
                prev_range = [int(feature.location.nofuzzy_start),int(feature.location.nofuzzy_end)]

                q_start_temp = feature.qualifiers['query_start']
                q_end_temp = feature.qualifiers['query_end']

                q_start = None
                if not isinstance(q_start_temp, int):
                    q_start = int(feature.qualifiers['query_start'][0])
                else:
                    q_start = int(q_start_temp)

                q_end = None
                if not isinstance(q_end_temp, int):
                    q_end = int(feature.qualifiers['query_start'][0])
                else:
                    q_end = int(q_end_temp)

                if feature.strand:
                    strand = int(feature.strand)

            if not krother.in_range(int(feature.location.nofuzzy_start),prev_range[0],prev_range[1],100):
                merged_features.append([start, prev_range[1], q_start, q_end, strand])

                q_start_temp = feature.qualifiers['query_start']
                q_end_temp = feature.qualifiers['query_end']

                q_start = None
                if not isinstance(q_start_temp, int):
                    q_start = int(feature.qualifiers['query_start'][0])
                else:
                    q_start = int(q_start_temp)

                q_end = None
                if not isinstance(q_end_temp, int):
                    q_end = int(feature.qualifiers['query_start'][0])
                else:
                    q_end = int(q_end_temp)

                start = int(feature.location.nofuzzy_start)
            prev_range = [int(feature.location.nofuzzy_start),max(prev_range[1],int(feature.location.nofuzzy_end))]

            q_start_temp = feature.qualifiers['query_start']
            q_end_temp = feature.qualifiers['query_end']

            q_start = None
            if not isinstance(q_start_temp, int):
                q_start = int(feature.qualifiers['query_start'][0])
            else:
                q_start = int(q_start_temp)

            q_end = None
            if not isinstance(q_end_temp, int):
                q_end = int(feature.qualifiers['query_start'][0])
            else:
                q_end = int(q_end_temp)

            if feature.strand:
                strand = int(feature.strand)
        if len(features) > 0:
            merged_features.append([start, prev_range[1], q_start, q_end, strand])
        for merged_feature in merged_features:
            default_qualifiers = {'query_start':merged_feature[2], 'query_end':merged_feature[3], 'label':merged_label}
            if qualifiers_dict:
                default_qualifiers = dict(default_qualifiers.items() + qualifiers_dict.items())
            alignment_feature = SeqFeature(FeatureLocation(merged_feature[0], merged_feature[1]), strand=merged_feature[4], type=annotation_type_merged, qualifiers=default_qualifiers)
            gb_record.features.append(alignment_feature)
            gb_records_dict[gb_record.annotations['gi']] = gb_record
    return gb_records_dict


if __name__ == '__main__':

    # Tests

    import os

    PS = os.path.sep

    pass
