from __future__ import print_function
#from __future__ import unicode_literals

def get_taxonomy(seq_record):
    return_value = seq_record.annotations['taxonomy']
    return return_value


def get_annotation(seq_record, annotation_label):
    return_value = None
    if annotation_label in seq_record.annotations:
        return_value = seq_record.annotations[annotation_label]
    return return_value


def get_features_with_qualifier_label(record, qualifier_label,
                                      feature_type=None):
    features = []
    for i, feature in enumerate(record.features):
        if ((feature_type is not None and feature.type == feature_type) or
                (feature_type is None)):
            for qualifier_key in feature.qualifiers.keys():
                if qualifier_key == qualifier_label:
                    features.append(i)
    return features


def get_features_with_qualifier(record, qualifier_label, qualifier,
                                feature_type=None, loose=False):
    result_features = []
    features = get_features_with_qualifier_label(record, qualifier_label,
                                                 feature_type)
    for feature_index in features:
        feature = record.features[feature_index]
        for q in feature.qualifiers[qualifier_label]:
            #if not loose and (q.lower().startswith(qualifier.lower())):
            if not loose and (q.lower() == qualifier.lower()):
                result_features.append(feature_index)
            if loose and qualifier.lower() in q.lower():
                result_features.append(feature_index)
    return result_features


def trim_residues(bio_seq_record, trim_length, right=False):
    '''
        Will trim trim_length residues from the left or right side of the
        sequence record.
    '''
    from Bio import SeqRecord
    sequence_record = bio_seq_record.seq.tomutable()
    keep_start = trim_length
    keep_end = len(sequence_record)
    if right:
        keep_start = 0
        keep_end = len(sequence_record) - trim_length
    sequence_record = sequence_record[keep_start:keep_end]
    letter_annotations = dict()
    if 'phred_quality' in bio_seq_record.letter_annotations:
        quality_scores = bio_seq_record.letter_annotations['phred_quality']
        letter_annotations['phred_quality'] = (
            quality_scores[keep_start:keep_end])
    bio_seq_record = SeqRecord.SeqRecord(
        seq=sequence_record.toseq(),
        id=bio_seq_record.id,
        name=bio_seq_record.name,
        description=bio_seq_record.description,
        dbxrefs=None,
        features=None,
        annotations=None,
        letter_annotations=letter_annotations)
    return bio_seq_record


def reverse_complement(bio_seq_record):
    result = bio_seq_record.reverse_complement(id=True, name=True,
                                               description=True,
                                               features=True, annotations=True,
                                               letter_annotations=True,
                                               dbxrefs=True)
    return result


def resolve_ambiguities(sequence):
    import numpy
    import kriupac
    for k in kriupac.IUPAC_AMBIGUOUS_DNA_DICT.keys():
        for i in range(0, sequence.count(kriupac.IUPAC_AMBIGUOUS_DNA_DICT[k])):
            rand = numpy.random.randint(0, len(k))
            sequence = sequence.replace(kriupac.IUPAC_AMBIGUOUS_DNA_DICT[k], k[rand], 1)
    return(sequence)


def location_from_string(location_string):

    from Bio.SeqFeature import CompoundLocation
    from Bio.SeqFeature import FeatureLocation
    from Bio.SeqFeature import BeforePosition, AfterPosition, ExactPosition

    operator_parts = location_string.split('{')
    location_strings = None
    operator = None
    if len(operator_parts) > 1:
        operator = operator_parts[0]
        location_strings = operator_parts[1].strip('}').split(', ')
    else:
        location_strings = operator_parts

    locations = list()

    for loc_str in location_strings:

        loc_str_split = loc_str.split('(')
        strand = None
        if len(loc_str_split) > 1:
            strand = int(loc_str_split[1].strip(')') + '1')
        bounds_list = loc_str_split[0].strip('[').strip(']').split(':')

        ########################################################################
        # Hack. Fixes things like: AY148051.1[0:362] in
        #   order{AY148051.1[0:362](+), [0:>355](+)}
        should_skip = False
        for b in bounds_list:
            if '[' in b:
                should_skip = True
                break
        if should_skip:
            continue
        ########################################################################

        l_bound = bounds_list[0].split('<')
        start = None
        if len(l_bound) == 2:
            start = BeforePosition(int(l_bound[1]))
        else:
            l_bound = bounds_list[0].split('>')
            if len(l_bound) == 2:
                start = AfterPosition(int(l_bound[1]))
            else:
                start = ExactPosition(int(l_bound[0]))

        r_bound = bounds_list[1].split('>')
        end = None
        if len(r_bound) == 2:
            end = AfterPosition(int(r_bound[1]))
        else:
            r_bound = bounds_list[1].split('<')
            if len(r_bound) == 2:
                end = BeforePosition(int(r_bound[1]))
            else:
                end = ExactPosition(int(r_bound[0]))

        loc = FeatureLocation(start=start, end=end, strand=strand)

        locations.append(loc)

    location = None
    if len(locations) > 1:
        location = CompoundLocation(parts=locations, operator=operator)
    else:
        location = locations[0]

    return location


def translate_cds(record, table):

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import generic_protein

    # Translation tables
    # http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi

    # Extract CDS
    feature_list = list()

    for f in record.features:
        if f.type.lower() == 'cds':
            feature_list.append(f)

    extraction_list = list()

    for f in feature_list:
        e = f.extract(record)
        start = int(f.qualifiers['codon_start'][0])
        if start > 1:
            e = e[start-1:len(e)]
        extraction_list.append(e)

    translation_list = list()

    for e in extraction_list:
        t = e.seq.translate(table=table)
        translation_list.append(t)

    seq = ''

    for t in translation_list:
        seq = seq + str(t)

    seq = Seq(seq, generic_protein)
    rec = SeqRecord(seq)

    rec.name = record.id.split('.')[0]
    rec.description = record.description
    rec.annotations['gi'] = record.annotations['gi']
    rec.annotations['organism'] = record.annotations['organism']
    rec.annotations['taxonomy'] = record.annotations['taxonomy']
    rec.id = record.id

    return rec


def merge_record_features(gb_records, annotation_type_to_merge, annotation_type_merged, merged_label, qualifiers_dict=None):

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
        # q_start = -1
        # q_end = -1
        strand = None
        for feature in features:
            if prev_range[0] == -1:
                start = int(feature.location.nofuzzy_start)
                prev_range = [int(feature.location.nofuzzy_start),int(feature.location.nofuzzy_end)]

                # q_start_temp = feature.qualifiers['query_start']
                # q_end_temp = feature.qualifiers['query_end']

                # q_start = None
                # if not isinstance(q_start_temp, int):
                #     q_start = int(feature.qualifiers['query_start'][0])
                # else:
                #     q_start = int(q_start_temp)

                # q_end = None
                # if not isinstance(q_end_temp, int):
                #     q_end = int(feature.qualifiers['query_start'][0])
                # else:
                #     q_end = int(q_end_temp)

                if feature.strand:
                    strand = int(feature.strand)

            if not krother.in_range(int(feature.location.nofuzzy_start),prev_range[0],prev_range[1],100):
                merged_features.append([start, prev_range[1], strand])

                # q_start_temp = feature.qualifiers['query_start']
                # q_end_temp = feature.qualifiers['query_end']

                # q_start = None
                # if not isinstance(q_start_temp, int):
                #     q_start = int(feature.qualifiers['query_start'][0])
                # else:
                #     q_start = int(q_start_temp)

                # q_end = None
                # if not isinstance(q_end_temp, int):
                #     q_end = int(feature.qualifiers['query_start'][0])
                # else:
                #     q_end = int(q_end_temp)

                start = int(feature.location.nofuzzy_start)
            prev_range = [int(feature.location.nofuzzy_start),max(prev_range[1],int(feature.location.nofuzzy_end))]

            # q_start_temp = feature.qualifiers['query_start']
            # q_end_temp = feature.qualifiers['query_end']

            # q_start = None
            # if not isinstance(q_start_temp, int):
            #     q_start = int(feature.qualifiers['query_start'][0])
            # else:
            #     q_start = int(q_start_temp)

            # q_end = None
            # if not isinstance(q_end_temp, int):
            #     q_end = int(feature.qualifiers['query_start'][0])
            # else:
            #     q_end = int(q_end_temp)

            if feature.strand:
                strand = int(feature.strand)
        if len(features) > 0:
            merged_features.append([start, prev_range[1], strand])
        for merged_feature in merged_features:
            # default_qualifiers = {'query_start':merged_feature[2], 'query_end':merged_feature[3], 'label':merged_label}
            # if qualifiers_dict:
            #     default_qualifiers = dict(default_qualifiers.items() + qualifiers_dict.items())
            # alignment_feature = SeqFeature(FeatureLocation(merged_feature[0], merged_feature[1]), strand=merged_feature[2], type=annotation_type_merged)
            alignment_feature = SeqFeature(FeatureLocation(merged_feature[0], merged_feature[1]), strand=1, type=annotation_type_merged)
            gb_record.features.append(alignment_feature)
            gb_records_dict[gb_record.annotations['gi']] = gb_record
    return gb_records_dict


# if __name__ == '__main__':

    # Tests

    # pass

    # from krpy import krbioio

    # records = krbioio.read_sequence_file(
    #     file_path='/Users/karolis/Desktop/strange_record.gb',
    #     file_format='genbank',
    #     ret_type='list')

    # record = records[0]
    # table = 1
    # extraction_list = translate_cds(record, table)

    # krbioio.write_sequence_file(
    #     records=extraction_list,
    #     file_path='/Users/karolis/Desktop/p_virginiana_e.gb',
    #     file_format='genbank')

    # location_from_string('join{[<109:199](+), [294:358](+), [444:545](+), [635:745](+), [829:>886](+)}')
    # location_from_string('[294:358](+)')

    # feats = records[0].features

    # for f in feats:
    #     # print(str(f.location))
    #     print(location_from_string(str(f.location)))
    #     print('--- --- ---')
