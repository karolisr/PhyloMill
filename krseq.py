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
