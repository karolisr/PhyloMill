from __future__ import print_function
#from __future__ import unicode_literals


def get_annotation(seq_record, annotation_label):
    return_value = None
    if annotation_label in seq_record.annotations:
        return_value = seq_record.annotations[annotation_label]
    return return_value


def get_features_with_qualifier_label(record, qualifier_label,
    feature_type=None):
    features = []
    for i, feature in enumerate(record.features):
        if ((feature_type != None and feature.type == feature_type) or
            (feature_type == None)):
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
            if not loose and q.startswith(qualifier):
                result_features.append(feature_index)
            if loose and qualifier in q:
                result_features.append(feature_index)
    return result_features
