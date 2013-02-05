from __future__ import print_function
#from __future__ import unicode_literals

def get_annotation(seq_record, annotation_label):
    return_value = None
    if seq_record.annotations.has_key(annotation_label):
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
    feature_type=None):
    result_features = []
    features = get_features_with_qualifier_label(record, qualifier_label,
        feature_type)
    for feature_index in features:
        feature = record.features[feature_index]
        for q in feature.qualifiers[qualifier_label]:
            if q.startswith(qualifier):
                result_features.append(feature_index)
    return result_features
