# -*- coding: utf-8 -*-

"""

This module deals with reading and writing of biological sequence and
alignment files.

"""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import with_statement

from xml.etree import ElementTree

from krpy import Error


def _parse_gbseq_xml_handle(gbseq_xml_handle):

    """

    Parses GBSeq XML handle.

    :param gbseq_xml_handle: A handle to parse.
    :type gbseq_xml_handle: handle

    :returns: A dictionary with these keys: accession, date_create,
        date_update, definition, division, features, length,
        mol_type, organism, seq, strandedness, taxid, taxonomy,
        topology, version
    :rtype: dict

    """

    tree = ElementTree.parse(gbseq_xml_handle)
    root = tree.getroot()

    return_value = list()

    for n_record in root.findall('GBSeq'):

        accession = None
        version = None
        # gi = None

        # n_seqids = n_record.find('GBSeq_other-seqids')
        # for n_seqid in n_seqids.findall('GBSeqid'):
        #     if n_seqid.text.startswith('gi|'):
        #         gi = n_seqid.text.strip('gi|')

        temp_acc_ver = n_record.find('GBSeq_accession-version')
        temp_acc_ver = temp_acc_ver.text
        accession = temp_acc_ver.split('.')[0]
        version = int(temp_acc_ver.split('.')[1])

        strandedness = n_record.find('GBSeq_strandedness')
        if strandedness is not None:
            strandedness = strandedness.text

        date_create = n_record.find('GBSeq_create-date')
        date_create = date_create.text

        date_update = n_record.find('GBSeq_update-date')
        date_update = date_update.text

        division = n_record.find('GBSeq_division')
        division = division.text

        topology = n_record.find('GBSeq_topology')
        topology = topology.text

        definition = n_record.find('GBSeq_definition')
        definition = definition.text

        mol_type = n_record.find('GBSeq_moltype')
        mol_type = mol_type.text

        organism = n_record.find('GBSeq_organism')
        organism = organism.text

        taxonomy = n_record.find('GBSeq_taxonomy')
        taxonomy = taxonomy.text.split('; ')

        seq = n_record.find('GBSeq_sequence')
        seq = seq.text

        length = n_record.find('GBSeq_length')
        length = int(length.text)

        if length != len(seq):
            message = (
                'Reported sequence length does not match the actual '
                'length of sequence. Reported: {r}; Actual: {a}.')
            message = message.format(r=length, a=len(seq))
            raise Error(message)

        features = dict()
        n_feature_table = n_record.find('GBSeq_feature-table')
        for n_feature in n_feature_table.findall('GBFeature'):
            n_feature_key = n_feature.find('GBFeature_key')
            fk = n_feature_key.text
            features[fk] = dict()

            features[fk]['intervals'] = list()
            n_intervals = n_feature.find('GBFeature_intervals')
            for n_interval in n_intervals.findall('GBInterval'):

                n_interval_from = n_interval.find('GBInterval_from')
                n_interval_to = n_interval.find('GBInterval_to')
                n_interval_point = n_interval.find('GBInterval_point')

                interval = None

                if (n_interval_from is not None) and \
                   (n_interval_to is not None):

                    start = int(n_interval_from.text)
                    end = int(n_interval_to.text)

                    # Make intervals zero-indexed.
                    if start < end:
                        start = start - 1
                    elif end < start:
                        end = end - 1

                    interval = [start, end]

                elif n_interval_point is not None:

                    interval = [int(n_interval_point.text) - 1]

                features[fk]['intervals'].append(interval)

            features[fk]['qualifiers'] = dict()
            qualifiers = features[fk]['qualifiers']
            n_qualifiers = n_feature.find('GBFeature_quals')
            for n_qualifier in n_qualifiers.findall('GBQualifier'):
                n_qualifier_name = n_qualifier.find('GBQualifier_name')
                n_qualifier_value = n_qualifier.find('GBQualifier_value')
                qualifiers[n_qualifier_name.text] = n_qualifier_value.text

        taxid = None
        taxid_temp = features['source']['qualifiers']['db_xref']
        if taxid_temp.startswith('taxon:'):
            taxid = taxid_temp.split('taxon:')[1]

        record_dict = dict()

        record_dict['seq'] = seq
        record_dict['mol_type'] = mol_type

        record_dict['accession'] = accession
        record_dict['version'] = version
        # record_dict['gi'] = gi
        record_dict['definition'] = definition

        record_dict['strandedness'] = strandedness
        record_dict['topology'] = topology
        record_dict['division'] = division

        record_dict['date_create'] = date_create
        record_dict['date_update'] = date_update

        record_dict['taxid'] = taxid
        record_dict['organism'] = organism
        record_dict['taxonomy'] = taxonomy

        record_dict['features'] = features

        record_dict['length'] = length

        return_value.append(record_dict)

    return return_value
