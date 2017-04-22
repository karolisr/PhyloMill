# -*- coding: utf-8 -*-

"""

This module contains classes that deal with amino acid and nucleotide
sequences.

"""

from __future__ import absolute_import
from __future__ import division
from __future__ import generators
from __future__ import nested_scopes
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import with_statement

from datetime import date

from krpy import Root
from krpy import Error
from krpy import STRING_TYPE
from krpy.iupac import *

SEQ_TYPE_NT = 'NT'
SEQ_TYPE_DNA = 'DNA'
SEQ_TYPE_RNA = 'RNA'
SEQ_TYPE_AA = 'AA'

SEQ_TYPES = [SEQ_TYPE_NT, SEQ_TYPE_DNA, SEQ_TYPE_RNA, SEQ_TYPE_AA]

MOL_TO_SEQ_TYPE_MAP = {
    'NT': SEQ_TYPE_NT,
    'DNA': SEQ_TYPE_DNA,
    'RNA': SEQ_TYPE_RNA,
    'mRNA': SEQ_TYPE_RNA,
    'rRNA': SEQ_TYPE_RNA,
    'tRNA': SEQ_TYPE_RNA,
    'AA': SEQ_TYPE_AA}


class Seq(Root):

    """

    This class can be used to create new amino acid and nucleotide
    sequence objects of types: :class:`NTSeq`,
    :class:`DNASeq`, :class:`RNASeq`, :class:`AASeq`.

    :param seq: Nucleotide or amino acid sequence.
    :type seq: str

    :param seq_type: The type of the sequence. One of the constants
        defined in module :mod:`krpy.seq` should be used:
        ``SEQ_TYPE_NT``, ``SEQ_TYPE_DNA``, ``SEQ_TYPE_RNA``,
        ``SEQ_TYPE_AA``. These constants have string values ``NT``,
        ``DNA``, ``RNA``, ``AA``, respectively.
    :type seq_type: str

    """

    def __new__(self, seq, seq_type):

        if seq_type in SEQ_TYPES:

            seq = seq.upper()

            if seq_type is SEQ_TYPE_NT:
                return NTSeq(seq=seq)
            elif seq_type is SEQ_TYPE_DNA:
                return DNASeq(seq=seq)
            elif seq_type is SEQ_TYPE_RNA:
                return RNASeq(seq=seq)
            elif seq_type is SEQ_TYPE_AA:
                return AASeq(seq=seq)

        else:
            # ToDo: Report error
            pass

    # __init__ declaration is exactly the same as __new__ so Sphinx
    # docstring parser picks it up.
    def __init__(self, seq, seq_type):
        pass


class _Seq(Root):

    """

    This is an abstract class.

    """

    def __init__(self, seq):
        self._seq = seq

    @property
    def seq(self):
        return self._seq

    @property
    def length(self):
        return len(self.seq)


class NTSeq(_Seq):

    """

    This class represents any nucleotide sequence and does not
    distinguish between RNA and DNA. Both thymine and uracil can occur
    in the instances of this class. :class:`DNASeq` and :class:`RNASeq`
    classes inherit from this class.

    :param seq: Nucleotide sequence.
    :type seq: str

    """

    def __init__(self, seq):
        if set(seq) <= NT_AMBIGUOUS:
            super(NTSeq, self).__init__(seq)
        else:
            message = ('NT sequence should not contain these characters: {s}.')
            message = message.format(s=', '.join(str(s) for s in set(seq) - NT_AMBIGUOUS))
            raise Error(message)


class DNASeq(NTSeq):

    """

    This class represents DNA sequence.

    :param seq: DNA sequence.
    :type seq: str

    """

    def __init__(self, seq):
        if set(seq) <= DNA_AMBIGUOUS:
            super(DNASeq, self).__init__(seq)
        elif set(seq) - DNA_AMBIGUOUS == RNA_ONLY_CHARS:
            seq = seq.replace('U', 'T')
            super(DNASeq, self).__init__(seq)
        else:
            message = ('DNA sequence should not contain these characters: {s}.')
            message = message.format(s=', '.join(str(s) for s in set(seq) - DNA_AMBIGUOUS))
            raise Error(message)


class RNASeq(NTSeq):

    """

    This class represents RNA sequence.

    :param seq: RNA sequence.
    :type seq: str

    """

    def __init__(self, seq):
        if set(seq) <= RNA_AMBIGUOUS:
            super(RNASeq, self).__init__(seq)
        elif set(seq) - RNA_AMBIGUOUS == DNA_ONLY_CHARS:
            seq = seq.replace('T', 'U')
            super(RNASeq, self).__init__(seq)
        else:
            message = ('RNA sequence should not contain these characters: {s}.')
            message = message.format(s=', '.join(str(s) for s in set(seq) - RNA_AMBIGUOUS))
            raise Error(message)


class AASeq(_Seq):

    """

    This class represents amino acid sequence.

    :param seq: Amino acid sequence.
    :type seq: str

    """

    def __init__(self, seq):
        if set(seq) <= AA_AMBIGUOUS:
            super(AASeq, self).__init__(seq)
        else:
            message = ('AA sequence should not contain these characters: {s}.')
            message = message.format(s=', '.join(str(s) for s in set(seq) - AA_AMBIGUOUS))
            raise Error(message)


class SeqRecord(Root):

    """

    This class stores a sequence and its associated meta-information. It
    is designed to accomodate GenBank records. For reference see:
    http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html

    """

    def __init__(self,
                 seq,
                 mol_type,
                 accession=None,
                 version=None,
                 # gi=None,
                 description=None,
                 strandedness=None,
                 topology=None,
                 division=None,
                 date_create=None,
                 date_update=None,
                 taxid=None,
                 organism=None,
                 taxonomy=None,
                 features=None):

        # immutable setters
        self._set_mol_type(mol_type)
        self._set_seq_type(self.mol_type)
        self._set_seq(seq)

        # mutable setters
        self.accession = accession
        self.version = version
        # self.gi = gi
        self.description = description

        self.strandedness = strandedness
        self.topology = topology
        self.division = division

        self.date_create = date_create
        self.date_update = date_update

        self.taxid = taxid
        self.organism = organism
        self.taxonomy = taxonomy

        self.features = features

    # seq_type and mol_type are not the same. mol_type is more specific.
    # For example: mol_type mRNA, seq_type RNA

    # mol_type
    @property
    def mol_type(self):
        return self._mol_type

    def _set_mol_type(self, value):
        self._mol_type = value

    # seq_type
    def _set_seq_type(self, value):

        if value in MOL_TO_SEQ_TYPE_MAP:
            self._seq_type = MOL_TO_SEQ_TYPE_MAP[value]
        elif value in SEQ_TYPES:
            self._seq_type = value
        else:
            message = ('Molecule type not supported: {s}.')
            message = message.format(s=value)
            raise Error(message)

    # seq
    @property
    def seq(self):
        return self._seq

    def _set_seq(self, value):

        if issubclass(type(value), _Seq):
            self._seq = value
        elif issubclass(type(value), STRING_TYPE):
            self._seq = Seq(seq=value, seq_type=self._seq_type)

    # accession
    @property
    def accession(self):
        return self._accession

    @accession.setter
    def accession(self, value):
        self._accession = value

    # version
    @property
    def version(self):
        return self._version

    @version.setter
    def version(self, value):
        if value is not None:
            try:
                self._version = int(value)
            except ValueError:
                print('Version should be an integer.')

    # # gi
    # @property
    # def gi(self):
    #     return self._gi

    # @gi.setter
    # def gi(self, value):
    #     try:
    #         self._gi = int(value)
    #     except ValueError:
    #         print('GI should be an integer.')

    # description (definition)
    @property
    def description(self):
        return self._description

    @description.setter
    def description(self, value):
        self._description = value

    # strandedness
    @property
    def strandedness(self):
        return self._strandedness

    @strandedness.setter
    def strandedness(self, value):
        self._strandedness = value

    # topology
    @property
    def topology(self):
        return self._topology

    @topology.setter
    def topology(self, value):
        self._topology = value

    # division
    @property
    def division(self):
        return self._division

    @division.setter
    def division(self, value):
        self._division = value

    # Deal with dates
    def _process_date_string(self, date_str):
        return_value = None
        split_date_str = date_str.split('-')
        if len(split_date_str) != 3:
            raise Error('Date should be a string of format: YYYY-MM-DD')
        else:
            try:
                return_value = date(
                    int(split_date_str[0]),
                    int(split_date_str[1]),
                    int(split_date_str[2]))
            except Exception as e:
                raise e

        return return_value

    # date_create
    @property
    def date_create(self):
        return self._date_create

    @date_create.setter
    def date_create(self, value):
        if value is not None:
            self._date_create = self._process_date_string(date_str=value)

    # date_update
    @property
    def date_update(self):
        return self._date_update

    @date_update.setter
    def date_update(self, value):
        if value is not None:
            self._date_update = self._process_date_string(date_str=value)

    # taxid
    @property
    def taxid(self):
        return self._taxid

    @taxid.setter
    def taxid(self, value):
        if value is not None:
            try:
                self._taxid = int(value)
            except ValueError:
                print('taxid should be an integer.')

    # organism
    @property
    def organism(self):
        return self._organism

    @organism.setter
    def organism(self, value):
        self._organism = value

    # taxonomy
    @property
    def taxonomy(self):
        return self._taxonomy

    @taxonomy.setter
    def taxonomy(self, value):
        self._taxonomy = value

    # features
    @property
    def features(self):
        return self._features

    @features.setter
    def features(self, value):
        self._features = value
