"""
IUPAC
Nucleotide Code:  Base:
----------------  -----
A.................Adenine
C.................Cytosine
G.................Guanine
T (or U)..........Thymine (or Uracil)
R.................A or G
Y.................C or T
S.................G or C
W.................A or T
K.................G or T
M.................A or C
B.................C or G or T
D.................A or G or T
H.................A or C or T
V.................A or C or G
N.................any base
. or -............gap
"""

IUPAC_AMBIGUOUS_DNA_STRING = 'RYSWKMBDHVN.-'

IUPAC_DNA_STRING = 'ACGT'

IUPAC_DNA_CHARACTERS = set(
    ['A',
     'C',
     'G',
     'T',
     'U',
     'N',
     '.',
     '-',
     'R',
     'Y',
     'S',
     'W',
     'K',
     'M',
     'B',
     'D',
     'H',
     'V'])

IUPAC_DNA_GAPS_STRING = '.-'

IUPAC_DNA_GAPS = set(
    ['.',
     '-'])

IUPAC_DNA_UNKNOWN = set(
    ['N'])

IUPAC_DOUBLE_DNA_DICT = {
    'AG': 'R',
    'CT': 'Y',
    'AC': 'M',
    'GT': 'K',
    'AT': 'W',
    'CG': 'S'
}

IUPAC_TRIPLE_DNA_DICT = {
    'CGT': 'B',
    'AGT': 'D',
    'ACT': 'H',
    'ACG': 'V'
}

IUPAC_AMBIGUOUS_DNA_DICT = {
    'AG': 'R',
    'CT': 'Y',
    'AC': 'M',
    'GT': 'K',
    'AT': 'W',
    'CG': 'S',
    'CGT': 'B',
    'AGT': 'D',
    'ACT': 'H',
    'ACG': 'V'
}

IUPAC_DNA_DICT = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'T',
    'AG': 'R',
    'CT': 'Y',
    'AC': 'M',
    'GT': 'K',
    'AT': 'W',
    'CG': 'S',
    'CGT': 'B',
    'AGT': 'D',
    'ACT': 'H',
    'ACG': 'V'
}

IUPAC_DOUBLE_DNA_DICT_REVERSE = {
    'R': 'AG',
    'Y': 'CT',
    'M': 'AC',
    'K': 'GT',
    'W': 'AT',
    'S': 'CG'
}

IUPAC_TRIPLE_DNA_DICT_REVERSE = {
    'B': 'CGT',
    'D': 'AGT',
    'H': 'ACT',
    'V': 'ACG'
}

IUPAC_AMBIGUOUS_DNA_DICT_REVERSE = {
    'R': 'AG',
    'Y': 'CT',
    'M': 'AC',
    'K': 'GT',
    'W': 'AT',
    'S': 'CG',
    'B': 'CGT',
    'D': 'AGT',
    'H': 'ACT',
    'V': 'ACG'
}

IUPAC_DNA_DICT_REVERSE = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'T',
    'R': 'AG',
    'Y': 'CT',
    'M': 'AC',
    'K': 'GT',
    'W': 'AT',
    'S': 'CG',
    'B': 'CGT',
    'D': 'AGT',
    'H': 'ACT',
    'V': 'ACG'
}
