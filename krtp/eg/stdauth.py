#!/usr/bin/python

import sys
import re
import codecs
from unidecode import unidecode

# Usage: std-auth.py auth.in > auth.out

# auth.in and auth.out: text files with one line per Authority string

# Also provide: authtable
#     Format = desired name | alternative 1 | alt 2 etc., e.g.:
#     Hammerst. | Hammersteinxx | Hammerstein
authtable = 'authority_alternates.dat'

# The input author names file can have spaces and punctuation within a name.
# These will be removed before forming the dictionary, however.

# Dealing with extended character sets and unicode (for, e.g., accented
# letters) turned out to be a big pain.  The shell and R commands to substitute
# relevant ascii characters are iconv with ASCII//TRANSLIT.  But in python, I
# had the best (only) success with a contributed library called Unidecode.
# http://pypi.python.org/pypi/Unidecode/
#
# Anyway, the overall plan is to allow non-ascii characters in the input files
# (author table and names to be standardized) and to return ascii output.

def clean_name(name):
    '''
        For the name of a single author.
        Replace special, non-ascii characters.
        Get rid of all spaces and the specified punctuation.
    '''
    name = unidecode(name)
    return (name.translate(None, " ,.-"))


def make_auth_dic(filename):
    '''
        Create the author table dictionary from the user-provided file.
        Note that 'desired names' get included as an 'alternative'.
    '''

    fp = codecs.open(filename, 'r', encoding='utf-8')
    dic = {}

    for line in fp:
        line = line.strip()
        if line and not line.startswith('#'):

            line_list = line.split('|')
            key = unidecode(line_list[0].strip())
            val = [ clean_name(x) for x in line_list ]

            if (dic.has_key(key)):
                dic[key] = dic[key] + val
            else:
                dic[key] = val

    fp.close()
    return(dic)


def get_std_auth(name, dic):
    '''
        Look in the author table dictionary.
        Return a standardized name if found.
        Otherwise, just return the name in ascii (without any accents).
    '''

    name2 = clean_name(name)

    auth = [k for k, v in dic.iteritems() if name2 in v]

    if len(auth) == 1:
        ans = auth[0]
        # print(name + "\tmatched to\t" + ans)

    else:
        ans = unidecode(name)
        # if len(auth) == 0:
        #     print(name + "\tunmatched")
        # else:
        #     print(name + "\tmultiply matched")

    return(ans)


def clean_spaces(line):
    '''
        Ensure:
            no space before .
            no space after . unless followed by ex, in, &, (
            no space before and one space after ,
            no space after ( or before )
            one space before ( and after )
            one space before and after &
            no multiple sequential spaces
            no spaces at beginning or end
    '''

    # could surely be more efficient with re.sub...
    line = '.'.join([x.strip() for x in line.split('.')])
    line = '. ex '.join([x.strip() for x in line.split('.ex ')])
    line = '. in '.join([x.strip() for x in line.split('.in ')])
    line = ' al. '.join([x.strip() for x in line.split(' al.')])
    line = ', '.join([x.strip() for x in line.split(',')])
    line = ') '.join([x.strip() for x in line.split(')')])
    line = ' ('.join([x.strip() for x in line.split('(')])
    line = ' & '.join([x.strip() for x in line.split('&')])
    line = re.sub('^\"\s*|\s*\"$', '"', line) # no space by surrounding quotes
    line = ' '.join([x.strip() for x in line.split()])
    return(line.strip())


def translate(authority, authority_alternates):
    '''
        This will translate a given authority name to a preffered spelling.

        Example:
            authority_alternates = make_auth_dic(authtable)
            translate('Linnaeus', authority_alternates)

        Note: (Karolis R.) I refactored this code from
            "if __name__ == '__main__'" block to make this file usable as a
            module, instead of supplying an input file, now we can translate
            individual authority alternates. I made adjustments so the original
            functionality remains unchanged.
    '''
    separators = ['', '(', ')', '&', ' ex ', ' in ', ', ', ' et al.', '"']
    splitter = re.compile('(\(|\)|\&| ex | in |\, | et al.|\")')

    authority = authority.replace(' Ex ', ' ex ')
    authority = authority.replace(' and ', ' & ')
    authority = authority.replace(' et ', ' & ')         # breaks 'et al.'
    authority = authority.replace(' & al.', ' et al.') # must follow 'et' sub

    line_list = splitter.split(authority)
    for (p, piece) in enumerate(line_list):

        if piece not in separators:

            # Set aside any leading and trailing whitespace.
            # The one useful element is the name to match.
            name_list = re.split('(^\s*|\s*$)', piece)
            for (i, name) in enumerate(name_list):
                if not (name.isspace() or name == ''):
                    name_list[i] = get_std_auth(name, authority_alternates)
                    # print(name_list[i])

            # Restore the leading/trailing whitespace to this name.
            line_list[p] = ''.join(name_list)

    # Join the pieces of this line and clean up the spacing.
    return clean_spaces(''.join(line_list))

if __name__ == '__main__':

    if (len(sys.argv) != 2):
        print('Usage: %s infile' % sys.argv[0])
        sys.exit(1)

    # A useful trick: /dev/stdin substitutes for an input file
    #     ./std-auth.py /dev/stdin < infile
    #     cat infile | ./std-auth.py /dev/stdin

    # The name standardization table.
    authtable = re.sub('stdauth.py$', authtable, sys.argv[0])
    authority_alternates = make_auth_dic(authtable)

    # The list of names to be standardized.
    fp = codecs.open(sys.argv[1], 'r', encoding='utf-8')

    for line in fp:

        print(translate(line, authority_alternates))

    fp.close()
