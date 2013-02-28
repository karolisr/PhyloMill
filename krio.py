from __future__ import print_function
#from __future__ import unicode_literals


def prepare_directory(path):

    '''
    Checks if directory at path exists and, if not, creates full path.
    '''

    import os
    if not os.path.exists(path):
        os.makedirs(path)
    return

# The commented file reading class from:
# http://www.mfasold.net/blog/2010/02/
#   python-recipe-read-csvtsv-textfiles-and-ignore-comment-lines


def num_lines_in_file(file_path):
    with open(file_path) as f:
        for i, l in enumerate(f):
            pass
        # ah... wonderful Python scope rules!
    return(i + 1)


class CommentedFile:

    '''
    Provide an open file handle with comments removed.
    '''

    def __init__(self, f, commentstring="#"):
        self.f = f
        self.commentstring = commentstring

    def next(self):
        line = self.f.next()
        while line.startswith(self.commentstring) or not line.strip():
            line = self.f.next()
        return line

    def __iter__(self):
        return self

    def close(self):
        self.f.close()


def read_table_file(
    path,
    has_headers=False,
    headers=None,
    delimiter=',',
    quotechar='"',
    stripchar='',
    rettype='dict'  # dict, list, set
):

    '''
    Reads a delimited text file.
    Returns:
        A list of dictionaries, one dictionary per row. Header names as keys.
    '''

    import krother

    handle = CommentedFile(open(path, 'rb'))
    if headers is not None:
        has_headers = False
    if has_headers:
        headers = handle.next()
        headers = krother.parse_line(headers, delimiter, quotechar, stripchar)

    return_value = list()

    if rettype.startswith('dict'):
        for l in handle:

            l_spl = krother.parse_line(l, delimiter, quotechar, stripchar)

            row_dict = dict()
            for i, h in enumerate(headers):
                row_dict[h] = l_spl[i]
            return_value.append(row_dict)

    if rettype.startswith('list') or rettype.startswith('set'):
        for l in handle:

            l_spl = krother.parse_line(l, delimiter, quotechar, stripchar)

            if rettype.startswith('set') and len(l_spl) == 1:
                l_spl = l_spl[0]

            return_value.append(l_spl)

    if rettype.startswith('set'):
        return_value = set(return_value)

    return(return_value)

if __name__ == '__main__':

    # Tests

    import os

    ps = os.path.sep

    # CommentedFile
    handle = CommentedFile(open('testdata' + ps + 'commented_file.csv', 'rb'))
    for i, line in enumerate(handle):
        print(i, repr(line))
    handle.close()

    # read_table_file
    table = read_table_file(path='testdata' + ps + 'commented_file.csv',
                            has_headers=True, headers=None, delimiter=',')
    for i, line in enumerate(table):
        print(i, repr(line))
