from __future__ import print_function
#from __future__ import unicode_literals


def parse_directory(path, file_name_sep, sort='forward'):

    '''
    Will parse a directory at a given path and return a list of dictionary
    objects with keys:
        name: name of a file without extension
        ext: file extension
        full: file name with extension
        path: full file path, relative to the input path
        split: file name split using file_name_sep input variable
    '''

    import os

    ps = os.path.sep
    path = path.rstrip(ps) + ps
    file_list = os.listdir(path)
    file_list.sort(reverse=False)
    if sort == 'reverse':
        file_list.sort(reverse=True)
    return_list = list()

    for f in file_list:
        # Mac hack
        if f == '.DS_Store':
            continue
        file_name = os.path.splitext(f)[0]
        file_ext = None
        if os.path.splitext(f)[1] != '':
            file_ext = os.path.splitext(f)[1].split('.')[1]
        isdir = os.path.isdir(path + f)
        file_name_split = file_name.split(file_name_sep)

        file_dict = dict()
        file_dict['name'] = file_name
        file_dict['ext'] = file_ext
        file_dict['full'] = f
        file_dict['path'] = path + f
        file_dict['basepath'] = path
        file_dict['split'] = file_name_split
        file_dict['isdir'] = isdir
        return_list.append(file_dict)

    return return_list


def prepare_directory(path):

    '''
    Checks if directory at path exists and, if not, creates full path.
    '''

    import os
    ps = os.path.sep
    if not os.path.exists(path):
        os.makedirs(path)
    path = path.rstrip(ps) + ps
    return(path)

# The commented file reading class from:
# http://www.mfasold.net/blog/2010/02/
#   python-recipe-read-csvtsv-textfiles-and-ignore-comment-lines


def num_lines_in_file(file_path, print_every=None):
    with open(file_path) as f:
        for i, l in enumerate(f):
            pass
            # if print_every:
                # if i % print_every == 0:
                    # print(i)
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
    handle=None,
    path=None,
    has_headers=False,
    headers=None,
    delimiter=',',
    quotechar='"',
    stripchar='',
    commentchar="#",
    rettype='dict'  # dict, list, set
):

    '''
    Reads a delimited text file.
    Returns:
        A list of dictionaries, one dictionary per row. Header names as keys.
    '''

    import os

    import krother

    if not os.path.exists(path):
        return None

    if not handle:
        handle = CommentedFile(open(path, 'rU'), commentstring=commentchar)
    else:
        handle = CommentedFile(handle, commentstring=commentchar)

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


# ToDo: read_tsv_file and write_tsv_file are special functions and may need to
# be renamed to a less generic name

def read_tsv_file(l_tsv_file_path, key_column):
    """

    """
    import csv
    csv_rows = csv.reader(open(l_tsv_file_path, 'rb'), delimiter='\t')
    csv_dict = {}
    for row in csv_rows:
        csv_dict[row[key_column]] = row
    return csv_dict


def write_tsv_file(l_tsv_dict, l_tsv_file_path):
    """

    """
    import csv
    rows = l_tsv_dict.values()
    rows.sort(key=lambda x: x[0], reverse=False)
    with open(l_tsv_file_path, 'wb') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerows(rows)


def replace_line_in_file(file_path, line_str, replace_str):

    import fileinput

    for line in fileinput.input(file_path, inplace=1):
        if line.startswith(line_str):
            print(replace_str)
        else:
            print(line.strip('\n'))


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
