from __future__ import print_function
from __future__ import unicode_literals

def prepare_directory(path):

    '''
    Checks if directory at path exists and, if not, creates full path.
    '''

    import os
    if not os.path.exists(path):
        os.makedirs(path)
    return

# The commented file reading class from:
# http://www.mfasold.net/blog/2010/02/python-recipe-read-csvtsv-textfiles-and-ignore-comment-lines
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

def read_table_file(path, has_headers=False, headers=None, delimiter=b',',
    iterator=True):

    '''
    Reads a delimited text file.
    Returns:
        A list of dictionaries, one dictionary per row. Header names as keys.
    '''

    import csv
    handle = CommentedFile(open(path, 'rb'))
    if headers != None:
        has_headers = False
    if has_headers:
        headers = handle.next()
        headers = headers.replace('"', '')
        headers = headers.split('\n')[0]
        headers = headers.split(delimiter)
    rows_iterator = csv.DictReader(handle, delimiter=delimiter, quotechar=b'"',
        fieldnames=headers)
    return_value = rows_iterator
    if not iterator:
        rows_list = list()
        for row in rows_iterator:
            rows_list.append(row)
        return_value = rows_list
    return return_value

if __name__ == '__main__':
    
    # Tests
    
    import os

    PS = os.path.sep

    # CommentedFile
    handle = CommentedFile(open('testdata'+PS+'commentedfile.csv', 'rb'))
    for i, line in enumerate(handle):
        print(i, repr(line))
    handle.close()

    # read_table_file
    table = read_table_file(path='testdata'+PS+'commentedfile.csv',
        has_headers=True, headers=None, delimiter=b',', iterator=True)
    for i, line in enumerate(table):
        print(i, repr(line))
