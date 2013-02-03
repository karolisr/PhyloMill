from __future__ import print_function
from __future__ import unicode_literals

def parse_directory(path, file_name_sep):

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
    
    file_list = os.listdir(path)
    
    return_list = list()
    
    for f in file_list:
    
        file_name = os.path.splitext(f)[0]
        file_ext = os.path.splitext(f)[1].split('.')[1]
        file_name_split = file_name.split(file_name_sep)
    
        file_dict = dict()
    
        file_dict['name'] = file_name
        file_dict['ext'] = file_ext
        file_dict['full'] = f
        file_dict['path'] = path + ps + f
        file_dict['split'] = file_name_split
    
        return_list.append(file_dict)
    
    return return_list

if __name__ == '__main__':
    
    # Tests
    
    import os

    PS = os.path.sep
    
    # parse_directory
    t_parse_directory = parse_directory('testdata'+PS+'parsedirectory', '$')
    for d in t_parse_directory:
        print(d)
