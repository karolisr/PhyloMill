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

if __name__ == '__main__':
    
    # Tests
    
    import os

    PS = os.path.sep

    pass
