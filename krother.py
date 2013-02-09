from __future__ import print_function
#from __future__ import unicode_literals


def attr(object):
    for property, value in vars(object).iteritems():
        print(property, " : ", value)
        print('---')
    return


def random_id(length):
    import os
    return str(abs(hash(str(os.urandom(200)))))[0:length]
