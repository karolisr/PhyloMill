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


def parse_line(line, delimiter, quotechar, stripchar):

    '''
    Will parse a line separated by delimiter and will return a list. Parts
    surrounded by quotechar will be grouped even if delimiter is found within.
    '''

    if not isinstance(line, basestring):
        return ''

    if not isinstance(delimiter, basestring):
        return ''

    #if not isinstance(quotechar, basestring):
    #    return ''

    #print(repr(line))

    nl = '\n'
    l_rem_nl = line.split(nl)[0]

    if stripchar:
        l_rem_nl = l_rem_nl.replace(stripchar, '')

    #print(repr(l_rem_nl))
    l_spl_q = [l_rem_nl]
    if quotechar:
        l_spl_q = l_rem_nl.split(quotechar)

        #print(l_spl_q)

        if l_rem_nl.startswith(quotechar):
            l_spl_q.remove('')
        if l_rem_nl.endswith(quotechar):
            l_spl_q.pop()

    #print(l_spl_q)

    final_list = list()

    for x in l_spl_q:
        if x.startswith(delimiter) or len(l_spl_q) == 1:
            if x is delimiter:
                continue
            else:
                x_spl_d = x.split(delimiter)
                if x.startswith(delimiter):
                    x_spl_d.remove('')
                if x.endswith(delimiter):
                    x_spl_d.pop()
                for y in x_spl_d:
                    y = y.strip()
                    final_list.append(y)
        else:
            final_list.append(x)

    #print(final_list)
    return(final_list)


def timestamp():
    import datetime
    return(str(datetime.datetime.now()).split('.')[0])


def compare_strings(l_strings):
    prev_s = None
    match = True
    for s in l_strings:
        if prev_s == None:
            prev_s = s
            continue
        if prev_s != s:
            match = False
            break
    return match


def in_range(x,a,b,percent):
    c = min(a,b)
    d = max(a,b)
    upper_limit = ((d-c)/(float(100)/float(percent)))+c
    if c<=x<=upper_limit:
        return True
    else:
        return False


def overlap(l_seg_a,l_seg_b):
    a1 = l_seg_a[0]
    a2 = l_seg_a[1]
    b1 = l_seg_b[0]
    b2 = l_seg_b[1]
    return_value = False
    if in_range(a1,b1,b2,100) or\
       in_range(a2,b1,b2,100) or\
       in_range(b1,a1,a2,100) or\
       in_range(b2,a1,a2,100):
        return_value = True
    return return_value
