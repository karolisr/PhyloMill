from __future__ import print_function
#from __future__ import unicode_literals

def clear_line():
    print(chr(27) + "[2K", end='\r')


def hide_cursor():
    import os
    os.system('setterm -cursor off')


def show_cursor():
    import os
    os.system('setterm -cursor on')


def print_progress(current, total, length, prefix, postfix, show_bar=True):
    import sys
    #print(chr(27) + "[2K", end='\r')
    completed = int((float(current) / total) * length)

    bar = ''
    if show_bar:
        left = length - completed
        bar = '|' + '='*completed + '.'*left + '| '

    print(
        chr(27) + "[2K",
        prefix,
        bar,
        current,
        '/',
        total,
        ' ',
        '%.2f' % round((float(current) / total) * 100, 2),
        '%',
        postfix,
        sep='',
        end='\r'
        )

    sys.stdout.flush()

if __name__ == '__main__':

    # Tests

    import os
    import time

    PS = os.path.sep

    # print_progress
    hide_cursor()
    for i in range(1, 1001):
        time.sleep(0.001)
        print_progress(i, 1000, 50, '')
    show_cursor()
