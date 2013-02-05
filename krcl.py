from __future__ import print_function
#from __future__ import unicode_literals

def hide_cursor():
    import os
    os.system('setterm -cursor off')

def show_cursor():
    import os
    os.system('setterm -cursor on')

def print_progress(current, total, length, prefix):
    #print(chr(27) + "[2K", end='\r')
    completed = int((float(current)/total) * length)
    left = length - completed
    print(
        prefix,
        '|',
        '=' * completed,
        '.' * left,
        '| ',
        current,
        '/',
        total,
        ' ',
        '%.2f' % round((float(current)/total) * 100, 2),
        '%',
        sep='',
        end='\r')

if __name__ == '__main__':
    
    # Tests
    
    import os
    import time

    PS = os.path.sep

    # print_progress
    hide_cursor()
    for i in range(1,1001):
        time.sleep(0.001)
        print_progress(i, 1000, 50, '')
    show_cursor()
        
