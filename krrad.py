from __future__ import print_function
#from __future__ import unicode_literals


def read_barcodes(file_path, delimiter, id_header, barcode_header):

    '''
    Read RAD barcodes and return a list of dictionaries with keys id and
    barcode.
    '''

    import krio

    barcodes = krio.read_table_file(path=file_path, has_headers=True,
        headers=None, delimiter=delimiter, quotechar='"', stripchar='',
        rettype='dict')

    id_header_new = 'id'
    barcode_header_new = 'barcode'

    for r in barcodes:
        if id_header_new != id_header:
            r[id_header_new] = r[id_header]
            del r[id_header]
        if barcode_header_new != barcode_header:
            r[barcode_header_new] = r[barcode_header]
            del r[barcode_header]

    return(barcodes)

if __name__ == '__main__':
    # Tests
    import os
    ps = os.path.sep

    # read_barcodes
    barcodes = read_barcodes(
        file_path='testdata' + ps + 'rad_barcodes.tsv',
        delimiter='\t',
        id_header='id',
        barcode_header='barcode')

    for r in barcodes:
        print(r)
