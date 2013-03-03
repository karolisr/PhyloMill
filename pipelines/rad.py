#!/usr/bin/env python

from __future__ import print_function
#from __future__ import unicode_literals

if __name__ == '__main__':

    import sys
    import os
    import argparse

    import krio
    import krrad

    parser = argparse.ArgumentParser()

    parser.add_argument('--commands', type=unicode,
                        help='Commands.')

    parser.add_argument('--output_dir', type=unicode,
                        help='Output directory path.')

    parser.add_argument('--forward_file', type=unicode,
                        help='FASTQ file with forward RAD reads.')

    parser.add_argument('--reverse_file', type=unicode,
                        help='FASTQ file with reverse_file RAD reads.')

    parser.add_argument('--split_pieces', type=int,
                        help='Number of pieces to split FASTQ file into.')

    args = parser.parse_args()

    output_dir = None

    ps = os.path.sep

    if args.output_dir:
        output_dir = args.output_dir.rstrip(ps) + ps

    if not output_dir:
        print('Output directory is required.')
        sys.exit(1)
    if not args.commands:
        print('Commands are required.')
        sys.exit(1)
    else:

        commands = None
        if args.commands:
            commands = set([x.strip() for x in args.commands.split(',')])

        if commands and ('split' in commands):

            if not args.forward_file:
                print('File with forward RAD reads is required.')
                sys.exit(1)
            elif not args.split_pieces:
                print('split_pieces is required.')
                sys.exit(1)

            # Create output directory
            krio.prepare_directory(output_dir)

            krrad.split_rad_fastq_file(
                pieces=args.split_pieces,
                output_dir=output_dir,
                forward_reads_file_path=args.forward_file,
                reverse_reads_file_path=args.reverse_file
            )

#./rad.py
#--output_dir /home/karolis/Dropbox/code/test
#--commands split
#--forward_file /home/karolis/Dropbox/code/krpy/testdata/rad_forward.fastq
#--reverse_file /home/karolis/Dropbox/code/krpy/testdata/rad_reverse.fastq
#--split_pieces 6
