#!/usr/bin/env python

from __future__ import print_function
#from __future__ import unicode_literals

if __name__ == '__main__':

    import sys
    import os
    import argparse

    from multiprocessing import Process
    from Queue import Queue
    from threading import Thread

    from Bio import SeqIO

    import krio
    import krbioio
    import krpipe
    import krnextgen

    ps = os.path.sep

    # Possible arguments ------------------------------------------------------
    parser = argparse.ArgumentParser()

    parser.add_argument('--commands', type=unicode,
                        help='Commands.')
    parser.add_argument('--output_dir', type=unicode,
                        help='Output directory path.')
    parser.add_argument('--barcodes_file', type=unicode,
                        help='Barcodes for RAD reads.')
    parser.add_argument('--forward_file', type=unicode,
                        help='FASTQ file with forward RAD reads.')
    parser.add_argument('--reverse_file', type=unicode,
                        help='FASTQ file with reverse_file RAD reads.')
    parser.add_argument('--threads', type=int,
                        help='Number of threads to use.')
    parser.add_argument('--output_file_format', type=unicode,
                        help='Output file format.')
    parser.add_argument('--max_barcode_mismatch_count', type=int,
                        help='How many mismatches will be allowed before \
                        barcode is considered bad.')
    parser.add_argument('--trim_barcode', default=False, action='store_true',
                        help='Should the barcode be trimmed.')
    parser.add_argument('--trim_extra', type=int,
                        help='How many extra bases will be trimmed after the \
                        the barcode was trimmed. This can be used to trim a \
                        restriction sites which will always be the same, etc.')
    parser.add_argument('--quality_score_treshold', type=int,
                        help='Minimum quality (phred) score required to \
                        accept a site.')
    parser.add_argument('--low_quality_residue', type=unicode,
                        help='Symbol to use when masking.')

    args = parser.parse_args()

    # Standardize output directory --------------------------------------------
    output_dir = None
    if args.output_dir:
        output_dir = args.output_dir.rstrip(ps) + ps

    # Check if prerequisites are met to run the pipeline ----------------------
    if not output_dir:
        print('Output directory is required.')
        sys.exit(1)
    if not args.commands:
        print('Commands are required.')
        sys.exit(1)
    else:

        # Determine which commands to run -------------------------------------
        commands = None
        if args.commands:
            commands = set([x.strip() for x in args.commands.split(',')])

        # Read_barcodes -------------------------------------------------------
        barcodes = None
        if args.barcodes_file:
            barcodes = krnextgen.read_barcodes(
                file_path=args.barcodes_file,
                delimiter='\t',
                id_header='id',
                barcode_header='barcode')

        split_raw_fastq_output_dir = output_dir + '01-raw-fastq-parts'
        dmltplx_output_dir_split = output_dir + '02-demultiplexed-fastq-parts'
        dmltplx_output_dir_combined = (output_dir +
                                       '03-demultiplexed-fastq-combined')
        masked_output_dir = (output_dir + '04-masked-fastq')

        # Split FASTQ files ---------------------------------------------------
        if commands and ('split' in commands):

            if not args.forward_file:
                print('File with forward RAD reads is required.')
                sys.exit(1)
            elif not args.threads:
                print('threads is required.')
                sys.exit(1)

            print('\n')

            krbioio.split_fastq_file(
                pieces=args.threads,
                output_dir=split_raw_fastq_output_dir,
                forward_reads_file_path=args.forward_file,
                reverse_reads_file_path=args.reverse_file
            )

        # Demultiplex split files ---------------------------------------------
        if commands and ('demultiplex' in commands):

            if not barcodes:
                print('Barcodes are required.')
                sys.exit(1)
            if not args.max_barcode_mismatch_count:
                print('max_barcode_mismatch_count is required.')
                sys.exit(1)
            if not args.trim_barcode:
                print('\nBarcodes will not be trimmed!')
            if not args.trim_extra:
                print('trim_extra is required.')
                sys.exit(1)

            split_file_list = krpipe.parse_directory(
                split_raw_fastq_output_dir, '_')
            split_file_list.sort(key=lambda x: x['name'], reverse=False)
            reverse = False
            for f in split_file_list:
                if f['split'][0] == 'r':
                    reverse = True
                    break

            processes = list()

            print('')

            for f in split_file_list:
                if f['split'][0] == 'f':
                    print('Demultiplexing file', f['split'][1])
                    input_file_format = f['ext']
                    output_dir_split = (dmltplx_output_dir_split + ps +
                                        f['split'][1])
                    reverse_reads_file_path = None
                    if reverse:
                        reverse_reads_file_path = (
                            split_raw_fastq_output_dir + ps +
                            'r_' +
                            f['split'][1] + '.' +
                            f['ext'])

                    p = Process(
                        target=krnextgen.demultiplex,
                        args=(
                            barcodes,
                            # forward_reads_file_path
                            f['path'],
                            reverse_reads_file_path,
                            input_file_format,
                            # output_file_format
                            'fastq',
                            args.max_barcode_mismatch_count,
                            # output_dir
                            output_dir_split,
                            args.trim_barcode,
                            args.trim_extra,
                            # write_every
                            1000
                        )
                    )

                    p.start()

                    processes.append(p)

            for p in processes:
                p.join()

            print('\nCombining demultiplexed results...')

            krnextgen.combine_demultiplexed_results(
                input_dir=dmltplx_output_dir_split,
                output_dir=dmltplx_output_dir_combined)

        # Mask low quality sites ----------------------------------------------
        if commands and ('mask' in commands):
            if not args.quality_score_treshold:
                print('quality_score_treshold is required.')
                sys.exit(1)
            if not args.low_quality_residue:
                print('low_quality_residue is required.')
                sys.exit(1)
            if not args.output_file_format:
                print('output_file_format is required.')
                sys.exit(1)

            masked_output_dir = masked_output_dir.rstrip(ps) + ps
            krio.prepare_directory(masked_output_dir)
            file_list = krpipe.parse_directory(dmltplx_output_dir_combined,
                                               ' ')

            print('\nMasking low quality sites.')

            queue = Queue()

            def t(q):
                while True:
                    f = q.get()
                    records = SeqIO.parse(f['path'], 'fastq')
                    output_file_path = masked_output_dir + f['full']
                    handle = open(output_file_path, 'w')
                    for r in records:
                        r_masked = krnextgen.mask_low_quality_sites(
                            bio_seq_record=r,
                            quality_score_treshold=args.quality_score_treshold,
                            low_quality_residue=str(args.low_quality_residue))
                        SeqIO.write(r_masked, handle, args.output_file_format)
                    handle.close()
                    q.task_done()

            for i in range(args.threads):
                worker = Thread(target=t, args=(queue,))
                worker.setDaemon(True)
                worker.start()

            for f in file_list:
                queue.put(f)

            queue.join()

# ./rad.py \
# --output_dir /home/karolis/Dropbox/code/test/rad \
# --commands split,demultiplex,mask \
# --forward_file /home/karolis/Dropbox/code/krpy/testdata/rad_forward.fastq \
# --reverse_file /home/karolis/Dropbox/code/krpy/testdata/rad_reverse.fastq \
# --threads 6 \
# --barcodes '/home/karolis/Dropbox/code/krpy/testdata/rad_barcodes.tsv' \
# --max_barcode_mismatch_count 1 \
# --trim_barcode \
# --output_file_format fastq \
# --trim_extra 5 \
# --quality_score_treshold 30 \
# --low_quality_residue N

# ./rad.py \
# --output_dir /data/gbs-new \
# --forward_file /data/gbs-andy-david/green_dzaya_DZAYA_R1.PF.fastq \
# --reverse_file /data/gbs-andy-david/green_dzaya_DZAYA_R2.PF.fastq \
# --threads 6 \
# --barcodes /data/gbs-andy-david/barcodes.tsv \
# --max_barcode_mismatch_count 1 \
# --trim_barcode \
# --output_file_format fastq \
# --trim_extra 5 \
# --quality_score_treshold 30 \
# --low_quality_residue N \
# --commands split,demultiplex,mask
