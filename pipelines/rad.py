#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
#from __future__ import unicode_literals

if __name__ == '__main__':

    import ConfigParser

    import sys
    import os
    import argparse
    import csv
    import subprocess
    from multiprocessing import Process
    from multiprocessing import JoinableQueue
    from multiprocessing import Manager

    from Bio.SeqIO.QualityIO import FastqGeneralIterator

    import krio
    import krbioio
    import krnextgen
    import krusearch

    ps = os.path.sep

    # Arguments ---------------------------------------------------------------

    parser = argparse.ArgumentParser()

    parser.add_argument('--commands', type=unicode, help='Commands.')

    parser.add_argument('--config', type=unicode,
                        help='Configuration file path.')

    # parser.add_argument('--output_dir', type=unicode,
    #                     help='Output directory path.')

#     parser.add_argument('--barcodes_file', type=unicode,
#                         help='Barcodes for RAD reads.')

#     parser.add_argument('--forward_file', type=unicode,
#                         help='FASTQ file with forward RAD reads.')

#     parser.add_argument('--reverse_file', type=unicode,
#                         help='FASTQ file with reverse_file RAD reads.')

#     parser.add_argument('--threads', type=int,
#                         help='Number of threads to use.')

#     # parser.add_argument('--output_file_format', type=unicode,
#     #                     help='Output file format.')

#     parser.add_argument('--max_barcode_mismatch_count', type=int,
#                         help='How many mismatches will be allowed before \
#                         barcode is considered bad.')

#     parser.add_argument('--trim_barcode', default=False, action='store_true',
#                         help='Should the barcode be trimmed.')

#     parser.add_argument('--trim_extra', type=int,
#                         help='How many extra bases will be trimmed after the \
#                         the barcode was trimmed. This can be used to trim a \
#                         restriction sites which will always be the same, etc.')

#     parser.add_argument('--quality_score_treshold', type=int,
#                         help='Minimum quality (phred) score required to \
#                         accept a site.')

#     parser.add_argument('--low_quality_residue', type=unicode,
#                         help='Symbol to use when masking.')

#     parser.add_argument('--max_prop_low_quality_sites', type=float,
#                         help='Maximum proportion of low quality sites in a \
#                         read for it to still be considered acceptable.')

#     parser.add_argument('--min_overlap', type=int,
#                         help='Minimum overlap required between forward and \
#                         reverse reads. If this overlap is not reached, the \
#                         reads will be concatenated.')

#     parser.add_argument('--mmmr_cutoff', type=float,
#                         help='When aligning forward and reverse reads to \
#                         check if they overlap, this value is used as a cutoff \
#                         when deciding if to accept or reject an alignment: \
#                         mmmr = match / (match + miss). miss does not include \
#                         ignored characters, by default \'N\'.')

#     parser.add_argument('--identity_threshold', type=float,
#                         help='Identity value for within-sample read \
#                         clustering.')

#     parser.add_argument('--min_seq_cluster', type=int,
#                         help='Minimum number of sequences in a cluster.')

#     parser.add_argument('--max_seq_cluster', type=int,
#                         help='Maximum number of sequences in a cluster.')

#     parser.add_argument('--error_rate_initial', type=float,
#                         help='Initial value of error rate (ε) for maximum \
#                         likelihood estimation.')

#     parser.add_argument('--heterozygosity_initial', type=float,
#                         help='Initial value of heterozygosity (π) for maximum \
#                         likelihood estimation.')

    args = parser.parse_args()

    # -------------------------------------------------------------------------

    # Check if prerequisites are met to run the pipeline ----------------------
    # if not output_dir:
    #     print('Output directory is required.')
    #     sys.exit(1)
    if not args.config:
        print('config_file is required.')
        sys.exit(1)
    if not args.commands:
        print('commands are required.')
        sys.exit(1)
    else:

        config = ConfigParser.SafeConfigParser()
        config.read(args.config)

        # Standardize output directory ----------------------------------------
        output_dir = config.get('General', 'output_directory')
        output_dir = output_dir.rstrip(ps) + ps

        # Determine which commands to run -------------------------------------
        commands = None
        if args.commands:
            commands = set([x.strip() for x in args.commands.split(',')])

        # Read_barcodes -------------------------------------------------------
        barcodes = config.get('General', 'barcodes_file')
        barcodes = krnextgen.read_barcodes(
            file_path=barcodes,
            delimiter='\t',
            id_header='id',
            barcode_header='barcode')

        # Number of CPU cores -------------------------------------------------
        cpu = config.getint('General', 'max_cpu_cores_available')

        split_raw_fastq_output_dir = output_dir + '01-raw-fastq-parts' + ps
        dmltplx_output_dir_split = (output_dir +
                                    '02-demultiplexed-fastq-parts' + ps)
        dmltplx_output_dir_combined = (output_dir +
                                       '03-demultiplexed-fastq-combined' + ps)
        masked_output_dir = output_dir + '04-masked-fastq' + ps
        binned_output_dir = output_dir + '05-binned-fasta' + ps
        clustered_output_dir = output_dir + '06-clustered' + ps
        sample_alignments_output_dir = output_dir + '07-sample-alignments' + ps
        analyzed_samples_output_dir = output_dir + '99-analyzed-samples' + ps
        krio.prepare_directory(analyzed_samples_output_dir)

        # Split FASTQ files ---------------------------------------------------
        if commands and ('split' in commands):

            # if not args.forward_file:
            #     print('File with forward RAD reads is required.')
            #     sys.exit(1)
            # if not args.threads:
            #     print('threads is required.')
            #     sys.exit(1)

            print()

            f_file = config.get('Demultiplex', 'forward_reads_file')
            r_file = config.get('Demultiplex', 'reverse_reads_file')

            krbioio.split_fastq_file(
                pieces=cpu,
                output_dir=split_raw_fastq_output_dir,
                forward_reads_file_path=f_file,
                reverse_reads_file_path=r_file
            )

        # Demultiplex split files ---------------------------------------------
        if commands and ('demultiplex' in commands):

            # if not barcodes:
            #     print('Barcodes are required.')
            #     sys.exit(1)
            # if not args.max_barcode_mismatch_count:
            #     print('max_barcode_mismatch_count is required.')
            #     sys.exit(1)
            # if not args.trim_barcode:
            #     print('\nBarcodes will not be trimmed!')
            # if not args.trim_extra:
            #     print('trim_extra is required.')
            #     sys.exit(1)

            file_list = krio.parse_directory(split_raw_fastq_output_dir, '_')
            file_list.sort(key=lambda x: x['name'], reverse=False)
            reverse = False
            for f in file_list:
                if f['split'][0] == 'r':
                    reverse = True
                    break

            print()

            processes = list()
            queue = JoinableQueue()

            def t(q):
                while True:
                    f = q.get()
                    print('Demultiplexing file', f['split'][1])
                    input_file_format = f['ext']
                    output_dir_split = (dmltplx_output_dir_split +
                                        f['split'][1])
                    reverse_reads_file_path = None
                    if reverse:
                        reverse_reads_file_path = (
                            split_raw_fastq_output_dir +
                            'r_' +
                            f['split'][1] + '.' +
                            f['ext'])
                    krnextgen.demultiplex(
                        barcodes=barcodes,
                        forward_reads_file_path=f['path'],
                        reverse_reads_file_path=reverse_reads_file_path,
                        input_file_format=input_file_format,
                        max_barcode_mismatch_count=config.getint(
                            'Demultiplex',
                            'max_bp_mismatch_in_barcode'),
                        output_dir=output_dir_split,
                        trim_barcode=True,
                        trim_extra=config.getint('Demultiplex', 'trim'),
                        write_every=1000
                    )
                    q.task_done()

            for f in file_list:
                if f['split'][0] == 'f':
                    queue.put(f)

            for i in range(cpu):
                worker = Process(target=t, args=(queue,))
                worker.start()
                processes.append(worker)

            queue.join()

            for p in processes:
                p.terminate()

            # Combine demultiplexed files
            print('\nCombining demultiplexed results...')
            krnextgen.combine_demultiplexed_results(
                input_dir=dmltplx_output_dir_split,
                output_dir=dmltplx_output_dir_combined)

            # Produce read lengths files
            combined_file_list = krio.parse_directory(
                dmltplx_output_dir_combined, '_')
            print('\nProducing read lengths files...\n')
            for f in combined_file_list:
                handle = open(f['path'], "rU")
                records = FastqGeneralIterator(handle)
                # records = SeqIO.parse(f['path'], 'fastq')
                lengths = list()
                for r in records:
                    lengths.append(len(r[1]))
                handle.close()
                print(f['split'][0], f['split'][-1], len(lengths), 'reads.')
                lengths_fp = (analyzed_samples_output_dir + f['split'][0] +
                              '_' + f['split'][-1] + '.lengths')
                handle = open(lengths_fp, 'wb')
                for l in lengths:
                    handle.write(str(l) + '\n')
                handle.close()

        # Mask low quality sites ----------------------------------------------
        if commands and ('mask' in commands):
            # if not args.quality_score_treshold:
            #     print('quality_score_treshold is required.')
            #     sys.exit(1)
            # if not args.low_quality_residue:
            #     print('low_quality_residue is required.')
            #     sys.exit(1)
            # if not args.threads:
            #     print('threads is required.')
            #     sys.exit(1)
            # if not args.output_file_format:
            #     print('output_file_format is required.')
            #     sys.exit(1)

            krio.prepare_directory(masked_output_dir)
            file_list = krio.parse_directory(dmltplx_output_dir_combined, ' ')

            print('\nMasking low quality sites...')

            processes = list()
            queue = JoinableQueue()

            def t(q):
                while True:
                    f = q.get()
                    handle_r = open(f['path'], "rU")

                    # records = FastqPhredIterator(handle_r)
                    records = FastqGeneralIterator(handle_r)
                    # records = SeqIO.parse(f['path'], 'fastq')
                    output_file_path = masked_output_dir + f['full']
                    handle_w = open(output_file_path, 'w')
                    for r in records:
                        masked = krnextgen.mask_low_quality_sites(
                            seq_str=r[1],
                            qual_str=r[2],
                            quality_score_treshold=(
                                config.getint('Mask', 'quality_score_treshold')
                            ),
                            low_quality_residue=(
                                config.get('General', 'low_quality_residue')))

                        # SeqIO.write(
                        #     sequences=r_masked,
                        #     handle=handle,
                        #     format='fastq'
                        # )

                        handle_w.write('@' + r[0] + '\n')
                        handle_w.write(masked + '\n')
                        handle_w.write('+\n')
                        handle_w.write(r[2] + '\n')

                    handle_r.close()
                    handle_w.close()
                    q.task_done()

            for f in file_list:
                if 'mismatch' not in f['name']:
                    queue.put(f)

            for i in range(cpu):
                worker = Process(target=t, args=(queue,))
                worker.start()
                processes.append(worker)

            queue.join()

            for p in processes:
                p.terminate()

        # Bin results by quality and produce forward and reverse consensus
        # sequences
        if commands and ('bin' in commands):
            # if not args.max_prop_low_quality_sites:
            #     print('max_prop_low_quality_sites is required.')
            #     sys.exit(1)
            # if not args.min_overlap:
            #     print('min_overlap is required.')
            #     sys.exit(1)
            # if not args.mmmr_cutoff:
            #     print('mmmr_cutoff is required.')
            #     sys.exit(1)
            # if not args.low_quality_residue:
            #     print('low_quality_residue is required.')
            #     sys.exit(1)
            # if not args.threads:
            #     print('threads is required.')
            #     sys.exit(1)

            krio.prepare_directory(binned_output_dir)
            file_list = krio.parse_directory(masked_output_dir, '_')

            print(('\nBinning results by quality and producing\n'
                   'forward and reverse consensus sequences...\n'))

            processes = list()
            queue = JoinableQueue()

            def t(q):
                while True:
                    f = q.get()

                    base_file_name = f['split'][0] + '_'
                    # + f['split'][1] + '_'

                    base_file_path = (masked_output_dir + base_file_name +
                                      f['split'][1] + '_')

                    print('Sample', f['split'][0])

                    # f_records = SeqIO.parse(base_file_path + 'f.fastq', 'fastq')
                    f_reads_handle = open(base_file_path + 'f.fastq', "rU")
                    f_reads = FastqGeneralIterator(f_reads_handle)

                    r_reads = None
                    if os.path.exists(base_file_path + 'r.fastq'):
                        # r_records = SeqIO.parse(base_file_path + 'r.fastq', 'fastq')
                        r_reads_handle = open(base_file_path + 'r.fastq', "rU")
                        r_reads = FastqGeneralIterator(r_reads_handle)

                    ofp_f_hq = (binned_output_dir + base_file_name +
                                'f_hq.fasta')
                    ofp_f_lq = (binned_output_dir + base_file_name +
                                'f_lq.fasta')

                    ofp_r_hq = None
                    ofp_r_lq = None

                    if r_reads:
                        ofp_r_hq = (binned_output_dir + base_file_name +
                                    'r_hq.fasta')
                        ofp_r_lq = (binned_output_dir + base_file_name +
                                    'r_lq.fasta')

                    ofp_all_hq = (binned_output_dir + base_file_name +
                                  'all_hq.fasta')

                    handle_f_hq = open(ofp_f_hq, 'w')
                    handle_f_lq = open(ofp_f_lq, 'w')

                    handle_r_hq = None
                    handle_r_lq = None

                    if r_reads:
                        handle_r_hq = open(ofp_r_hq, 'w')
                        handle_r_lq = open(ofp_r_lq, 'w')

                    handle_all_hq = open(ofp_all_hq, 'w')

                    lengths_path = (analyzed_samples_output_dir +
                                    f['split'][0] +  # '_' + f['split'][1] +
                                    '.binned')
                    handle_lengths = open(lengths_path, 'w')

                    reads_string = (
                        'fhq' + '\t' +
                        'rhq' + '\t' +
                        'cns' + '\t' +
                        'msg' + '\t' +
                        'flq' + '\t' +
                        'rlq' + '\n')

                    handle_lengths.write(reads_string)

                    for f_title, f_seq, f_qual in f_reads:
                        r_title = None
                        r_seq = None
                        r_qual = None
                        if r_reads:
                            r_title, r_seq, r_qual = r_reads.next()

                        binned = krnextgen.bin_reads(
                            title=f_title,
                            f_seq_str=f_seq,
                            r_seq_str=r_seq,
                            max_prop_low_quality_sites=
                            config.getfloat('Bin',
                                            'max_prop_low_quality_sites'),
                            min_overlap=config.getint('Bin', 'min_overlap'),
                            mmmr_cutoff=config.getfloat('Bin', 'mmmr_cutoff'),
                            low_quality_residue=
                            config.get('General', 'low_quality_residue')
                        )

                        f_hq = binned[0]
                        r_hq = binned[1]
                        consensus = binned[2]
                        consensus_title = binned[3]
                        consensus_message = binned[4]

                        ### COMMON FIND IN F READS
                        ### TGCAA GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG
                        ### FIND IN REVERSE READS
                        ### ACACTCTTTCCCTACACGACGCTCTTCCGATCT barcode TGCAG

                        if (consensus and len(consensus) <
                                config.getint('Bin', 'min_read_length')):
                            f_hq = False
                            r_hq = False
                            consensus = False
                            consensus_message = ''

                        str_fhq = ''
                        str_rhq = ''
                        str_cons = ''
                        str_msg = str(consensus_message)
                        str_flq = ''
                        str_rlq = ''

                        if f_hq:
                            # SeqIO.write(sequences=fr, handle=handle_f_hq,
                            #             format='fasta')
                            handle_f_hq.write('>'+f_title+'\n'+f_seq+'\n')
                            str_fhq = str(len(f_seq))
                        else:
                            # SeqIO.write(sequences=fr, handle=handle_f_lq,
                            #             format='fasta')
                            handle_f_lq.write('>'+f_title+'\n'+f_seq+'\n')
                            str_flq = str(len(f_seq))

                        if r_hq:
                            # SeqIO.write(sequences=rr, handle=handle_r_hq,
                            #             format='fasta')
                            handle_r_hq.write('>'+r_title+'\n'+r_seq+'\n')
                            str_rhq = str(len(r_seq))
                        elif r_reads:
                            # SeqIO.write(sequences=rr, handle=handle_r_lq,
                            #             format='fasta')
                            handle_r_lq.write('>'+r_title+'\n'+r_seq+'\n')
                            str_rlq = str(len(r_seq))

                        if f_hq and r_hq and consensus:
                            # SeqIO.write(sequences=consensus,
                            #             handle=handle_all_hq,
                            #             format='fasta')
                            handle_all_hq.write('>'+consensus_title+'\n' +
                                                consensus+'\n')
                            str_cons = str(len(consensus))

                        if f_hq and not consensus:
                            # SeqIO.write(sequences=fr, handle=handle_all_hq,
                            #             format='fasta')
                            handle_all_hq.write('>'+f_title+'\n'+f_seq+'\n')

                        if r_hq and not consensus:
                            # SeqIO.write(sequences=rr, handle=handle_all_hq,
                            #             format='fasta')
                            handle_all_hq.write('>'+r_title+'\n'+r_seq+'\n')

                        reads_string = (
                            str_fhq + '\t' +
                            str_rhq + '\t' +
                            str_cons + '\t' +
                            str_msg + '\t' +
                            str_flq + '\t' +
                            str_rlq + '\n')

                        handle_lengths.write(reads_string)

                    handle_lengths.close()

                    handle_f_hq.close()
                    handle_f_lq.close()

                    if r_reads:
                        handle_r_hq.close()
                        handle_r_lq.close()

                    handle_all_hq.close()
                    q.task_done()

            for f in file_list:
                if f['name'] != 'mismatch_f' and f['split'][-1] == 'f':
                    queue.put(f)

            for i in range(cpu):
                worker = Process(target=t, args=(queue,))
                worker.start()
                processes.append(worker)

            queue.join()

            for p in processes:
                p.terminate()

        # Cluster
        if commands and ('cluster' in commands):
            # if not args.threads:
            #     print('threads is required.')
            #     sys.exit(1)
            # if not args.identity_threshold:
            #     print('identity_threshold is required.')
            #     sys.exit(1)

            krio.prepare_directory(clustered_output_dir)
            file_list = krio.parse_directory(binned_output_dir, '_')

            print('\nSorting before clustering...\n')

            # Sort files, this will take memory, so we will not parallelize
            # this
            for f in file_list:
                if f['split'][-1] == 'hq' and f['split'][-2] == 'all':
                    print(f['full'])
                    ifp_split = os.path.splitext(f['path'])
                    subprocess.call(
                        ('usearch6' + ' -quiet' +
                            ' -sortbylength ' + f['path'] +
                            ' -output ' + ifp_split[0] + '_sorted' +
                            ifp_split[1]), shell=True)

            # We need to do this again as now we have new (sorted) files in the
            # directory
            file_list = krio.parse_directory(binned_output_dir, '_')

            print('\nClustering reads within samples...\n')

            processes = list()
            queue = JoinableQueue()

            def t(q):
                while True:
                    f = q.get()
                    # ifp_split = os.path.splitext(f['path'])
                    # input_file_path = ifp_split[0] + '_sorted' + ifp_split[1]
                    output_file = clustered_output_dir + f['name'] + '.uc'
                    print(f['full'])
                    krusearch.cluster_file(
                        input_file_path=f['path'],
                        output_file_path=output_file,
                        identity_threshold=
                        config.getfloat('Cluster', 'identity_threshold'),
                        sorted_input=True,
                        algorithm='smallmem',
                        strand='both',
                        threads=1,
                        quiet=True,
                        program=config.get('Cluster', 'usearch6_executable')
                    )
                    q.task_done()

            for f in file_list:
                if (f['split'][-1] == 'sorted' and
                    f['split'][-2] == 'hq' and
                        f['split'][-3] == 'all'):
                    queue.put(f)

            for i in range(cpu):
                worker = Process(target=t, args=(queue,))
                worker.start()
                processes.append(worker)

            queue.join()

            for p in processes:
                p.terminate()

        # Produce sample alignments and nucleotide counts per site
        if commands and ('align_samples' in commands):
            # if not args.min_seq_cluster:
            #     print('min_seq_cluster is required.')
            #     sys.exit(1)
            # if not args.max_seq_cluster:
            #     print('max_seq_cluster is required.')
            #     sys.exit(1)
            # if not args.threads:
            #     print('threads is required.')
            #     sys.exit(1)

            krio.prepare_directory(sample_alignments_output_dir)
            file_list = krio.parse_directory(clustered_output_dir, '_')

            print('\nProducing sample alignments and\n'
                  'nucleotide counts per site...\n')

            def t(f):
                fasta_file_path = binned_output_dir + f['name'] + '.fasta'
                aln_output_file_path = (
                    sample_alignments_output_dir +
                    f['name'] + '.alignment')
                counts_output_file_path = (
                    sample_alignments_output_dir +
                    f['name'] + '.counts')
                print(f['full'])
                cluster_depths = krnextgen.align_clusters(
                    min_seq_cluster=config.getint('General',
                                                  'min_seq_cluster'),
                    max_seq_cluster=config.getint('General',
                                                  'max_seq_cluster'),
                    uc_file_path=f['path'],
                    fasta_file_path=fasta_file_path,
                    aln_output_file_path=aln_output_file_path,
                    counts_output_file_path=counts_output_file_path,
                    threads=cpu)

                handle = open((analyzed_samples_output_dir +
                               f['split'][0] +
                               #  '_' + f['split'][1] +
                               '.clusters'), 'wb')
                for c in cluster_depths:
                    handle.write(str(c) + '\n')
                handle.close()

                ns = krnextgen.nt_site_counts(counts_output_file_path, 1, 0)
                coverage_list = [sum(x) for x in ns]
                handle = open((analyzed_samples_output_dir +
                               f['split'][0] +
                               #  '_' + f['split'][1] +
                               '.coverage'), 'wb')
                for c in coverage_list:
                    handle.write(str(c) + '\n')
                handle.close()

            for f in file_list:
                if (f['split'][-1] == 'sorted' and
                    f['split'][-2] == 'hq' and
                        f['split'][-3] == 'all'):
                    t(f)

#             # processes = list()
#             # queue = JoinableQueue()

#             # def t(q):
#             #     while True:
#             #         q_input = q.get()
#             #         f = q_input[0]
#             #         # q_id = str(q_input[1])
#             #         fasta_file_path = binned_output_dir + f['name'] + '.fasta'
#             #         aln_output_file_path = (
#             #             sample_alignments_output_dir +
#             #             f['name'] + '.alignment')
#             #         counts_output_file_path = (
#             #             sample_alignments_output_dir +
#             #             f['name'] + '.counts')
#             #         print('File:', f['full'])
#             #         cluster_depths = krnextgen.align_clusters(
#             #             min_seq_cluster=args.min_seq_cluster,
#             #             max_seq_cluster=args.max_seq_cluster,
#             #             uc_file_path=f['path'],
#             #             fasta_file_path=fasta_file_path,
#             #             aln_output_file_path=aln_output_file_path,
#             #             counts_output_file_path=counts_output_file_path,
#             #             # temp_dir_path=sample_alignments_output_dir,
#             #             # temp_file_id=q_id,
#             #             threads=args.threads)

#             #         handle = open((analyzed_samples_output_dir +
#             #                        f['split'][0] + '_' +
#             #                        f['split'][1] +
#             #                        '_clusters.csv'), 'wb')
#             #         for c in cluster_depths:
#             #             handle.write(str(c) + '\n')
#             #         handle.close()

#             #         q.task_done()

#             # for i, f in enumerate(file_list):
#             #     if (f['split'][-1] == 'sorted' and
#             #         f['split'][-2] == 'hq' and
#             #             f['split'][-3] == 'all'):
#             #         queue.put([f, i])

#             ### THREADS HARD-CODED ###
#             # for i in range(1):
#             #     worker = Process(target=t, args=(queue,))
#             #     worker.start()
#             #     processes.append(worker)

#             # queue.join()

#             # for p in processes:
#             #     p.terminate()

        # Estimate error rate, heterozygosity. Produce some statistics
        if commands and ('analyze_samples' in commands):
            # if not args.min_seq_cluster:
            #     print('min_seq_cluster is required.')
            #     sys.exit(1)
            # if not args.max_seq_cluster:
            #     print('max_seq_cluster is required.')
            #     sys.exit(1)
            # if not args.error_rate_initial:
            #     print('error_rate_initial is required.')
            #     sys.exit(1)
            # if not args.heterozygosity_initial:
            #     print('heterozygosity_initial is required.')
            #     sys.exit(1)
            # if not args.threads:
            #     print('threads is required.')
            #     sys.exit(1)

            file_list = krio.parse_directory(sample_alignments_output_dir, '_')

            print('\nAnalyzing samples...\n')

            manager = Manager()
            results = manager.list()

            processes = list()
            queue = JoinableQueue()

            def t(q):
                while True:
                    q_input = q.get()
                    f = q_input[0]
                    results = q_input[1]
                    print(f['split'][0], 'starting...')
                    p = krnextgen.nt_freq(f['path'])
                    # print(f['split'][0], p)
                    ns = krnextgen.nt_site_counts(
                        f['path'],
                        config.getint('General', 'min_seq_cluster'),
                        config.getint('General', 'max_seq_cluster'))
                    mle = krnextgen.mle_e_and_pi(
                        ns=ns,
                        p=p,
                        e0=config.getfloat('e and pi', 'error_rate_initial'),
                        pi0=config.getfloat('e and pi',
                                            'heterozygosity_initial'))

                    results_dict = dict()

                    results_dict['sample'] = str(f['split'][0])
                    results_dict['barcode'] = str(f['split'][1])
                    results_dict['fA'] = p[0]
                    results_dict['fC'] = p[1]
                    results_dict['fG'] = p[2]
                    results_dict['fT'] = p[3]
                    results_dict['e'] = mle[0]
                    results_dict['pi'] = mle[1]
                    results_dict['negll'] = mle[2]

#                     # rps = float(sum(coverage_list)) / float(len(coverage_list))
#                     # results_dict['rps'] = rps
#                     # results_dict['sites'] = len(coverage_list)

#                     # cluster_list = list()
#                     # handle = open((analyzed_samples_output_dir +
#                     #                results_dict['sample'] + '_' +
#                     #                results_dict['barcode'] +
#                     #                '_clusters.csv'), 'rb')
#                     # for c in handle:
#                     #     cluster_list.append(int(c))
#                     # handle.close()
#                     # rpc = float(sum(cluster_list)) / float(len(cluster_list))
#                     # results_dict['rpc'] = rpc
#                     # results_dict['clusters'] = len(cluster_list)

                    handle = open((analyzed_samples_output_dir +
                                   results_dict['sample'] +
                                   # '_' +
                                   # results_dict['barcode'] +
                                   '.stats'), 'wb')

                    handle.write('sample\t' + str(results_dict['sample']) +
                                 '\n')
                    handle.write('barcode\t' + str(results_dict['barcode']) +
                                 '\n')
#                     # reads per cluster (rpc)
#                     # handle.write('rpc\t' + str(rpc) + '\n')
#                     # handle.write('clusters\t' +
#                     #              str(results_dict['clusters']) + '\n')
#                     # reads per site (rps)
#                     # handle.write('rps\t' + str(rps) + '\n')
#                     # handle.write('sites\t' + str(results_dict['sites']) + '\n')
                    handle.write('fA\t' + str(results_dict['fA']) + '\n')
                    handle.write('fC\t' + str(results_dict['fC']) + '\n')
                    handle.write('fG\t' + str(results_dict['fG']) + '\n')
                    handle.write('fT\t' + str(results_dict['fT']) + '\n')
                    handle.write('e\t' + str(results_dict['e']) + '\n')
                    handle.write('pi\t' + str(results_dict['pi']) + '\n')
                    handle.write('negll\t' + str(results_dict['negll']) + '\n')

                    handle.close()

                    results.append(results_dict)

                    print(f['split'][0], 'done.')

                    q.task_done()

            for f in file_list:
                if f['ext'] == 'counts':
                    queue.put([f, results])

            for i in range(cpu):
                worker = Process(target=t, args=(queue,))
                worker.start()
                processes.append(worker)

            queue.join()

            for p in processes:
                p.terminate()

            results_new = list()
            for r in results:
                results_new.append(r)

            results_new.sort(key=lambda x: x['sample'], reverse=False)

            with open(analyzed_samples_output_dir+'stats.csv', 'wb') as f:
                # fieldnames = ['sample', 'barcode', 'rpc', 'clusters', 'rps',
                #               'sites', 'fA', 'fC', 'fG', 'fT', 'e', 'pi',
                #               'negll']
                fieldnames = ['sample', 'barcode', 'fA', 'fC', 'fG', 'fT', 'e',
                              'pi', 'negll']
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writerow({
                    'sample': 'sample',
                    'barcode': 'barcode',
                    # 'rpc': 'rpc',
                    # 'clusters': 'clusters',
                    # 'rps': 'rps',
                    # 'sites': 'sites',
                    'fA': 'fA',
                    'fC': 'fC',
                    'fG': 'fG',
                    'fT': 'fT',
                    'e': 'e',
                    'pi': 'pi',
                    'negll': 'negll',
                })
                writer.writerows(results_new)


# # p = nt_freq('/home/karolis/Dropbox/code/krpy/testdata/nt.counts')
# # print(p)
# # ns = nt_site_counts('/home/karolis/Dropbox/code/krpy/testdata/nt.counts')
# # mle = mle_e_and_pi(ns, p, e0=0.001, pi0=0.001)
# # print(mle)

# # --commands split,demultiplex,mask,bin,cluster,align_samples,analyze_samples

# # time ./rad.py \
# # --output_dir /Users/karolis/Dropbox/code/test/rad \
# # --forward_file /Users/karolis/Dropbox/code/krpy/testdata/rad_forward.fastq \
# # --reverse_file /Users/karolis/Dropbox/code/krpy/testdata/rad_reverse.fastq \
# # --threads 4 \
# # --barcodes '/Users/karolis/Dropbox/code/krpy/testdata/rad_barcodes.tsv' \
# # --max_barcode_mismatch_count 1 \
# # --trim_barcode \
# # --trim_extra 5 \
# # --quality_score_treshold 30 \
# # --low_quality_residue N \
# # --max_prop_low_quality_sites 0.10 \
# # --min_overlap 5 \
# # --mmmr_cutoff 0.85 \
# # --identity_threshold 0.90 \
# # --min_seq_cluster 2 \
# # --max_seq_cluster 1000 \
# # --error_rate_initial 0.001 \
# # --heterozygosity_initial 0.001 \
# # --commands split

# # time ./rad.py \
# # --output_dir /data/gbs-new \
# # --forward_file /data/gbs-andy-david/green_dzaya_DZAYA_R1.PF.fastq \
# # --reverse_file /data/gbs-andy-david/green_dzaya_DZAYA_R2.PF.fastq \
# # --threads 6 \
# # --barcodes /data/gbs-andy-david/barcodes.tsv \
# # --max_barcode_mismatch_count 1 \
# # --trim_barcode \
# # --trim_extra 5 \
# # --quality_score_treshold 30 \
# # --low_quality_residue N \
# # --max_prop_low_quality_sites 0.10 \
# # --min_overlap 5 \
# # --mmmr_cutoff 0.85 \
# # --identity_threshold 0.90 \
# # --min_seq_cluster 10 \
# # --max_seq_cluster 1000 \
# # --error_rate_initial 0.0001 \
# # --heterozygosity_initial 0.001 \
# # --commands analyze_samples
