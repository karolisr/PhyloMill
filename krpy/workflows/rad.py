#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
# from __future__ import unicode_literals


def write_log(msg, log_file_path, append=True):
    log_handle = None
    if append:
        log_handle = open(log_file_path, 'a')
    else:
        log_handle = open(log_file_path, 'wb')
    log_handle.write(msg)
    log_handle.write('\n')
    log_handle.close()
    return()


if __name__ == '__main__':

    import ConfigParser

    import sys
    import os
    import argparse
    import csv
    import subprocess
    import string
    import datetime
    import random
    import shutil

    from multiprocessing import Process
    from multiprocessing import JoinableQueue
    from multiprocessing import Manager

    import numpy

    import datrie

    from Bio.SeqIO.QualityIO import FastqGeneralIterator
    from Bio import AlignIO

    import krio
    import krbioio
    import krnextgen
    import krusearch
    import kriupac
    import kralign
    import krother

    # Alignment concatenation function is recursive.
    sys.setrecursionlimit(500000)

    ps = os.path.sep

    # Arguments ---------------------------------------------------------------

    parser = argparse.ArgumentParser()

    parser.add_argument('--commands', type=unicode, help='Commands.')

    parser.add_argument('--config', type=unicode,
                        help='Configuration file path.')

    parser.add_argument('--group', type=unicode,
                        help='Sample group to analyze.')

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
        print('--config is required.')
        sys.exit(1)
    if not args.commands:
        print('--commands are required.')
        sys.exit(1)
    else:

        config = ConfigParser.SafeConfigParser()
        config.read(args.config)

        numpy.random.seed(config.getint('General', 'random_seed'))

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
        consensus_output_dir = output_dir + '08-consensus' + ps
        grouped_consensus_output_dir = output_dir + '09-consensus-grouped' + ps
        between_sample_clusters_output_dir = output_dir + '10-between-sample-clusters' + ps
        between_sample_alignments_output_dir = output_dir + '11-between-sample-alignments' + ps
        analyzed_samples_output_dir = output_dir + '99-analyzed-samples' + ps
        krio.prepare_directory(analyzed_samples_output_dir)

        # Group samples
        sample_groups_dict = dict()
        # sample_groups_dict = datrie.Trie(string.printable)
        sample_groups = config.items('Sample Groups')
        for group in sample_groups:
            if group[0] != 'outgroup_samples':
                samples = group[1].replace('\n', '').split(',')
                sample_groups_dict[group[0]] = samples
        outgroups = config.get('Sample Groups', 'outgroup_samples').replace('\n', '').split(',')

        # print(outgroups)
        # for k in sample_groups_dict.keys():
        #     print(k, sample_groups_dict[k])

        # Map to Reference options
        groups_map_ref_loci_dict = dict()
        # groups_map_ref_loci_dict = datrie.Trie(string.printable)
        groups_map_to_ref_step_dict = dict()
        # groups_map_to_ref_step_dict = datrie.Trie(string.printable)
        map_to_reference_options = config.items('Map to Reference')
        for mtro in map_to_reference_options:
            mtro_name = mtro[0].split('.')
            if mtro_name[1] == 'infer_loci_using_reference':
                map_bool = config.getboolean('Map to Reference', mtro[0])
                groups_map_ref_loci_dict[mtro_name[0]] = map_bool
            else:
                map_to_ref_options = mtro_name[1].split('_')
                map_to_ref_options[0] = int(map_to_ref_options[0])
                map_to_ref_options.append(mtro[1])
                if mtro_name[0] not in groups_map_to_ref_step_dict.keys():
                    groups_map_to_ref_step_dict[mtro_name[0]] = list()
                groups_map_to_ref_step_dict[mtro_name[0]].append(map_to_ref_options)
        for key in groups_map_to_ref_step_dict.keys():
            groups_map_to_ref_step_dict[key].sort(key=lambda x: x[0], reverse=False)

        lfp = output_dir + 'log.txt'

        start_time = datetime.datetime.now()

        print()
        msg = krother.timestamp()
        print(msg)
        if os.path.exists(lfp):
            write_log('', lfp, append=True)
            write_log('--- Logging started - ' + msg + ' ------------------------------------', lfp, append=True)
        else:
            write_log('--- Logging started - ' + msg + ' ------------------------------------', lfp, append=False)

        print()
        write_log('', lfp)

        msg = 'Using ' + str(cpu) + ' CPU cores.\n'
        print(msg)
        write_log(msg, lfp)

        # Split FASTQ files ---------------------------------------------------
        if commands and ('split' in commands):

            # if not args.forward_file:
            #     print('File with forward RAD reads is required.')
            #     sys.exit(1)
            # if not args.threads:
            #     print('threads is required.')
            #     sys.exit(1)

            f_file = config.get('General', 'forward_reads_file')
            r_file = config.get('General', 'reverse_reads_file')

            krbioio.split_fastq_file(
                pieces=cpu,
                output_dir=split_raw_fastq_output_dir,
                forward_reads_file_path=f_file,
                reverse_reads_file_path=r_file,
                log_func=write_log,
                log_file_path=lfp
            )

            print()
            write_log('', lfp)

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

            processes = list()
            queue = JoinableQueue()

            trim_extra = len(config.get('General', 'f_sticky'))

            def t(q):
                while True:
                    f = q.get()
                    # lock.acquire()
                    msg = krother.timestamp() + ' - Demultiplexing File ' + f['split'][1]
                    print(msg)
                    write_log(msg, lfp)
                    # lock.release()
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
                        trim_extra=trim_extra,
                        write_every=1000
                    )
                    q.task_done()

            for f in file_list:
                if f['split'][0] == 'f':
                    queue.put(f)

            for i in range(cpu):
                # lock = Lock()
                worker = Process(target=t, args=(queue,))
                worker.start()
                processes.append(worker)

            queue.join()

            for p in processes:
                p.terminate()

            # Combine demultiplexed files
            msg = '\nCombining demultiplexed results...'
            print(msg)
            write_log(msg, lfp)
            krnextgen.combine_demultiplexed_results(
                input_dir=dmltplx_output_dir_split,
                output_dir=dmltplx_output_dir_combined)

            # Produce read lengths files
            combined_file_list = krio.parse_directory(
                dmltplx_output_dir_combined, '_')
            msg = '\nProducing read lengths files...\n'
            print(msg)
            write_log(msg, lfp)
            mismatch = False
            for f in combined_file_list:
                if f['isdir']:
                    continue
                handle = open(f['path'], "rU")
                records = FastqGeneralIterator(handle)
                # records = SeqIO.parse(f['path'], 'fastq')
                lengths = list()
                for r in records:
                    lengths.append(len(r[1]))
                handle.close()
                sample = 'Sample '
                if f['split'][0] == 'Mismatch':
                    if not mismatch:
                        mismatch = True
                        sample = '\n'
                    else:
                        sample = ''

                msg = sample + f['split'][0] + ' ' + f['split'][-1].upper() + ', ' + str(len(lengths)) + ' reads.'
                print(msg)
                write_log(msg, lfp)
                lengths_fp = (analyzed_samples_output_dir + f['split'][0] +
                              '_' + f['split'][-1] + '.lengths')
                handle = open(lengths_fp, 'wb')
                for l in lengths:
                    handle.write(str(l) + '\n')
                handle.close()

            # Move demultiplexed files to directories by sample groups
            msg = '\nMoving demultiplexed results to directories by group...'
            print(msg)
            write_log(msg, lfp)

            file_list = krio.parse_directory(
                path=dmltplx_output_dir_combined,
                file_name_sep='_',
                sort='forward'
            )

            for group in sample_groups_dict.keys():
                group_samples = sample_groups_dict[group]

                output_dir = dmltplx_output_dir_combined + group + ps
                krio.prepare_directory(output_dir)

                for sample_name in group_samples:
                    for f in file_list:
                        output_file_path = output_dir + f['full']
                        if f['name'].startswith(sample_name):
                            shutil.move(f['path'], output_file_path)

            print()
            write_log('', lfp)

        # Mask low quality sites ----------------------------------------------
        if commands and ('mask' in commands):
            if not args.group:
                print('Sample group is required.')
                sys.exit(1)

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

            masked_output_dir_sample = masked_output_dir + args.group + ps
            dmltplx_output_dir_combined_sample = dmltplx_output_dir_combined + args.group + ps

            krio.prepare_directory(masked_output_dir_sample)
            file_list = krio.parse_directory(dmltplx_output_dir_combined_sample, ' ')

            msg = 'Masking low quality sites...\n'
            print(msg)
            write_log(msg, lfp)

            processes = list()
            queue = JoinableQueue()

            def t(q):
                while True:
                    f = q.get()
                    handle_r = open(f['path'], "rU")

                    f_name = f['split'][0].split('_')

                    msg = krother.timestamp() + ' - Sample ' + f_name[0] + ' ' + f_name[2].upper() + ' starting...'
                    print(msg)
                    write_log(msg, lfp)

                    # records = FastqPhredIterator(handle_r)
                    records = FastqGeneralIterator(handle_r)
                    # records = SeqIO.parse(f['path'], 'fastq')
                    output_file_path = masked_output_dir_sample + f['full']
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

                    msg = krother.timestamp() + ' - Sample ' + f_name[0] + ' ' + f_name[2].upper() + ' done.'
                    print(msg)
                    write_log(msg, lfp)

                    q.task_done()

            for f in file_list:
                if 'Mismatch' not in f['name']:
                    queue.put(f)

            for i in range(cpu):
                worker = Process(target=t, args=(queue,))
                worker.start()
                processes.append(worker)

            queue.join()

            for p in processes:
                p.terminate()

            print()
            write_log('', lfp)

        # Bin results by quality and produce forward and reverse consensus
        # sequences
        if commands and ('bin' in commands):
            if not args.group:
                print('Sample group is required.')
                sys.exit(1)

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

            masked_output_dir_sample = masked_output_dir + args.group + ps
            binned_output_dir_sample = binned_output_dir + args.group + ps

            krio.prepare_directory(binned_output_dir_sample)
            file_list = krio.parse_directory(masked_output_dir_sample, '_')

            msg = ('Binning results by quality and producing\n'
                   'forward and reverse consensus sequences...\n')
            print(msg)
            write_log(msg, lfp)

            processes = list()
            queue = JoinableQueue()

            def t(q):
                while True:
                    f = q.get()

                    base_file_name = f['split'][0] + '_'
                    # + f['split'][1] + '_'

                    base_file_path = (masked_output_dir_sample + base_file_name +
                                      f['split'][1] + '_')

                    msg = krother.timestamp() + ' - Sample ' + f['split'][0] + ' starting...'
                    print(msg)
                    write_log(msg, lfp)

                    # f_records = SeqIO.parse(base_file_path + 'f.fastq', 'fastq')
                    f_reads_handle = open(base_file_path + 'f.fastq', "rU")
                    f_reads = FastqGeneralIterator(f_reads_handle)

                    r_reads = None
                    if os.path.exists(base_file_path + 'r.fastq'):
                        # r_records = SeqIO.parse(base_file_path + 'r.fastq', 'fastq')
                        r_reads_handle = open(base_file_path + 'r.fastq', "rU")
                        r_reads = FastqGeneralIterator(r_reads_handle)

                    ofp_f_hq = (binned_output_dir_sample + base_file_name +
                                'f_hq.fasta')
                    ofp_f_lq = (binned_output_dir_sample + base_file_name +
                                'f_lq.fasta')

                    ofp_r_hq = None
                    ofp_r_lq = None

                    if r_reads:
                        ofp_r_hq = (binned_output_dir_sample + base_file_name +
                                    'r_hq.fasta')
                        ofp_r_lq = (binned_output_dir_sample + base_file_name +
                                    'r_lq.fasta')

                    ofp_all_hq = (binned_output_dir_sample + base_file_name +
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

                    # Forward read oligo components
                    # These can be found in overreaching reverse reads
                    # barcode_adapter-barcode-f_sticky
                    barcode_adapter = config.get('General', 'barcode_adapter')
                    barcode = f['split'][1]
                    f_sticky = config.get('General', 'f_sticky')
                    f_oligo = barcode_adapter + barcode.upper() + f_sticky

                    # Reverse read oligo components
                    # These can be found in overreaching forward reads
                    # r_sticky-common_adapter
                    common_adapter = config.get('General', 'common_adapter')
                    r_sticky = config.get('General', 'r_sticky')
                    r_oligo = r_sticky + common_adapter

                    max_prop_low_quality_sites = config.getfloat(
                        'Bin', 'max_prop_low_quality_sites')
                    min_overlap = config.getint('Bin', 'min_overlap')
                    mmmr_cutoff = config.getfloat('Bin', 'mmmr_cutoff')
                    concatenate = config.getboolean('Bin', 'concatenate')
                    low_quality_residue = config.get(
                        'General', 'low_quality_residue')
                    min_read_length = config.getint('Bin', 'min_read_length')

                    for f_title, f_seq, f_qual in f_reads:
                        r_title = None
                        r_seq = None
                        r_qual = None
                        if r_reads:
                            r_title, r_seq, r_qual = r_reads.next()

                        # print(concatenate)

                        binned = krnextgen.bin_reads(
                            title=f_title,
                            f_seq_str=f_seq,
                            r_seq_str=r_seq,
                            max_prop_low_quality_sites=max_prop_low_quality_sites,
                            min_overlap=min_overlap,
                            mmmr_cutoff=mmmr_cutoff,
                            concatenate=concatenate,
                            low_quality_residue=low_quality_residue,
                            f_oligo=f_oligo,
                            r_oligo=r_oligo,
                            min_read_length=min_read_length
                        )

                        f_hq = binned[0]
                        r_hq = binned[1]
                        f_seq = binned[2]
                        r_seq = binned[3]
                        consensus = binned[4]
                        consensus_title = binned[5]
                        consensus_message = binned[6]
                        # cons_alignment = binned[7]

                        # print(f_hq, f_seq)
                        # print(r_hq, r_seq)

                        # if cons_alignment and cons_alignment[0] > 0:
                        #     print()
                        #     print(cons_alignment[4][0] * ' ' + f_seq)
                        #     print(cons_alignment[3][0] * ' ' + r_seq)

                        # print(consensus)
                        # print(consensus_message)
                        # print('----------------------------------------------')

                        str_fhq = ''
                        str_rhq = ''
                        str_cons = ''
                        str_msg = str(consensus_message)
                        str_flq = ''
                        str_rlq = ''

                        if f_hq:
                            # SeqIO.write(sequences=fr, handle=handle_f_hq,
                            #             format='fasta')
                            handle_f_hq.write('>' + f_title + '\n' + f_seq + '\n')
                            str_fhq = str(len(f_seq))
                        else:
                            # SeqIO.write(sequences=fr, handle=handle_f_lq,
                            #             format='fasta')
                            handle_f_lq.write('>' + f_title + '\n' + f_seq + '\n')
                            str_flq = str(len(f_seq))

                        if r_hq:
                            # SeqIO.write(sequences=rr, handle=handle_r_hq,
                            #             format='fasta')
                            handle_r_hq.write('>' + r_title + '\n' + r_seq + '\n')
                            str_rhq = str(len(r_seq))
                        elif r_reads:
                            # SeqIO.write(sequences=rr, handle=handle_r_lq,
                            #             format='fasta')
                            handle_r_lq.write('>' + r_title + '\n' + r_seq + '\n')
                            str_rlq = str(len(r_seq))

                        if f_hq and r_hq and consensus:
                            # SeqIO.write(sequences=consensus,
                            #             handle=handle_all_hq,
                            #             format='fasta')
                            handle_all_hq.write('>' + consensus_title + '\n' +
                                                consensus + '\n')
                            str_cons = str(len(consensus))

                        if f_hq and not consensus:
                            # SeqIO.write(sequences=fr, handle=handle_all_hq,
                            #             format='fasta')
                            handle_all_hq.write('>' + f_title + '\n' + f_seq + '\n')

                        if r_hq and not consensus:
                            # SeqIO.write(sequences=rr, handle=handle_all_hq,
                            #             format='fasta')
                            handle_all_hq.write('>' + r_title + '\n' + r_seq + '\n')

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

                    msg = krother.timestamp() + ' - Sample ' + f['split'][0] + ' done.'
                    print(msg)
                    write_log(msg, lfp)

                    q.task_done()

            for f in file_list:
                if f['name'] != 'Mismatch_f' and f['split'][-1] == 'f':
                    queue.put(f)

            for i in range(cpu):
                worker = Process(target=t, args=(queue,))
                worker.start()
                processes.append(worker)

            queue.join()

            for p in processes:
                p.terminate()

            print()
            write_log('', lfp)

        # Map to Reference
        if commands and ('map_to_reference' in commands):
            if not args.group:
                print('Sample group is required.')
                sys.exit(1)

            if args.group not in groups_map_to_ref_step_dict.keys():
                print("Group '" + args.group + "' does not have any options set in the 'Map to Reference' section of the configuration file.")
                sys.exit(1)

            msg = 'Mapping to reference...\n'
            print(msg)
            write_log(msg, lfp)

            ###KRTEMP

            binned_output_dir_sample = binned_output_dir + args.group + ps
            file_list = krio.parse_directory(binned_output_dir_sample, '_')

            steps = groups_map_to_ref_step_dict[args.group]

            # last_action = ''

            for i, step in enumerate(steps):

                order = step[0]
                program = step[1]
                action = step[2]
                index_path = step[3]

                if i > 0:
                    print()
                    write_log('', lfp)

                msg = 'Step ' + str(order) + ' ' + program + ' ' + action + '\n'
                print(msg)
                write_log(msg, lfp)

                mapped_output_dir_sample = binned_output_dir + args.group + '-' + str(order) + '-' + program + '-' + action + ps
                krio.prepare_directory(mapped_output_dir_sample)

                for f in file_list:
                    if f['ext'] == 'fasta' and f['split'][-1] == 'hq' and f['split'][-2] == 'all':

                        msg = krother.timestamp() + ' - Sample ' + f['split'][0]
                        print(msg)
                        write_log(msg, lfp)

                        bowtie2_out_fasta_option = ''
                        if program == 'bowtie2' and action == 'reject':
                            bowtie2_out_fasta_option = ' --un '
                        elif program == 'bowtie2' and action == 'accept':
                            bowtie2_out_fasta_option = ' --al '

                        subprocess.call(
                            (config.get('General', program+'_executable') +
                                ' --quiet' +
                                ' --threads ' + str(cpu) +
                                # ' --non-deterministic ' +
                                ' --seed ' + config.get('General', 'random_seed') +
                                ' --no-unal' +
                                ' --no-hd' +
                                ' --no-sq' +
                                ' -f' +
                                ' --gbar 2' +
                                ' -k 50' +
                                ' --np 0' +
                                ' --n-ceil L,0,0.5' +
                                ' --end-to-end' +
                                # ' --local' +
                                # ' --ma 0' +
                                ' --mp 6,2' +
                                ' --rdg 6,3' +
                                ' --rfg 6,3' +
                                ' -D 25 -R 3 -N 0 -L 20 -i S,1,0.25' +
                                # ' -D 20 -R 3 -N 0 -L 20 -i S,1,0.50' +
                                ' --met-file ' + mapped_output_dir_sample + f['name'] + '_bowtie2_metrics.tsv' +
                                ' --met 1' +
                                ' -x ' + index_path +
                                ' -U ' + f['path'] +
                                bowtie2_out_fasta_option + mapped_output_dir_sample + f['name'] + '.fasta' +
                                ' -S ' + mapped_output_dir_sample + f['name'] + '.sam'), shell=True)

                file_list = krio.parse_directory(mapped_output_dir_sample, '_')
                # last_action = action

            ###KRTEMPEND

            # # Infer loci using a reference
            # if groups_map_ref_loci_dict[args.group] == True:

            #     msg = '\nInferring loci using a reference...\n'
            #     print(msg)
            #     write_log(msg, lfp)

            #     ###KRTEMP
            #     # last_action = 'accept'
            #     # file_list = krio.parse_directory('/Users/karolis/Dropbox/prj/rad-test/out-05-400/05-binned-fasta/andy-2-bowtie2-accept', '_')
            #     ###KRTEMPEND

            #     if last_action != 'accept':

            #         msg = "In order to infer loci using a reference, last reference mapping action should be 'accept' not 'reject'.\n"
            #         print(msg)
            #         write_log(msg, lfp)

            #     else:

            #         sample_alignments_output_dir_sample = sample_alignments_output_dir + args.group + ps
            #         krio.prepare_directory(sample_alignments_output_dir_sample)

            #         aln_program = config.get('Align Within Samples', 'program')
            #         aln_program_exe = config.get('General', aln_program + '_executable')
            #         aln_program_options = config.get('Align Within Samples', 'options')

            #         processes = list()
            #         queue = JoinableQueue()

            #         def t(q):
            #             while True:
            #                 f = q.get()
            #                 msg = krother.timestamp() + ' - Sample ' + f['split'][0]
            #                 print(msg)
            #                 write_log(msg, lfp)

            #                 # WRITE ALIGNMENTS

            #                 aln_output_file_path = (
            #                     sample_alignments_output_dir_sample +
            #                     f['name'] + '.alignment')
            #                 counts_output_file_path = (
            #                     sample_alignments_output_dir_sample +
            #                     f['name'] + '.counts')

            #                 cluster_depths = krnextgen.alignments_from_sam_file(
            #                     min_seq_cluster=config.getint('General',
            #                                                   'min_seq_cluster'),
            #                     max_seq_cluster=config.getint('General',
            #                                                   'max_seq_cluster'),
            #                     sam_file_path=f['path'],
            #                     aln_output_file_path=aln_output_file_path,
            #                     counts_output_file_path=counts_output_file_path,
            #                     program=aln_program,
            #                     options=aln_program_options,
            #                     program_executable=aln_program_exe)

            #                 # This is the same as in "Produce sample alignments
            #                 # and nucleotide counts per site"

            #                 handle = open((analyzed_samples_output_dir +
            #                                f['split'][0] +
            #                                #  '_' + f['split'][1] +
            #                                '.clusters'), 'wb')
            #                 for c in cluster_depths:
            #                     handle.write(str(c) + '\n')
            #                 handle.close()

            #                 ns = krnextgen.nt_site_counts(counts_output_file_path, 1, 0)
            #                 coverage_list = [sum(x) for x in ns]
            #                 handle = open((analyzed_samples_output_dir +
            #                                f['split'][0] +
            #                                #  '_' + f['split'][1] +
            #                                '.coverage'), 'wb')
            #                 for c in coverage_list:
            #                     handle.write(str(c) + '\n')
            #                 handle.close()

            #                 msg = krother.timestamp() + ' - Sample ' + f['split'][0] + ' done.'
            #                 print(msg)
            #                 write_log(msg, lfp)

            #                 q.task_done()

            #                 ####################################################

            #         for f in file_list:
            #             ###KRTEMP
            #             if f['ext'] == 'sam' and f['split'][-1] == 'hq' and f['split'][-2] == 'all':
            #             # if f['ext'] == 'sam' and f['name'] == '01.S.chilense.Timar_all_hq':
            #             ###KRTEMPEND
            #                 queue.put(f)

            #         threads_local = cpu

            #         for i in range(threads_local):
            #             worker = Process(target=t, args=(queue,))
            #             worker.start()
            #             processes.append(worker)

            #         queue.join()

            #         for p in processes:
            #             p.terminate()

            print()
            write_log('', lfp)

        # Cluster
        if commands and ('cluster_samples' in commands):
            if not args.group:
                print('Sample group is required.')
                sys.exit(1)

            if groups_map_ref_loci_dict[args.group] == True:
                print("Group '" + args.group + "' is set to infer loci using reference.")
                sys.exit(1)

            # if not args.threads:
            #     print('threads is required.')
            #     sys.exit(1)
            # if not args.identity_threshold:
            #     print('identity_threshold is required.')
            #     sys.exit(1)

            clustered_output_dir_sample = clustered_output_dir + args.group + ps

            folders_for_group = list()
            folder_list = krio.parse_directory(
                path=binned_output_dir,
                file_name_sep='-',
                sort='reverse')
            for f in folder_list:
                if f['split'][0] == args.group:
                    folders_for_group.append(f)
            binned_output_dir_sample = folders_for_group[0]['path'] + ps

            krio.prepare_directory(clustered_output_dir_sample)
            file_list = krio.parse_directory(binned_output_dir_sample, '_')

            msg = 'Sorting before clustering...\n'
            print(msg)
            write_log(msg, lfp)

            # Sort files, this will take memory, so we will not parallelize
            # this
            for f in file_list:
                if f['ext'] == 'fasta' and f['split'][-1] == 'hq' and f['split'][-2] == 'all':
                    msg = krother.timestamp() + ' - Sample ' + f['split'][0]
                    print(msg)
                    write_log(msg, lfp)
                    ifp_split = os.path.splitext(f['path'])
                    subprocess.call(
                        (config.get('General', 'usearch_executable') + ' -quiet' +
                            ' -sortbylength ' + f['path'] +
                            ' -output ' + ifp_split[0] + '_sorted' +
                            ifp_split[1]), shell=True)

            # We need to do this again as now we have new (sorted) files in the
            # directory
            file_list = krio.parse_directory(binned_output_dir_sample, '_')

            print()
            write_log('', lfp)

            msg = 'Clustering reads within samples...\n'
            print(msg)
            write_log(msg, lfp)

            processes = list()
            queue = JoinableQueue()

            def t(q):
                while True:
                    f = q.get()
                    # ifp_split = os.path.splitext(f['path'])
                    # input_file_path = ifp_split[0] + '_sorted' + ifp_split[1]

                    msg = krother.timestamp() + ' - Sample ' + f['split'][0] + ' starting...'
                    print(msg)
                    write_log(msg, lfp)

                    first_uc = clustered_output_dir_sample + f['name'] + '_1_1_uc'
                    first_cons_path = clustered_output_dir_sample + f['name'] + '_1_2_cons'

                    krusearch.cluster_file(
                        input_file_path=f['path'],
                        output_file_path=first_uc,
                        identity_threshold=
                        config.getfloat('Cluster Within Samples', 'identity_threshold'),
                        consensus_file_path=first_cons_path,
                        sorted_input=True,
                        algorithm='smallmem',
                        strand='both',
                        threads=1,
                        quiet=True,
                        program=config.get('General', 'usearch_executable'),
                        heuristics=config.getboolean('Cluster Within Samples', 'heuristics'),
                        query_coverage=config.getfloat('Cluster Within Samples', 'query_coverage'),
                        target_coverage=config.getfloat('Cluster Within Samples', 'target_coverage'),
                        sizeout=True,
                        sizein=False,
                        usersort=False
                    )

                    subprocess.call(
                        (config.get('General', 'usearch_executable') + ' -quiet' +
                            ' -sortbysize ' + first_cons_path +
                            ' -output ' + first_cons_path + '_sorted'), shell=True)

                    second_uc = clustered_output_dir_sample + f['name'] + '_2_1_uc'

                    krusearch.cluster_file(
                        input_file_path=first_cons_path + '_sorted',
                        output_file_path=second_uc,
                        identity_threshold=
                        config.getfloat('Cluster Within Samples', 'identity_threshold'),
                        consensus_file_path=False,
                        sorted_input=True,
                        algorithm='smallmem',
                        strand='both',
                        threads=1,
                        quiet=True,
                        program=config.get('General', 'usearch_executable'),
                        heuristics=False,
                        query_coverage=config.getfloat('Cluster Within Samples', 'query_coverage'),
                        target_coverage=config.getfloat('Cluster Within Samples', 'target_coverage'),
                        sizeout=True,
                        sizein=True,
                        usersort=True
                    )

                    # Produce a final cluster file

                    final_uc = clustered_output_dir_sample + f['split'][0] + '.uc'
                    first_dict = krusearch.parse_uc_file(first_uc, 'centroid')
                    second_dict = krusearch.parse_uc_file(second_uc, 'centroid')
                    new_dict = datrie.Trie(string.printable)

                    for k2 in second_dict.keys():
                        cluster2 = second_dict[k2]
                        centroid = unicode(cluster2[0][1].split('=')[1].split(';')[0])
                        # print(centroid)
                        if len(cluster2) > 1:
                            for i2, r2 in enumerate(cluster2):
                                if i2 == 0:
                                    new_dict[centroid] = first_dict[centroid]
                                    # print(i2, centroid)
                                else:
                                    minor_centroid = unicode(r2[1].split('=')[1].split(';')[0])
                                    if r2[0] == '-':
                                        for mr in first_dict[minor_centroid]:
                                            if mr[0] == '-':
                                                mr[0] = '+'
                                            else:
                                                mr[0] = '-'
                                    new_dict[centroid] = new_dict[centroid] + first_dict[minor_centroid]
                                    # print(i2, minor_centroid)
                            # print('--- --- ---')
                        else:
                            new_dict[centroid] = first_dict[centroid]

                    krusearch.write_uc_file(new_dict, final_uc)

                    msg = krother.timestamp() + ' - Sample ' + f['split'][0] + ' done.'
                    print(msg)
                    write_log(msg, lfp)

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

            print()
            write_log('', lfp)

        # Produce sample alignments and nucleotide counts per site
        if commands and ('align_samples' in commands):
            if not args.group:
                print('Sample group is required.')
                sys.exit(1)

            # if groups_map_ref_loci_dict[args.group] == True:
            #     print("Group '" + args.group + "' is set to infer loci using reference.")
            #     sys.exit(1)

            # if not args.min_seq_cluster:
            #     print('min_seq_cluster is required.')
            #     sys.exit(1)
            # if not args.max_seq_cluster:
            #     print('max_seq_cluster is required.')
            #     sys.exit(1)
            # if not args.threads:
            #     print('threads is required.')
            #     sys.exit(1)

            clustered_output_dir_sample = clustered_output_dir + args.group + ps
            sample_alignments_output_dir_sample = sample_alignments_output_dir + args.group + ps

            folders_for_group = list()
            folder_list = krio.parse_directory(
                path=binned_output_dir,
                file_name_sep='-',
                sort='reverse')
            for f in folder_list:
                if f['split'][0] == args.group:
                    folders_for_group.append(f)
            binned_output_dir_sample = folders_for_group[0]['path'] + ps

            krio.prepare_directory(sample_alignments_output_dir_sample)

            msg = ('Producing sample alignments and\n'
                   'read counts per site...\n')
            print(msg)
            write_log(msg, lfp)

            file_list = list()

            # Infer loci using a reference
            if groups_map_ref_loci_dict[args.group] == True:
                msg = 'Inferring loci using a reference...\n'
                print(msg)
                write_log(msg, lfp)

                last_action = binned_output_dir_sample.split('-')[-1].strip(ps)

                if last_action != 'accept':
                    msg = "In order to infer loci using a reference, last reference mapping action should be 'accept' not 'reject'.\n"
                    print(msg)
                    write_log(msg, lfp)
                    sys.exit(1)

                file_list = krio.parse_directory(binned_output_dir_sample, '_')
            else:
                file_list = krio.parse_directory(clustered_output_dir_sample, '_')

            processes = list()
            queue = JoinableQueue()

            aln_program = config.get('Align Within Samples', 'program')
            aln_program_exe = config.get('General', aln_program + '_executable')
            aln_program_options = config.get('Align Within Samples', 'options')

            def t(q):
                while True:
                    f = q.get()
                    fasta_file_path = binned_output_dir_sample + f['name'] + '_all_hq.fasta'
                    aln_output_file_path = (
                        sample_alignments_output_dir_sample +
                        f['name'] + '.alignment')
                    counts_output_file_path = (
                        sample_alignments_output_dir_sample +
                        f['name'] + '.counts')

                    msg = krother.timestamp() + ' - Sample ' + f['split'][0] + ' starting...'
                    print(msg)
                    write_log(msg, lfp)

                    cluster_depths = None

                    if groups_map_ref_loci_dict[args.group] == False:

                        cluster_depths = krnextgen.align_clusters(
                            min_seq_cluster=config.getint('Align Within Samples',
                                                          'min_seq_cluster'),
                            max_seq_cluster=config.getint('Align Within Samples',
                                                          'max_seq_cluster'),
                            uc_file_path=f['path'],
                            fasta_file_path=fasta_file_path,
                            aln_output_file_path=aln_output_file_path,
                            counts_output_file_path=counts_output_file_path,
                            program=aln_program,
                            options=aln_program_options,
                            # options='--retree 1 --thread '+str(cpu)
                            program_executable=aln_program_exe
                        )

                    else:

                        cluster_depths = krnextgen.alignments_from_sam_file(
                            min_seq_cluster=config.getint('Align Within Samples',
                                                          'min_seq_cluster'),
                            # max_seq_cluster=0,
                            max_seq_cluster=config.getint('Align Within Samples',
                                                          'max_seq_cluster'),
                            sam_file_path=f['path'],
                            aln_output_file_path=aln_output_file_path,
                            counts_output_file_path=counts_output_file_path,
                            program=aln_program,
                            options=aln_program_options,
                            program_executable=aln_program_exe
                        )

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

                    msg = krother.timestamp() + ' - Sample ' + f['split'][0] + ' done.'
                    print(msg)
                    write_log(msg, lfp)

                    q.task_done()

            for f in file_list:
                if groups_map_ref_loci_dict[args.group] == False:
                    if f['ext'] == 'uc':
                        queue.put(f)
                else:
                    if f['ext'] == 'sam' and f['split'][-1] == 'hq' and f['split'][-2] == 'all':
                        queue.put(f)

            threads_local = cpu
            # if aln_program == 'muscle':
            #     threads_local = cpu * 2

            for i in range(threads_local):
                worker = Process(target=t, args=(queue,))
                worker.start()
                processes.append(worker)

            queue.join()

            for p in processes:
                p.terminate()

            print()
            write_log('', lfp)

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
            if not args.group:
                print('Sample group is required.')
                sys.exit(1)

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

            sample_alignments_output_dir_sample = sample_alignments_output_dir + args.group + ps

            file_list = krio.parse_directory(sample_alignments_output_dir_sample, '_')

            msg = ('Estimating error rate and within-sample heterozygosity...\n')
            print(msg)
            write_log(msg, lfp)

            manager = Manager()
            results = manager.list()

            processes = list()
            queue = JoinableQueue()

            def t(q):
                while True:
                    q_input = q.get()
                    f = q_input[0]
                    results = q_input[1]
                    msg = krother.timestamp() + ' - Sample ' + f['split'][0] + ' starting...'
                    print(msg)
                    write_log(msg, lfp)
                    p = krnextgen.nt_freq(f['path'])
                    # print(f['split'][0], p)
                    ns = krnextgen.nt_site_counts(
                        f['path'],
                        config.getint('e and pi', 'min_seq_cluster'),
                        config.getint('e and pi', 'max_seq_cluster'))
                    mle = krnextgen.mle_e_and_pi(
                        ns=ns,
                        p=p,
                        e0=config.getfloat('e and pi', 'error_rate_initial'),
                        pi0=config.getfloat('e and pi',
                                            'heterozygosity_initial'))

                    results_dict = dict()
                    # results_dict = datrie.Trie(string.printable)

                    results_dict['sample'] = str(f['split'][0])
                    # results_dict['barcode'] = str(f['split'][1])
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
                    # handle.write('barcode\t' + str(results_dict['barcode']) +
                                 # '\n')
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

                    msg = krother.timestamp() + ' - Sample ' + f['split'][0] + ' done: e=' + str(results_dict['e']) + ', pi=' + str(results_dict['pi'])
                    print(msg)
                    write_log(msg, lfp)

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

            with open(analyzed_samples_output_dir + 'stats.csv', 'wb') as f:
                # fieldnames = ['sample', 'barcode', 'rpc', 'clusters', 'rps',
                #               'sites', 'fA', 'fC', 'fG', 'fT', 'e', 'pi',
                #               'negll']
                fieldnames = ['sample', 'fA', 'fC', 'fG', 'fT', 'e',
                              'pi', 'negll']
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writerow({
                    'sample': 'sample',
                    # 'barcode': 'barcode',
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

            print()
            write_log('', lfp)

        # Call consensus bases
        if commands and ('consensus' in commands):
            if not args.group:
                print('Sample group is required.')
                sys.exit(1)

            group = args.group

            consensus_output_dir_sample = consensus_output_dir + group + ps
            sample_alignments_output_dir_sample = sample_alignments_output_dir + group + ps

            krio.prepare_directory(consensus_output_dir_sample)
            file_list = krio.parse_directory(sample_alignments_output_dir_sample, '_')

            msg = 'Calling consensus bases...\n'
            print(msg)
            write_log(msg, lfp)

            processes = list()
            queue = JoinableQueue()

            # Values are tuples: (e, pi)
            group_stats = dict()
            # group_stats = datrie.Trie(string.printable)
            sample_stats = dict()
            # sample_stats = datrie.Trie(string.printable)

#             for group in sample_groups_dict.keys():
            group_samples = sample_groups_dict[group]
            e_list = list()
            h_list = list()
            for sample in group_samples:
                stats_path = analyzed_samples_output_dir + str(sample) + '.stats'
                samples_stats_handle = open(stats_path, 'rb')
                lines = samples_stats_handle.readlines()
                e = float("inf")
                pi = 0.0
                for l in lines:
                    if l.startswith('e'):
                        e = float(l.split('\t')[1])
                        e_list.append(e)
                    if l.startswith('pi'):
                        pi = float(l.split('\t')[1])
                        h_list.append(pi)
                sample_stats[sample] = (e, pi)
            mean_error = sum(e_list) / float(len(e_list))
            mean_heter = sum(h_list) / float(len(h_list))
            group_stats[group] = (mean_error, mean_heter)
            # print(group)
            # print('e', mean_error)
            # print('pi', mean_heter)

            use_mean_e_and_pi = config.getboolean('Consensus', 'use_mean_e_and_pi')

            # Print log messages
            if use_mean_e_and_pi:
                msg = 'Using mean group error rate and within-sample heterozygosity:\n\n'
#                 for group in sample_groups_dict.keys():
                error = group_stats[group][0]
                heter = group_stats[group][1]
                msg = msg + group + ' e=' + str(error) + ', pi=' + str(heter) + '\n'
                print(msg)
                write_log(msg, lfp)
            else:
                msg = 'Using per-sample error rate and within-sample heterozygosity:\n\n'
#                 for group in sample_groups_dict.keys():
                for sample in sample_groups_dict[group]:
                    error = sample_stats[sample][0]
                    heter = sample_stats[sample][1]
                    msg = msg + sample + ' e=' + str(error) + ', pi=' + str(heter) + '\n'
                print(msg)
                write_log(msg, lfp)
            # End print log messages

            def t(q):
                while True:
                    q_input = q.get()
                    f = q_input[0]
                    group_stats = q_input[1]
                    sample_stats = q_input[2]
                    use_mean_e_and_pi = q_input[3]

                    msg = krother.timestamp() + ' - Sample ' + f['split'][0] + ' starting...'
                    print(msg)
                    write_log(msg, lfp)

                    ns = krnextgen.nt_site_counts(
                        f['path'],
                        1,
                        0,
                        rettype='dict')

                    handle = open((consensus_output_dir_sample +
                                   f['split'][0] +
                                   '_' +
                                   'consensus' +
                                   '.fasta'), 'wb')

                    current_sample = f['split'][0]
                    error = float("inf")
                    heter = 0.0
                    if use_mean_e_and_pi:
#                         for group in sample_groups_dict.keys():
                        group_samples = sample_groups_dict[group]
                        if current_sample in group_samples:
                            error = group_stats[group][0]
                            heter = group_stats[group][1]
#                             break
                    else:
                        error = sample_stats[current_sample][0]
                        heter = sample_stats[current_sample][1]

                    # print('e', error)
                    # print('pi', heter)

                    threshold_probability = config.getfloat(
                        'Consensus', 'threshold_probability')
                    low_quality_residue = config.get(
                        'General', 'low_quality_residue')
                    min_seq_cluster = config.getint('Consensus', 'min_seq_cluster')
                    max_seq_cluster = config.getint('Consensus', 'max_seq_cluster')

                    for k in ns.keys():
                        handle.write('>' + k + '\n')
                        sites = ns[k]
                        sequence = ''
                        for site in sites:
                            c = krnextgen.consensus_base(
                                site,
                                e=error,
                                pi=heter,
                                p=threshold_probability,
                                low_quality_residue=low_quality_residue,
                                min_total_per_site=min_seq_cluster,
                                max_total_per_site=max_seq_cluster)
                            # print(site, c)
                            sequence = sequence + c[3]
                        handle.write(sequence + '\n')

                    handle.close()

                    msg = krother.timestamp() + ' - Sample ' + f['split'][0] + ' done.'
                    print(msg)
                    write_log(msg, lfp)

                    q.task_done()

            for f in file_list:
                if f['ext'] == 'counts':
                    queue.put((f, group_stats, sample_stats, use_mean_e_and_pi))

            for i in range(cpu):
                worker = Process(target=t, args=(queue,))
                worker.start()
                processes.append(worker)

            queue.join()

            for p in processes:
                p.terminate()

            # Create a master consensus file
            msg = '\nCreating grouped-consensus files...'
            print(msg)
            write_log(msg, lfp)
            krio.prepare_directory(grouped_consensus_output_dir)
            # file_list = krio.parse_directory(consensus_output_dir, '_')

            iupac = kriupac.IUPAC_DOUBLE_DNA_DICT

            low_quality_residue = config.get('General', 'low_quality_residue')
            min_read_length = config.getint('Consensus', 'min_read_length')
            max_prop_low_quality_sites = config.getfloat(
                'Consensus', 'max_prop_low_quality_sites')

#             for group in sample_groups_dict.keys():

            group_samples = sample_groups_dict[group]

            consensus_handle = open((grouped_consensus_output_dir + group
                                    + '_consensus' + '.fasta'), 'wb')

            consensus_handle_masked = open((grouped_consensus_output_dir
                                           + group + '_consensus_masked'
                                           + '.fasta'), 'wb')

            # for f in file_list:
            for sample in group_samples:
                # sample = f['split'][0]
                # f_handle = open(f['path'], 'rb')
                f_handle = open(consensus_output_dir_sample + sample + '_consensus.fasta', 'rb')
                lines = f_handle.readlines()
                label = ''
                for l in lines:
                    if l.startswith('>'):
                        label = l.split('>')[1]
                    else:
                        l = l.strip()
                        seq = l.strip(low_quality_residue)

                        prop_lq = krnextgen.proportion_low_quality_sites(
                            seq, low_quality_residue=low_quality_residue)

                        # print(str(len(seq)) + ' ' + str(prop_lq) + ' ' + seq)

                        if len(seq) >= min_read_length and prop_lq <= max_prop_low_quality_sites:
                            consensus_handle.write('>' + sample + '_' + label)
                            consensus_handle.write(seq + '\n')

                            consensus_handle_masked.write('>' + sample + '_' + label)
                            # seq = re.sub('[RYMKWS]', low_quality_residue, seq)
                            # print(seq)
                            for k in iupac.keys():
                                for i in range(0, seq.count(iupac[k])):
                                    # rand = binom.rvs(1, 0.5)
                                    # rand = binomial(1, 0.5)
                                    rand = numpy.random.randint(0, 2)
                                    seq = seq.replace(iupac[k], k[rand], 1)
                            # print(seq)
                            consensus_handle_masked.write(seq + '\n')

            consensus_handle.close()
            consensus_handle_masked.close()

            print()
            write_log('', lfp)

        # Cluster sequences between samples
        if commands and ('cluster_between' in commands):
            if not args.group:
                print('Sample group is required.')
                sys.exit(1)

            krio.prepare_directory(between_sample_clusters_output_dir)
            file_list = krio.parse_directory(grouped_consensus_output_dir, '_')

            msg = 'Clustering loci between samples...\n'
            print(msg)
            write_log(msg, lfp)

            for f in file_list:
                if f['split'][-1] != 'masked':
                    continue

                if not f['name'].startswith(args.group):
                    continue

                # f_path = consensus_output_dir + 'consensus.fasta'
                msg = krother.timestamp() + ' - ' + f['split'][0]
                print(msg)
                write_log(msg, lfp)
                f_path = f['path']
                ifp_split = os.path.splitext(f_path)
                f_path_sorted = ifp_split[0] + '_sorted' + ifp_split[1]
                subprocess.call((config.get('General', 'usearch_executable') + ' -quiet' + ' -sortbylength ' + f_path + ' -output ' + f_path_sorted), shell=True)

                # Put outgroups at the end of file
                f_path_outgroup_sorted = ifp_split[0] + '_sorted_outgroup' + ifp_split[1]
                f_path_sorted_handle = open(f_path_sorted, "rU")
                f_outgroup_sorted_handle = open(f_path_outgroup_sorted, "w")

                outgroup_lines = list()
                next_line_out = False
                next_line_in = False
                for l in f_path_sorted_handle:
                    if next_line_out and not l.startswith('>'):
                        outgroup_lines.append(l)
                        continue
                    elif next_line_in and not l.startswith('>'):
                        f_outgroup_sorted_handle.write(l)
                        continue
                    else:
                        next_line_out = False
                        next_line_in = False
                    for o in outgroups:
                        if l.startswith('>' + o):
                            outgroup_lines.append(l)
                            next_line_out = True
                            break
                    if not next_line_out:
                        f_outgroup_sorted_handle.write(l)
                        next_line_in = True

                f_outgroup_sorted_handle.writelines(outgroup_lines)

                f_outgroup_sorted_handle.close()
                f_path_sorted_handle.close()

                # Cluster
                output_file = between_sample_clusters_output_dir + f['split'][0] + '_clustered' + '.uc'
                krusearch.cluster_file(
                    input_file_path=f_path_outgroup_sorted,
                    output_file_path=output_file,
                    identity_threshold=
                    config.getfloat('Cluster Between Samples', 'identity_threshold'),
                    sorted_input=True,
                    algorithm='smallmem',
                    strand='both',
                    threads=1,
                    quiet=True,
                    program=config.get('General', 'usearch_executable'),
                    heuristics=config.getboolean('Cluster Between Samples', 'heuristics'),
                    query_coverage=config.getfloat('Cluster Between Samples', 'query_coverage'),
                    target_coverage=config.getfloat('Cluster Between Samples', 'target_coverage'),
                    sizein=False,
                    sizeout=False,
                    usersort=True
                )

            print()
            write_log('', lfp)

        # Align loci between samples
        if commands and ('align_between' in commands):
            if not args.group:
                print('Sample group is required.')
                sys.exit(1)

            min_seq_locus = config.getint('Align Between Samples', 'min_seq_locus')
            max_hetero_per_column = config.getint('Align Between Samples', 'max_het_per_column')
            between_sample_alignments_output_dir = between_sample_alignments_output_dir.rstrip(ps) + '-min_seq_locus_' + str(min_seq_locus) + '-max_hetero_per_column_' + str(max_hetero_per_column) + ps

            krio.prepare_directory(between_sample_alignments_output_dir)
            file_list = krio.parse_directory(between_sample_clusters_output_dir, '_')

            msg = 'Aligning loci between samples...\n'
            print(msg)
            write_log(msg, lfp)

            aln_program = config.get('Align Between Samples', 'program')
            aln_program_exe = config.get('General', aln_program + '_executable')
            aln_program_options = config.get('Align Between Samples', 'options')

            for f in file_list:

                if not f['name'].startswith(args.group):
                    continue

                msg = krother.timestamp() + ' - ' + f['split'][0]
                print(msg)
                write_log(msg, lfp)

                fasta_file_path = grouped_consensus_output_dir + f['split'][0] + '_consensus.fasta'
                aln_clustal_phylip_file_path = between_sample_alignments_output_dir + f['split'][0] + '.phy'
                aln_output_file_path = between_sample_alignments_output_dir + f['split'][0] + '.alignment'
                counts_output_file_path = between_sample_alignments_output_dir + f['split'][0] + '.counts'

                cluster_depths = krnextgen.align_clusters(
                    min_seq_cluster=min_seq_locus,
                    max_seq_cluster=0,
                    uc_file_path=f['path'],
                    fasta_file_path=fasta_file_path,
                    aln_clustal_phylip_file_path=aln_clustal_phylip_file_path,
                    aln_output_file_path=aln_output_file_path,
                    counts_output_file_path=counts_output_file_path,
                    program=aln_program,
                    options=aln_program_options,
                    program_executable=aln_program_exe
                )

            file_list = krio.parse_directory(between_sample_alignments_output_dir, '_')

            msg = '\nCreating concatenated alignment...\n'
            print(msg)
            write_log(msg, lfp)

            iupac_double_dna_set = set(kriupac.IUPAC_DOUBLE_DNA_DICT.values())
            max_alignment_length = config.getint('Align Between Samples', 'max_alignment_length')

            for f in file_list:
                # print(f['full'])
                if f['ext'] == 'phy' and len(f['split']) == 1:
                    msg = f['split'][0] + '\n'
                    print(msg)
                    write_log(msg, lfp)
                    all_aln = list(AlignIO.parse(f['path'], "phylip-relaxed"))
                    # all_aln = list()
                    # for aln in all_aln_generator:
                    #     all_aln.append(aln)

                    # Filter loci that have too many heterozygosities per column
                    # or are too long
                    number_accepted = 0
                    all_aln_filtered = list()
                    number_of_aln = len(all_aln)
                    for aln in all_aln:
                        # krcl.print_progress(current=cur, total=number_of_aln, length=50, prefix='\t')
                        accept = True
                        column_count = aln.get_alignment_length()
                        if column_count > max_alignment_length:
                            accept = False
                            # print('Reject: ' + str(column_count) + ' > ' + str(max_alignment_length))
                            continue
                        for column in range(0, column_count):
                            counts = dict()
                            # counts = datrie.Trie(string.printable)
                            for c in aln[:, column]:
                                c = c.upper()
                                if c in iupac_double_dna_set:
                                    counts[c] = counts.get(c, 0) + 1
                            if sum(counts.values()) > max_hetero_per_column:
                                accept = False
                                # print('Reject:', counts)
                                break
                        if accept:
                            all_aln_filtered.append(aln)
                            number_accepted = number_accepted + 1

                    msg = 'Accepted ' + str(number_accepted) + ' out of ' + str(number_of_aln) + ' loci.\n'
                    print(msg)
                    write_log(msg, lfp)
                    ############################################################

                    partitions_output_file = between_sample_alignments_output_dir + f['split'][0] + '_concatenated_partitions.csv'
                    raxml_partitions_output_file = between_sample_alignments_output_dir + f['split'][0] + '_concatenated_partitions_raxml'
                    f_part = open(partitions_output_file, 'wb')
                    f_part_raxml = open(raxml_partitions_output_file, 'wb')

                    concatenated = kralign.concatenate(all_aln_filtered, 10)

                    concatenated_aln = concatenated[0]
                    cat_partitions = concatenated[1]

                    f_part.write('locus,start,end\n')
                    for i, part in enumerate(cat_partitions):
                        raxml_part_line = 'DNA, ' + str(i) + ' = ' + str(part[0]) + '-' + str(part[1]) + '\n'
                        f_part_raxml.write(raxml_part_line)
                        part_line = str(i) + ',' + str(part[0]) + ',' + str(part[1]) + '\n'
                        f_part.write(part_line)

                    cat_aln_output_file_path = between_sample_alignments_output_dir + f['split'][0] + '_concatenated.phy'
                    AlignIO.write(concatenated_aln, cat_aln_output_file_path, "phylip-relaxed")

                    f_part.close()
                    f_part_raxml.close()

                    msg = '\nProducing RAxML commands...\n'
                    print(msg)
                    write_log(msg, lfp)

                    raxml_commands_file = between_sample_alignments_output_dir + f['split'][0] + '_raxml_commands.txt'
                    f_raxml = open(raxml_commands_file, 'wb')

                    raxml_line_1 = 'raxml \\\n'
                    f_raxml.write(raxml_line_1)

                    raxml_line_2 = '-s ' + cat_aln_output_file_path + ' \\\n'
                    f_raxml.write(raxml_line_2)

                    # raxml_line_3 = '-q ' + raxml_partitions_output_file + ' \\\n'
                    # f_raxml.write(raxml_line_3)

                    raxml_output_dir = between_sample_alignments_output_dir + f['split'][0] + '_raxml'
                    krio.prepare_directory(raxml_output_dir)

                    raxml_line_4 = '-w ' + raxml_output_dir + ' \\\n'
                    f_raxml.write(raxml_line_4)

                    outgroups_for_tree = list()
                    for a in concatenated_aln:
                        if a.id in outgroups:
                            outgroups_for_tree.append(a.id)

                    raxml_line_5 = '-o ' + ','.join(outgroups_for_tree) + ' \\\n'
                    f_raxml.write(raxml_line_5)

                    raxml_line_6 = '-m GTRCAT \\\n-j \\\n-T 4 \\\n-N 1 \\\n'
                    f_raxml.write(raxml_line_6)

                    raxml_line_7 = '-p ' + str(random.randrange(0, 1000000000)) + ' \\\n'
                    f_raxml.write(raxml_line_7)

                    raxml_line_8 = '-n ' + f['split'][0] + '\n'
                    f_raxml.write(raxml_line_8)

                    f_raxml.close()

        msg = 'Started: ' + str(start_time).split('.')[0]
        print(msg)
        write_log(msg, lfp)

        msg = 'Ended: ' + krother.timestamp()
        print(msg)
        write_log(msg, lfp)
