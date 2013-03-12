#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
#from __future__ import unicode_literals

if __name__ == '__main__':

    import sys
    import os
    import argparse
    import csv

    from multiprocessing import Process
    from multiprocessing import JoinableQueue
    from multiprocessing import Manager
    import subprocess

    from Bio import SeqIO

    import krio
    import krbioio
    import krpipe
    import krnextgen
    import krusearch

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

    # parser.add_argument('--output_file_format', type=unicode,
    #                     help='Output file format.')

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

    parser.add_argument('--max_prop_low_quality_sites', type=float,
                        help='Maximum proportion of low quality sites in a \
                        read for it to still be considered acceptable.')

    parser.add_argument('--min_overlap', type=int,
                        help='Minimum overlap required between forward and \
                        reverse reads. If this overlap is not reached, the \
                        reads will be concatenated.')

    parser.add_argument('--mmmr_cutoff', type=float,
                        help='When aligning forward and reverse reads to \
                        check if they overlap, this value is used as a cutoff \
                        when deciding if to accept or reject an alignment: \
                        mmmr = match / (match + miss). miss does not include \
                        ignored characters, by default \'N\'.')

    parser.add_argument('--identity_threshold', type=float,
                        help='Identity value for within-sample read \
                        clustering.')

    parser.add_argument('--min_seq_cluster', type=int,
                        help='Minimum number of sequences in a cluster.')

    parser.add_argument('--max_seq_cluster', type=int,
                        help='Maximum number of sequences in a cluster.')

    parser.add_argument('--error_rate_initial', type=float,
                        help='Initial value of error rate (ε) for maximum \
                        likelihood estimation.')

    parser.add_argument('--heterozygosity_initial', type=float,
                        help='Initial value of heterozygosity (π) for maximum \
                        likelihood estimation.')

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
        masked_output_dir = output_dir + '04-masked-fastq'
        binned_output_dir = output_dir + '05-binned-fasta'
        clustered_output_dir = output_dir + '06-clustered'
        sample_alignments_output_dir = output_dir + '07-sample-alignments'
        analyzed_samples_output_dir = output_dir + '08-analyzed-samples'

        # Split FASTQ files ---------------------------------------------------
        if commands and ('split' in commands):

            if not args.forward_file:
                print('File with forward RAD reads is required.')
                sys.exit(1)
            if not args.threads:
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
            if not args.threads:
                print('threads is required.')
                sys.exit(1)
            # if not args.output_file_format:
            #     print('output_file_format is required.')
            #     sys.exit(1)

            masked_output_dir = masked_output_dir.rstrip(ps) + ps
            krio.prepare_directory(masked_output_dir)
            file_list = krpipe.parse_directory(dmltplx_output_dir_combined,
                                               ' ')

            print('\nMasking low quality sites...')

            processes = list()
            queue = JoinableQueue()

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
                        SeqIO.write(
                            sequences=r_masked,
                            handle=handle,
                            format='fastq'
                        )
                    handle.close()
                    q.task_done()

            for f in file_list:
                if 'mismatch' not in f['name']:
                    queue.put(f)

            for i in range(args.threads):
                worker = Process(target=t, args=(queue,))
                worker.start()
                processes.append(worker)

            queue.join()

            for p in processes:
                p.terminate()

        # Bin results by quality and produce forward and reverse consensus
        # sequences
        if commands and ('bin' in commands):
            if not args.max_prop_low_quality_sites:
                print('max_prop_low_quality_sites is required.')
                sys.exit(1)
            if not args.min_overlap:
                print('min_overlap is required.')
                sys.exit(1)
            if not args.mmmr_cutoff:
                print('mmmr_cutoff is required.')
                sys.exit(1)
            if not args.low_quality_residue:
                print('low_quality_residue is required.')
                sys.exit(1)
            if not args.threads:
                print('threads is required.')
                sys.exit(1)

            masked_output_dir = masked_output_dir.rstrip(ps) + ps
            binned_output_dir = binned_output_dir.rstrip(ps) + ps
            krio.prepare_directory(binned_output_dir)
            file_list = krpipe.parse_directory(masked_output_dir, '_')

            print(('\nBinning results by quality and producing\n'
                   'forward and reverse consensus sequences...\n'))

            processes = list()
            queue = JoinableQueue()

            def t(q):
                while True:
                    f = q.get()

                    base_file_name = f['split'][0] + '_' + f['split'][1] + '_'
                    base_file_path = masked_output_dir + base_file_name

                    print(f['split'][0] + '_' + f['split'][1])

                    f_records = SeqIO.parse(base_file_path + 'f.fastq',
                                            'fastq')

                    r_records = None
                    if os.path.exists(base_file_path + 'r.fastq'):
                        r_records = SeqIO.parse(base_file_path + 'r.fastq',
                                                'fastq')

                    ofp_f_hq = (binned_output_dir + base_file_name +
                                'f_hq.fasta')
                    ofp_f_lq = (binned_output_dir + base_file_name +
                                'f_lq.fasta')

                    ofp_r_hq = None
                    ofp_r_lq = None

                    if r_records:
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

                    if r_records:
                        handle_r_hq = open(ofp_r_hq, 'w')
                        handle_r_lq = open(ofp_r_lq, 'w')

                    handle_all_hq = open(ofp_all_hq, 'w')

                    for fr in f_records:
                        rr = None
                        if r_records:
                            rr = r_records.next()

                        binned = krnextgen.bin_reads(
                            f_record=fr,
                            r_record=rr,
                            max_prop_low_quality_sites=
                            args.max_prop_low_quality_sites,
                            min_overlap=args.min_overlap,
                            mmmr_cutoff=args.mmmr_cutoff,
                            low_quality_residue=str(args.low_quality_residue)
                        )

                        f_hq = binned[0]
                        r_hq = binned[1]
                        consensus = binned[2]

                        if f_hq:
                            SeqIO.write(sequences=fr, handle=handle_f_hq,
                                        format='fasta')
                        else:
                            SeqIO.write(sequences=fr, handle=handle_f_lq,
                                        format='fasta')

                        if r_hq:
                            SeqIO.write(sequences=rr, handle=handle_r_hq,
                                        format='fasta')
                        elif r_records:
                            SeqIO.write(sequences=rr, handle=handle_r_lq,
                                        format='fasta')

                        if f_hq and r_hq and consensus:
                            SeqIO.write(sequences=consensus,
                                        handle=handle_all_hq,
                                        format='fasta')

                        if f_hq and not consensus:
                            SeqIO.write(sequences=fr, handle=handle_all_hq,
                                        format='fasta')

                        if r_hq and not consensus:
                            SeqIO.write(sequences=rr, handle=handle_all_hq,
                                        format='fasta')

                    handle_f_hq.close()
                    handle_f_lq.close()

                    if r_records:
                        handle_r_hq.close()
                        handle_r_lq.close()

                    handle_all_hq.close()
                    q.task_done()

            for f in file_list:
                if f['name'] != 'mismatch_f' and f['split'][-1] == 'f':
                    queue.put(f)

            for i in range(args.threads):
                worker = Process(target=t, args=(queue,))
                worker.start()
                processes.append(worker)

            queue.join()

            for p in processes:
                p.terminate()

        # Cluster
        if commands and ('cluster' in commands):
            if not args.threads:
                print('threads is required.')
                sys.exit(1)
            if not args.identity_threshold:
                print('identity_threshold is required.')
                sys.exit(1)

            binned_output_dir = binned_output_dir.rstrip(ps) + ps
            clustered_output_dir = clustered_output_dir.rstrip(ps) + ps
            krio.prepare_directory(clustered_output_dir)
            file_list = krpipe.parse_directory(binned_output_dir, '_')

            print('\nSorting before clustering...\n')

            # Sort files, this will take memory, so we will not parallelize
            # this
            for f in file_list:
                if f['split'][-1] == 'hq' and f['split'][-2] == 'all':
                    print('Sorting file:', f['full'])
                    ifp_split = os.path.splitext(f['path'])
                    subprocess.call(
                        ('usearch6' + ' -quiet' +
                            ' -sortbylength ' + f['path'] +
                            ' -output ' + ifp_split[0] + '_sorted' +
                            ifp_split[1]), shell=True)

            # We need to do this again as now we have new (sorted) files in the
            # directory
            file_list = krpipe.parse_directory(binned_output_dir, '_')

            print('\nClustering reads within samples...\n')

            processes = list()
            queue = JoinableQueue()

            def t(q):
                while True:
                    f = q.get()
                    # ifp_split = os.path.splitext(f['path'])
                    # input_file_path = ifp_split[0] + '_sorted' + ifp_split[1]
                    output_file = clustered_output_dir + f['name'] + '.uc'
                    print('Clustering file:', f['full'])
                    krusearch.cluster_file(
                        input_file_path=f['path'],
                        output_file_path=output_file,
                        identity_threshold=args.identity_threshold,
                        sorted_input=True,
                        algorithm='smallmem',
                        strand='both',
                        threads=1,
                        quiet=True,
                        program='usearch6'
                    )
                    q.task_done()

            for f in file_list:
                if (f['split'][-1] == 'sorted' and
                    f['split'][-2] == 'hq' and
                        f['split'][-3] == 'all'):
                    queue.put(f)

            for i in range(args.threads):
                worker = Process(target=t, args=(queue,))
                worker.start()
                processes.append(worker)

            queue.join()

            for p in processes:
                p.terminate()

        # Produce sample alignments and nucleotide counts per site
        if commands and ('align_samples' in commands):
            if not args.min_seq_cluster:
                print('min_seq_cluster is required.')
                sys.exit(1)
            if not args.max_seq_cluster:
                print('max_seq_cluster is required.')
                sys.exit(1)
            if not args.threads:
                print('threads is required.')
                sys.exit(1)

            analyzed_samples_output_dir = (
                analyzed_samples_output_dir.rstrip(ps) + ps)

            krio.prepare_directory(analyzed_samples_output_dir)

            binned_output_dir = binned_output_dir.rstrip(ps) + ps

            clustered_output_dir = clustered_output_dir.rstrip(ps) + ps

            sample_alignments_output_dir = (
                sample_alignments_output_dir.rstrip(ps) + ps)

            krio.prepare_directory(sample_alignments_output_dir)

            file_list = krpipe.parse_directory(clustered_output_dir, '_')

            print('\nProducing sample alignments and\n'
                  'nucleotide counts per site...\n')

            processes = list()
            queue = JoinableQueue()

            def t(q):
                while True:
                    q_input = q.get()
                    f = q_input[0]
                    q_id = str(q_input[1])
                    fasta_file_path = binned_output_dir + f['name'] + '.fasta'
                    aln_output_file_path = (
                        sample_alignments_output_dir +
                        f['name'] + '.alignment')
                    counts_output_file_path = (
                        sample_alignments_output_dir +
                        f['name'] + '.counts')
                    print('File:', f['full'])
                    cluster_depths = krnextgen.align_clusters(
                        min_seq_cluster=args.min_seq_cluster,
                        max_seq_cluster=args.max_seq_cluster,
                        uc_file_path=f['path'],
                        fasta_file_path=fasta_file_path,
                        aln_output_file_path=aln_output_file_path,
                        counts_output_file_path=counts_output_file_path,
                        temp_dir_path=sample_alignments_output_dir,
                        temp_file_id=q_id,
                        threads=args.threads)

                    handle = open((analyzed_samples_output_dir +
                                   f['split'][0] + '_' +
                                   f['split'][1] +
                                   '_clusters.csv'), 'wb')
                    for c in cluster_depths:
                        handle.write(str(c) + '\n')
                    handle.close()

                    q.task_done()

            for i, f in enumerate(file_list):
                if (f['split'][-1] == 'sorted' and
                    f['split'][-2] == 'hq' and
                        f['split'][-3] == 'all'):
                    queue.put([f, i])

            ### THREADS HARD-CODED ###
            for i in range(1):
                worker = Process(target=t, args=(queue,))
                worker.start()
                processes.append(worker)

            queue.join()

            for p in processes:
                p.terminate()

        # Estimate error rate, heterozygosity. Produce some statistics
        if commands and ('analyze_samples' in commands):
            if not args.error_rate_initial:
                print('error_rate_initial is required.')
                sys.exit(1)
            if not args.heterozygosity_initial:
                print('heterozygosity_initial is required.')
                sys.exit(1)
            if not args.threads:
                print('threads is required.')
                sys.exit(1)

            sample_alignments_output_dir = (
                sample_alignments_output_dir.rstrip(ps) + ps)
            analyzed_samples_output_dir = (
                analyzed_samples_output_dir.rstrip(ps) + ps)
            krio.prepare_directory(analyzed_samples_output_dir)
            file_list = krpipe.parse_directory(
                sample_alignments_output_dir, '_')

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
                    print('File:', f['name'], 'starting...')
                    p = krnextgen.nt_freq(f['path'])
                    print(p)
                    ns = krnextgen.nt_site_counts(
                        f['path'], args.min_seq_cluster)
                    mle = krnextgen.mle_e_and_pi(
                        ns=ns,
                        p=p,
                        e0=args.error_rate_initial,
                        pi0=args.heterozygosity_initial)

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

                    coverage_list = [sum(x) for x in ns]
                    handle = open((analyzed_samples_output_dir +
                                   results_dict['sample'] + '_' +
                                   results_dict['barcode'] +
                                   '_coverage.csv'), 'wb')
                    for c in coverage_list:
                        handle.write(str(c) + '\n')
                    handle.close()
                    rps = float(sum(coverage_list)) / float(len(coverage_list))
                    results_dict['rps'] = rps
                    results_dict['sites'] = len(coverage_list)

                    cluster_list = list()
                    handle = open((analyzed_samples_output_dir +
                                   results_dict['sample'] + '_' +
                                   results_dict['barcode'] +
                                   '_clusters.csv'), 'rb')
                    for c in handle:
                        cluster_list.append(int(c))
                    handle.close()
                    rpc = float(sum(cluster_list)) / float(len(cluster_list))
                    results_dict['rpc'] = rpc
                    results_dict['clusters'] = len(cluster_list)

                    handle = open((analyzed_samples_output_dir +
                                   results_dict['sample'] + '_' +
                                   results_dict['barcode'] +
                                   '_stats.txt'), 'wb')

                    handle.write('sample\t' + str(results_dict['sample']) +
                                 '\n')
                    handle.write('barcode\t' + str(results_dict['barcode']) +
                                 '\n')
                    # reads per cluster (rpc)
                    handle.write('rpc\t' + str(rpc) + '\n')
                    handle.write('clusters\t' +
                                 str(results_dict['clusters']) + '\n')
                    # reads per site (rps)
                    handle.write('rps\t' + str(rps) + '\n')
                    handle.write('sites\t' + str(results_dict['sites']) + '\n')
                    handle.write('fA\t' + str(results_dict['fA']) + '\n')
                    handle.write('fC\t' + str(results_dict['fC']) + '\n')
                    handle.write('fG\t' + str(results_dict['fG']) + '\n')
                    handle.write('fT\t' + str(results_dict['fT']) + '\n')
                    handle.write('e\t' + str(results_dict['e']) + '\n')
                    handle.write('pi\t' + str(results_dict['pi']) + '\n')
                    handle.write('negll\t' + str(results_dict['negll']) + '\n')

                    handle.close()

                    results.append(results_dict)

                    print('File:', f['name'], 'done.')

                    q.task_done()

            for f in file_list:
                if f['ext'] == 'counts':
                    queue.put([f, results])

            for i in range(args.threads):
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
                fieldnames = ['sample', 'barcode', 'rpc', 'clusters', 'rps',
                              'sites', 'fA', 'fC', 'fG', 'fT', 'e', 'pi',
                              'negll']
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writerow({
                    'sample': 'sample',
                    'barcode': 'barcode',
                    'rpc': 'rpc',
                    'clusters': 'clusters',
                    'rps': 'rps',
                    'sites': 'sites',
                    'fA': 'fA',
                    'fC': 'fC',
                    'fG': 'fG',
                    'fT': 'fT',
                    'e': 'e',
                    'pi': 'pi',
                    'negll': 'negll',
                })
                writer.writerows(results_new)


# p = nt_freq('/home/karolis/Dropbox/code/krpy/testdata/nt.counts')
# print(p)
# ns = nt_site_counts('/home/karolis/Dropbox/code/krpy/testdata/nt.counts')
# mle = mle_e_and_pi(ns, p, e0=0.001, pi0=0.001)
# print(mle)

# --commands split,demultiplex,mask,bin,cluster,align_samples,analyze_samples

# time ./rad.py \
# --output_dir /home/karolis/Dropbox/code/test/rad \
# --forward_file /home/karolis/Dropbox/code/krpy/testdata/rad_forward.fastq \
# --reverse_file /home/karolis/Dropbox/code/krpy/testdata/rad_reverse.fastq \
# --threads 6 \
# --barcodes '/home/karolis/Dropbox/code/krpy/testdata/rad_barcodes.tsv' \
# --max_barcode_mismatch_count 1 \
# --trim_barcode \
# --trim_extra 5 \
# --quality_score_treshold 30 \
# --low_quality_residue N \
# --max_prop_low_quality_sites 0.10 \
# --min_overlap 5 \
# --mmmr_cutoff 0.85 \
# --identity_threshold 0.90 \
# --min_seq_cluster 4 \
# --max_seq_cluster 1000 \
# --error_rate_initial 0.001 \
# --heterozygosity_initial 0.001 \
# --commands split,demultiplex,mask,bin,cluster,align_samples,analyze_samples

# time ./rad.py \
# --output_dir /data/gbs-new \
# --forward_file /data/gbs-andy-david/green_dzaya_DZAYA_R1.PF.fastq \
# --reverse_file /data/gbs-andy-david/green_dzaya_DZAYA_R2.PF.fastq \
# --threads 6 \
# --barcodes /data/gbs-andy-david/barcodes.tsv \
# --max_barcode_mismatch_count 1 \
# --trim_barcode \
# --trim_extra 5 \
# --quality_score_treshold 30 \
# --low_quality_residue N \
# --max_prop_low_quality_sites 0.10 \
# --min_overlap 5 \
# --mmmr_cutoff 0.85 \
# --identity_threshold 0.90 \
# --min_seq_cluster 10 \
# --max_seq_cluster 1000 \
# --error_rate_initial 0.0001 \
# --heterozygosity_initial 0.001 \
# --commands analyze_samples
