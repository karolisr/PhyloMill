#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
# from __future__ import unicode_literals

if __name__ == '__main__':

    import os
    import sys
    import argparse
    import inspect
    import shutil
    import ConfigParser

    from krpy import krncbi

    ps = os.path.sep

    ############################################################################

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-p',
        '--project_dir',
        type=unicode,
        help="Prepares clean project directory if it doesn't exist. Otherwise \
              sets the project directory.")

    parser.add_argument(
        '-c',
        '--commands',
        type=unicode,
        help='Commands to run.')

    args = parser.parse_args()

    ############################################################################

    commands = None
    if args.commands:
        commands = set([x.strip() for x in args.commands.split(',')])

    ############################################################################

    prj_dir_path = None

    # Prepare clean project directory
    if args.project_dir:

        prj_dir_path = args.project_dir.rstrip(ps)

        if os.path.exists(prj_dir_path):

            print('Using project directory at', prj_dir_path)

            prj_dir_path = prj_dir_path + ps

        else:

            print('Creating project directory at', prj_dir_path)

            prj_dir_path = prj_dir_path + ps

            # Script filename
            script_file_path = inspect.getfile(inspect.currentframe())

            # Script directory path
            script_dir_path = os.path.dirname(os.path.abspath(script_file_path))

            prj_template_dir_path = script_dir_path.strip('workflows') + \
                'data' + ps + 'phylo-prj-tempate'

            shutil.copytree(prj_template_dir_path, prj_dir_path, symlinks=False,
                            ignore=None)

            # wd = os.getcwd()
            # os.chdir(prj_dir_path+'03-name_resolution')
            # call(['./get_ncbi_data.sh'], stdout=open(os.devnull, 'wb'))
            # os.chdir(wd)

            sys.exit(0)

    ############################################################################

    if not commands:
        print('No commands given.')
        sys.exit(0)

    ############################################################################

    # Get configuration information
    config_file_path = prj_dir_path + 'config'
    config = ConfigParser.SafeConfigParser(allow_no_value=True)
    config.optionxform=str
    config.read(config_file_path)

    email = config.get('General', 'email')

    # Parse taxa
    tax_temp = config.items('Taxa')
    tax_temp = [x[0] for x in tax_temp]

    tax_ids = list()

    for tax in tax_temp:
        if tax.isalpha():
            tax_id = list(krncbi.esearch(tax, 'taxonomy', email))[0]
            tax_ids.append(str(tax_id))
        else:
            tax_ids.append(tax)

    tax_ncbi_query_strings = list()
    for t in tax_ids:
        tnqs = 'txid' + str(t) + '[Organism]'
        tax_ncbi_query_strings.append(tnqs)
    taxa_query_str = ' OR '.join(tax_ncbi_query_strings)

    # Parse loci
    loci_temp = config.items('Loci')
    loci = dict()
    for l in loci_temp:
        loci[l[0]] = [x.strip() for x in l[1].split(',')]

    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
    ############################################################################
