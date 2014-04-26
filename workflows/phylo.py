#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import unicode_literals

if __name__ == '__main__':

    import os
    import sys
    import argparse
    import inspect
    import shutil

    ps = os.path.sep

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-p',
        '--project_dir',
        type=unicode,
        help="Prepares clean project directory if it doesn't exist. Otherwise \
              sets the project directory.")

    args = parser.parse_args()

    # Prepare clean project directory
    if args.project_dir:

        prj_dir_path = args.project_dir.rstrip(ps)

        if os.path.exists(prj_dir_path):
            print('Directory at', prj_dir_path, 'already exists.')
            sys.exit(0)

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
