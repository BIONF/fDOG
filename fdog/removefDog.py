# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2020 Vinh Tran
#
#  This script is used to uninstall fdog and its data
#
#  This script is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License <http://www.gnu.org/licenses/> for
#  more details
#
#  Contact: tran@bio.uni-frankfurt.de
#
#######################################################################

import sys
import os
import argparse
import subprocess
from pathlib import Path
import shutil

def query_yes_no(question, default='yes'):
    valid = {'yes': True, 'y': True, 'ye': True,
             'no': False, 'n': False}
    if default is None:
        prompt = ' [y/n] '
    elif default == 'yes':
        prompt = ' [Y/n] '
    elif default == 'no':
        prompt = ' [y/N] '
    else:
        raise ValueError('invalid default answer: "%s"' % default)
    while True:
        # sys.stdout.write(question + prompt)
        choice = sys.stdin.readline().rstrip().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write('Please respond with "yes" or "no" '
                             '(or "y" or "n").\n')

def main():
    version = '0.0.1'
    parser = argparse.ArgumentParser(description='You are running fdog.remove version ' + str(version) + '.')
    parser.add_argument('--data', help='Remove fdog together with all files/data within the installed fdog directory', action='store_true', default=False)
    args = parser.parse_args()
    data = args.data

    fdogPath = os.path.realpath(__file__).replace('/removefDog.py','')
    pathconfigFile = fdogPath + '/bin/pathconfig.txt'
    if not os.path.exists(pathconfigFile):
        sys.exit('No pathconfig.txt found. Please run fdog.setup (https://github.com/BIONF/fDOG/wiki/Installation#setup-fdog).')
    with open(pathconfigFile) as f:
        dataPath = f.readline().strip()

    if data:
        print('All files and folders in %s will be removed! Enter to continue' % fdogPath)
    else:
        print('fdog will be uninstalled. Some files/data still can be found in %s! Enter to continue' % fdogPath)
    if query_yes_no('Are you sure?'):
        if data:
            folders = ['bin', 'core_orthologs', 'taxonomy', 'data']
            for f in folders:
                dirPath = fdogPath+'/'+f
                if os.path.exists(os.path.abspath(dirPath)):
                    print('removing %s...' % f)
                    shutil.rmtree(dirPath)
        uninstallCmd = 'pip uninstall fdog'
        try:
            subprocess.call([uninstallCmd], shell = True)
        except:
            print('Error by uninstalling fdog. Please manually uninstall it using pip uninstall fdog')
        if data:
            if os.path.exists(os.path.abspath(fdogPath)):
                shutil.rmtree(fdogPath)

    print('NOTE: fdog genome data are still available at %s.' % dataPath)

if __name__ == '__main__':
    main()
