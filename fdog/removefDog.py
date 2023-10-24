# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2022 Vinh Tran
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
import shutil
from pkg_resources import get_distribution

import fdog.setupfDog as setupfDog_fn

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
    version = get_distribution('fdog').version
    parser = argparse.ArgumentParser(description='You are running fDOG version ' + str(version) + '.')
    parser.add_argument('--all', help='Remove fdog together with all files/data within the installed fdog directory', action='store_true', default=False)
    args = parser.parse_args()
    data = args.all

    fdogPath = os.path.realpath(__file__).replace('/removefDog.py','')
    dataPath = setupfDog_fn.get_data_path(fdogPath)

    if data:
        print('All files and folders in %s will be removed! Enter to continue' % fdogPath)
    else:
        print('fdog will be uninstalled. Some files/data still can be found in %s! Enter to continue' % fdogPath)
    if query_yes_no('Are you sure?'):
        if data:
            folders = ['bin', 'data']
            for f in folders:
                dirPath = fdogPath+'/'+f
                if os.path.exists(os.path.abspath(dirPath)):
                    print('removing %s...' % f)
                    shutil.rmtree(dirPath)
        uninstallCmd = 'pip uninstall fdog'
        try:
            subprocess.call([uninstallCmd], shell = True)
        except:
            print('Error by uninstalling fdog. Please manually uninstall it using <pip uninstall fdog>')
        if os.path.exists(os.path.abspath(fdogPath)):
            shutil.rmtree(fdogPath)

    print('NOTE: fdog data are still available at\n%s.' % dataPath)


if __name__ == '__main__':
    main()
