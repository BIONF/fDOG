# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2022 Vinh Tran
#
#  This script is used to set default data directories for fdog.
#  These include the path to the core taxa (coreTaxa_dir),
#  search taxa (searchTaxa_dir) and FAS annotation json files (annotation_dir).
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

import os
import argparse
from pkg_resources import get_distribution

import fdog.libs.zzz as general_fn
import fdog.checkData as check_data_fn


def set_data_path(searchpath, corepath, annopath, checkOff):
    """ Set default fDOG data path to pathconfig.yml file """
    fdogPath = os.path.realpath(__file__).replace('/setPaths.py','')
    pathconfigFile = fdogPath + '/bin/pathconfig.yml'
    flag = 0
    if os.path.exists(pathconfigFile):
        print('Default fDOG data paths in %s will be overwritten! Enter to continue.' % pathconfigFile)
        if general_fn.query_yes_no(''):
            flag = 1
    else:
        flag = 0

    if not checkOff:
        caution = check_data(searchpath, corepath, annopath)
        if caution == 1:
            print('Check done! Data are ready to use WITH CAUTION! Are you sure to add these paths as default? (Y/N)')
            if general_fn.query_yes_no('', default='no'):
                flag = 1
            else:
                flag = 0
    else:
        print('WARNING: Data will not be checked! Run fdog.checkData if you encounter any problems!')

    if flag == 1:
        with open(pathconfigFile, 'w') as cf:
            cf.write('corepath: \'%s\'\n' % corepath)
            cf.write('searchpath: \'%s\'\n' % searchpath)
            cf.write('annopath: \'%s\'\n' % annopath)
        print('Finished! New data paths have been saved in %s' % pathconfigFile)
    else:
        print(f'{pathconfigFile} remains unchanged!')


def check_data(searchpath, corepath, annopath):
    """ Perform data check """
    caution = check_data_fn.run_check([searchpath, corepath, annopath, False, False, False, False, False])
    return(caution)


def main():
    version = get_distribution('fdog').version
    parser = argparse.ArgumentParser(description='You are running fDOG version ' + str(version) + '.')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('--searchpath', help='Path to search taxa folder (e.g. fdog_data/searchTaxa_dir)', action='store', default='', required=True)
    required.add_argument('--corepath', help='Path to core taxa folder (e.g. fdog_data/coreTaxa_dir)', action='store', default='', required=True)
    required.add_argument('--annopath', help='Path to annotation folder (e.g. fdog_data/annotation_dir)', action='store', default='', required=True)
    optional.add_argument('--checkOff', help='Turn off checking for valid data', action='store_true', default=False)

    args = parser.parse_args()
    searchpath = os.path.abspath(args.searchpath)
    corepath = os.path.abspath(args.corepath)
    annopath = os.path.abspath(args.annopath)
    checkOff = args.checkOff

    set_data_path(searchpath, corepath, annopath, checkOff)

if __name__ == '__main__':
    main()
