# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2022 Vinh Tran
#
#  This script is used to list all available taxa of the installed fdog
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
from ete3 import NCBITaxa


def getNcbiName(taxonName):
    ncbi = NCBITaxa()
    taxId = taxonName.split('@')[1]
    try:
        name = ncbi.get_taxid_translator([taxId])[int(taxId)]
    except:
        name = taxonName
    return(name)


def getTaxa():
    # get data path
    fdogPath = os.path.realpath(__file__).replace('/showTaxa.py','')
    pathconfigFile = fdogPath + '/bin/pathconfig.txt'
    if not os.path.exists(pathconfigFile):
        sys.exit('No pathconfig.txt found. Please run fdog.setup (https://github.com/BIONF/fDOG/wiki/Installation#setup-fdog).')
    with open(pathconfigFile) as f:
        dataPath = f.readline().strip()

    # print taxa in coreTaxa_dir
    print('##### Data found at %s' % dataPath)
    print('\n##### Taxa in the core sets, which can be used as reference species #####\n')
    for taxon in sorted(os.listdir(dataPath + '/coreTaxa_dir/')):
        if os.path.isdir(dataPath + '/coreTaxa_dir/' + taxon):
            print('%s\t%s' % (taxon, getNcbiName(taxon)))

    # print taxa in searchTaxa_dir
    print('\n##### Search taxa. in which you can search orthologs #####\n')
    for taxon in sorted(os.listdir(dataPath + '/searchTaxa_dir/')):
        if os.path.isdir(dataPath + '/searchTaxa_dir/' + taxon):
            print('%s\t%s' % (taxon, getNcbiName(taxon)))


def main():
    getTaxa()

if __name__ == '__main__':
    main()
