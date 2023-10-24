# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2022 Vinh Tran
#
#  This script is used to prepare data for fdog.
#  For each given genome FASTA file, It will create a folder within searchTaxa_dir
#  with the naming scheme of fdog ([Species acronym]@[NCBI ID]@[Proteome version]
#  e.g HUMAN@9606@3), a annotation file in JSON format in annotation_dir and
#  a blast DB in coreTaxa_dir folder (optional).
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
from pathlib import Path
import shutil
import multiprocessing as mp
from datetime import datetime
from pkg_resources import get_distribution

import fdog.libs.zzz as general_fn
import fdog.libs.tree as tree_fn
import fdog.libs.addtaxon as add_taxon_fn


def main():
    version = get_distribution('fdog').version
    parser = argparse.ArgumentParser(description='You are running fDOG version ' + str(version) + '.')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-f', '--fasta', help='FASTA file of input taxon', action='store', default='', required=True)
    required.add_argument('-i', '--taxid', help='Taxonomy ID of input taxon', action='store', default='', required=True, type=int)
    optional.add_argument('-o', '--outPath', help='Path to output directory', action='store', default='')
    optional.add_argument('--searchpath', help='Path to search taxa folder (e.g. fdog_data/searchTaxa_dir)', action='store', default='')
    optional.add_argument('--corepath', help='Path to core taxa folder (e.g. fdog_data/coreTaxa_dir)', action='store', default='')
    optional.add_argument('--annopath', help='Path to annotation folder (e.g. fdog_data/annotation_dir)', action='store', default='')
    optional.add_argument('-n', '--name', help='Acronym name of input taxon', action='store', default='', type=str)
    optional.add_argument('-v', '--verProt', help='Proteome version', action='store', default='', type=str)
    optional.add_argument('-c', '--coreTaxa', help='Include this taxon to core taxa (i.e. taxa in coreTaxa_dir folder)', action='store_true', default=False)
    optional.add_argument('-a', '--noAnno', help='Do NOT annotate this taxon using fas.doAnno', action='store_true', default=False)
    optional.add_argument('--cpus', help='Number of CPUs used for annotation. Default = available cores - 1', action='store', default=0, type=int)
    optional.add_argument('--replace', help='Replace special characters in sequences by "X"', action='store_true', default=False)
    optional.add_argument('--delete', help='Delete special characters in sequences', action='store_true', default=False)
    optional.add_argument('--force', help='Force overwrite existing data', action='store_true', default=False)

    args = parser.parse_args()

    general_fn.check_file_exist(args.fasta)
    faIn = args.fasta
    taxId = str(args.taxid)
    outPath = args.outPath
    searchpath = args.searchpath
    corepath = args.corepath
    annopath = args.annopath
    name = args.name.upper()
    ver = str(args.verProt)
    if ver == '':
        ver = datetime.today().strftime('%y%m%d')
    noAnno = args.noAnno
    coreTaxa = args.coreTaxa
    cpus = args.cpus
    if cpus == 0:
        cpus = mp.cpu_count()-2
    replace = args.replace
    delete = args.delete
    add_taxon_fn.check_conflict_opts(replace, delete)
    force = args.force

    ### species name after fdog naming scheme
    spec_name = add_taxon_fn.generate_spec_name(taxId, name, ver)
    print('Species name\t%s' % spec_name)

    ### get paths
    fdogPath = os.path.realpath(__file__).replace('/addTaxon.py','')
    (outPath, searchpath, corepath, annopath) = add_taxon_fn.get_paths(outPath, fdogPath, searchpath, corepath, annopath)

    ### remove old folder if force is set
    if force == True:
        if os.path.exists('%s/%s' % (searchpath, spec_name)):
            shutil.rmtree('%s/%s' % (searchpath, spec_name))
        if os.path.exists('%s/%s' % (corepath, spec_name)):
            shutil.rmtree('%s/%s' % (corepath, spec_name))

    ### initiate paths
    genome_path = add_taxon_fn.create_folders(searchpath, corepath, annopath, spec_name, coreTaxa, noAnno)

    ### create file in searchTaxa_dir
    print('Parsing FASTA file...')
    genome_file = add_taxon_fn.create_genome([faIn, genome_path, spec_name, force, replace, delete])
    out_msg = 'Output for %s can be found in %s'  % (spec_name, searchpath)

    ### create blast db
    if coreTaxa:
        print('\nCreating Blast DB...')
        add_taxon_fn.create_blastdb([searchpath, corepath, outPath, spec_name, genome_file, force, False])
        out_msg = '%s, %s' % (out_msg, corepath)

    ### create annotation
    if not noAnno:
        add_taxon_fn.create_annoFile(annopath, genome_file, cpus, force)
        if os.path.exists('%s/tmp' % annopath):
            if not os.listdir('%s/tmp' % annopath):
                shutil.rmtree('%s/tmp' % annopath)
        out_msg = '%s, %s' % (out_msg, annopath)

    print('\n==> %s' % out_msg)

if __name__ == '__main__':
    main()
