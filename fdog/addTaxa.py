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
from Bio import SeqIO
import multiprocessing as mp
from tqdm import tqdm
from ete3 import NCBITaxa
import re
import shutil
from datetime import datetime
import time
from pkg_resources import get_distribution
from collections import OrderedDict

import fdog.libs.zzz as general_fn
import fdog.libs.tree as tree_fn
import fdog.libs.addtaxon as add_taxon_fn


def parse_map_file(mapping_file, folIn):
    """ Create spec name from mapping file
    And also check if given input files in mapping file exist
    """
    name_dict = {}
    with open(mapping_file) as f:
      for line in f:
        if not '#' in line:
            tmp = line.split('\t')
            file_name = tmp[0]
            file_in = '%s/%s' % (folIn, file_name)
            general_fn.check_file_exist(file_in)
            tax_id = tmp[1].strip()
            try:
                tax_name = tmp[2].strip()
            except:
                tax_name = ''
            try:
                ver = tmp[3].strip()
            except:
                ver = datetime.today().strftime('%y%m%d')
            spec_name = add_taxon_fn.generate_spec_name(tax_id, tax_name, ver)
            name_dict[file_in] = spec_name
    return(name_dict)


def main():
    version = get_distribution('fdog').version
    parser = argparse.ArgumentParser(description='You are running fDOG version ' + str(version) + '.')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i', '--input', help='Path to input folder', action='store', default='', required=True)
    required.add_argument('-m', '--mapping',
                            help='Tab-delimited text file containing <fasta_file_name>tab<taxonID>tab<taxonName>tab<genome_version>. The last 2 columns are optional.',
                            action='store', default='', required=True)
    optional.add_argument('-o', '--outPath', help='Path to output directory', action='store', default='')
    optional.add_argument('--searchpath', help='Path to search taxa folder (e.g. fdog_data/searchTaxa_dir)', action='store', default='')
    optional.add_argument('--corepath', help='Path to core taxa folder (e.g. fdog_data/coreTaxa_dir)', action='store', default='')
    optional.add_argument('--annopath', help='Path to annotation folder (e.g. fdog_data/annotation_dir)', action='store', default='')
    optional.add_argument('-c', '--coreTaxa', help='Include this taxon to core taxa (i.e. taxa in coreTaxa_dir folder)', action='store_true', default=False)
    optional.add_argument('-a', '--noAnno', help='Do NOT annotate this taxon using fas.doAnno', action='store_true', default=False)
    optional.add_argument('--cpus', help='Number of CPUs used for annotation. Default = available cores - 1', action='store', default=0, type=int)
    optional.add_argument('--replace', help='Replace special characters in sequences by "X"', action='store_true', default=False)
    optional.add_argument('--delete', help='Delete special characters in sequences', action='store_true', default=False)
    optional.add_argument('--force', help='Force overwrite existing data', action='store_true', default=False)

    args = parser.parse_args()
    folIn = args.input
    folIn = os.path.abspath(folIn)
    mapping = args.mapping
    general_fn.check_file_exist(mapping)
    outPath = args.outPath
    searchpath = args.searchpath
    corepath = args.corepath
    annopath = args.annopath
    noAnno = args.noAnno
    coreTaxa = args.coreTaxa
    cpus = args.cpus
    if cpus == 0:
        cpus = mp.cpu_count()-2
    replace = args.replace
    delete = args.delete
    add_taxon_fn.check_conflict_opts(replace, delete)
    force = args.force

    start = time.time()
    ### parse mapping file
    name_dict = parse_map_file(mapping, folIn)

    ### initiate paths
    fdogPath = os.path.realpath(__file__).replace('/addTaxa.py','')
    (outPath, searchpath, corepath, annopath) = add_taxon_fn.get_paths(outPath, fdogPath, searchpath, corepath, annopath)
    Path(searchpath).mkdir(parents = True, exist_ok = True)

    ### create file in searchTaxa_dir [and coreTaxa_dir]
    genome_jobs = []
    blast_jobs = []
    for f in name_dict:
        spec_name = name_dict[f]
        ## remove old folder if force is set
        if force == True:
            if os.path.exists('%s/%s' % (searchpath, spec_name)):
                shutil.rmtree('%s/%s' % (searchpath, spec_name))
            if os.path.exists('%s/%s' % (corepath, spec_name)):
                shutil.rmtree('%s/%s' % (corepath, spec_name))
        ## create jobs
        genome_path = '%s/%s' % (searchpath, spec_name)
        Path(genome_path).mkdir(parents = True, exist_ok = True)
        genome_jobs.append([f, genome_path, spec_name, force, replace, delete])
        if coreTaxa:
            genome_file = '%s/%s.fa' % (genome_path, spec_name)
            blast_jobs.append([searchpath, corepath, outPath, spec_name, genome_file, force, True])
    pool = mp.Pool(cpus)

    print('Parsing genome for %s species...' % len(genome_jobs))
    genome_out = []
    for _ in tqdm(pool.imap_unordered(add_taxon_fn.create_genome, genome_jobs),
            total=len(genome_jobs)):
        genome_out.append(_)
    out_msg = 'Output for %s can be found in %s'  % (spec_name, searchpath)
    if len(blast_jobs) > 0:
        print('\nCreating Blast DB for %s species...' % len(blast_jobs))
        blast_out = []
        for _ in tqdm(pool.imap_unordered(add_taxon_fn.create_blastdb, blast_jobs),
                total=len(blast_jobs)):
            blast_out.append(_)
        out_msg = '%s, %s' % (out_msg, corepath)

    ### create annotation
    if not noAnno:
        Path(annopath).mkdir(parents = True, exist_ok = True)
        for f in name_dict:
            genome_file = '%s/%s/%s.fa' % (searchpath, name_dict[f], name_dict[f])
            add_taxon_fn.create_annoFile(annopath, genome_file, cpus, force)
        if os.path.exists('%s/tmp' % annopath):
            shutil.rmtree('%s/tmp' % annopath)
        out_msg = '%s, %s' % (out_msg, annopath)

    end = time.time()
    print('==> Adding %s taxa finished in %s'  % (len(name_dict), '{:5.3f}s'.format(end - start)))
    print('==> %s' % out_msg)

if __name__ == '__main__':
    main()
