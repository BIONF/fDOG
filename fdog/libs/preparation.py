# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2022 Vinh Tran
#
#  This file is part of fDOG tool https://github.com/BIONF/fDOG
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
from pathlib import Path
from Bio import SeqIO

import fdog.libs.zzz as general_fn
import fdog.libs.fasta as fasta_fn
import fdog.libs.blast as blast_fn
import fdog.libs.output as output_fn


##### FUNCTIONS FOR DATA/INPUT PREPARATION #####

def parsing_paths(args):
    """ Getting path to hmm core set, blast_dir, genome_dir and weight_dir"""
    (pathFile, outpath, hmmpath, blastpath, searchpath, weightpath) = args
    ### get fdog and data path
    data_path = ''
    fdog_path = os.path.realpath(__file__).replace('/libs/preparation.py','')
    pathconfig_file = fdog_path + '/bin/pathconfig.txt'
    # pathconfig_file = '/home/vinh/anaconda3/envs/test_fas/lib/python3.9/site-packages/fdog/bin/pathconfig.txt' ######################################  REMOVE THIS
    if not os.path.exists(pathconfig_file):
        sys.exit(
            'No pathconfig.txt found at %s. Please run fdog.setup '
            + '(https://github.com/BIONF/fDOG/wiki/Installation#setup-fdog).'
            % pathconfig_file)
    if pathFile == '':
        with open(pathconfig_file) as f:
            data_path = f.readline().strip()
    else:
        cfg = general_fn.load_config(pathFile)
        try:
            data_path = cfg['dataPath']
        except:
            data_path = 'config'

    if hmmpath == '':
        hmmpath = outpath + '/core_orthologs'
        Path(hmmpath).mkdir(parents = True, exist_ok = True)

    if blastpath == '':
        blastpath = data_path + '/blast_dir'
        if data_path == 'config':
            try:
                blastpath = cfg['blastpath']
            except:
                sys.exit('blastpath not found in %s' % pathFile)
    if searchpath == '':
        searchpath = data_path + '/genome_dir'
        if data_path == 'config':
            try:
                searchpath = cfg['searchpath']
            except:
                sys.exit('searchpath not found in %s' % pathFile)
    if weightpath == '':
        weightpath = data_path + '/weight_dir'
        if data_path == 'config':
            try:
                weightpath = cfg['weightpath']
            except:
                sys.exit('weightpath not found in %s' % pathFile)
    return(hmmpath, blastpath, searchpath, weightpath)


def check_input(args):
    (seqFile, refspec, outpath, hmmpath, blastpath,
        searchpath, weightpath, pathFile) = args
    fdog_path = os.path.realpath(__file__).replace('/libs/preparation.py','')
    # create output directory
    Path(outpath).mkdir(parents = True, exist_ok = True)
    Path(hmmpath).mkdir(parents = True, exist_ok = True)
    # check path existing
    hmmpath, blastpath, searchpath, weightpath = parsing_paths(
        [pathFile, outpath, hmmpath, blastpath, searchpath, weightpath])
    for path in [hmmpath, blastpath, searchpath, weightpath]:
        general_fn.check_file_exist(path)
    # check for seqFile
    if not os.path.exists(os.path.abspath(seqFile)):
        if not os.path.exists(fdog_path + '/data/' + seqFile):
            sys.exit(
                '%s not found in %s or %s'
                % (seqFile, os.getcwd(), fdog_path + '/data/'))
        else:
            seqFile = fdog_path + '/data/' + seqFile
    else:
        seqFile = os.path.abspath(seqFile)
    # check refspec
    if not os.path.exists(os.path.abspath(blastpath+'/'+refspec)):
        exit('Reference taxon %s not found in %s' % (refspec, blastpath))
    return (seqFile, hmmpath, blastpath, searchpath, weightpath)


def identify_seed_id(seqFile, refspec, blastpath, debug):
    refspec_db = '%s/%s/%s' % (blastpath, refspec, refspec)
    # first check if input seed ID existiert in refspec genome
    refspec_fa = fasta_fn.read_fasta('%s.fa' % refspec_db)
    seed_fa = SeqIO.parse(open(seqFile),'fasta')
    for seed in seed_fa:
        try:
            if len(refspec_fa.fetch(seed.id)) == len(seed.seq):
                return(seed.id)
        except:
            output_fn.print_debug(debug, 'Identify seed ID', 'Input seed ID not found!')
    # otherwise, perform blast search
    blast_xml = blast_fn.do_blastsearch(seqFile, refspec_db)
    blast_out = blast_fn.parse_blast_xml(blast_xml)
    for hit in blast_out['hits']:
        if blast_out['hits'][hit]['align_len'] == blast_out['query_len']:
            return(hit)
        else:
            sys.exit('ERROR: Cannot find seed sequence in genome of reference species!')