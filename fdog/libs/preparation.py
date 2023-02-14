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
    """ Getting path to hmm core set, coreTaxa_dir, searchTaxa_dir and annotation_dir"""
    (pathFile, outpath, hmmpath, corepath, searchpath, annopath) = args
    ### get fdog and data path
    data_path = ''
    fdog_path = os.path.realpath(__file__).replace('/libs/preparation.py','')
    pathconfigFile = fdog_path + '/bin/pathconfig.yml'
    if not os.path.exists(pathconfigFile):
        sys.exit(
            f'No pathconfig.txt found at {pathconfigFile}. Please run fdog.setup '
            + '(https://github.com/BIONF/fDOG/wiki/Installation#setup-fdog).')

    if pathFile:
        pathconfigFile = os.path.abspath(pathFile)

    cfg = general_fn.load_config(pathconfigFile)
    try:
        data_path = cfg['dataPath']
    except:
        data_path = os.getcwd()

    if hmmpath == '':
        hmmpath = outpath + '/core_orthologs'
        Path(hmmpath).mkdir(parents = True, exist_ok = True)

    if corepath == '':
        try:
            corepath = cfg['corepath']
        except:
            corepath = data_path + '/coreTaxa_dir'
        general_fn.check_file_exist(corepath)
    if searchpath == '':
        try:
            searchpath = cfg['searchpath']
        except:
            searchpath = data_path + '/searchTaxa_dir'
        general_fn.check_file_exist(searchpath)
    if annopath == '':
        try:
            annopath = cfg['annopath']
        except:
            annopath = data_path + '/annotation_dir'
        general_fn.check_file_exist(annopath)
    return(hmmpath, corepath, searchpath, annopath)


def check_input(args):
    (seqFile, refspec, outpath, hmmpath, corepath,
        searchpath, annopath, pathFile) = args
    fdog_path = os.path.realpath(__file__).replace('/libs/preparation.py','')
    # create output directory
    Path(outpath).mkdir(parents = True, exist_ok = True)
    Path(hmmpath).mkdir(parents = True, exist_ok = True)
    # check path existing
    hmmpath, corepath, searchpath, annopath = parsing_paths(
        [pathFile, outpath, hmmpath, corepath, searchpath, annopath])
    for path in [hmmpath, corepath, searchpath, annopath]:
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
    if not os.path.exists(os.path.abspath(corepath+'/'+refspec)):
        exit('Reference taxon %s not found in %s' % (refspec, corepath))
    return (seqFile, hmmpath, corepath, searchpath, annopath)


def get_seed_id_from_fa(core_fa, refspec):
    """ Get seed ID from core ortholog fasta file
    (used if --reuseCore option is specified)
    """
    core_seqs = SeqIO.to_dict((SeqIO.parse(open(core_fa), 'fasta')))
    core_ids = core_seqs.keys()
    seed_id = [s for s in core_ids if refspec in s][0].split('|')[-1]
    return(seed_id)


def identify_seed_id(seqFile, refspec, corepath, debug, silentOff):
    """ Identify seed ID in reference protein set using BLAST """
    refspec_db = '%s/%s/%s' % (corepath, refspec, refspec)
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
    blast_xml = blast_fn.do_blastsearch(seqFile, refspec_db, evalBlast = 0.001)
    blast_out = blast_fn.parse_blast_xml(blast_xml)
    for hit in blast_out['hits']:
        if blast_out['hits'][hit]['align_len'] == blast_out['query_len']:
            return(hit)
        elif abs(int(blast_out['hits'][hit]['align_len']) - int(blast_out['query_len'])) < 10:
            output_fn.print_stdout(silentOff, 'WARNING: Found seed sequence shorter than input!')
            return(hit)
        else:
            sys.exit('ERROR: Cannot find seed sequence in genome of reference species for %s!' % blast_out['query'])
