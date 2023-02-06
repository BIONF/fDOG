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
import subprocess
from ete3 import NCBITaxa
import re
from datetime import datetime
from collections import OrderedDict

import fdog.libs.zzz as general_fn
import fdog.libs.blast as blast_fn
import fdog.libs.fasta as fasta_fn
import fdog.libs.tree as tree_fn

##### FUNCTIONS RELATED TO ADDING NEW TAXON TO FDOG DATABASE #####

def check_conflict_opts(replace, delete):
    """ Check if both replace and delete option are specified """
    if delete:
        if replace:
            sys.exit('*** ERROR: only one option can be choose between "--replace" and "--delete"')
    if replace:
        if delete:
            sys.exit('*** ERROR: only one option can be choose between "--replace" and "--delete"')


def get_paths(outPath, fdogPath, searchpath, corepath, annopath):
    """ Get path to searchTaxa_dir, coreTaxa_dir and annotation_dir """
    if outPath == '':
        pathconfigFile = fdogPath + '/bin/pathconfig.yml'
        if not os.path.exists(pathconfigFile):
            sys.exit('No pathconfig.yml found. Please run fdog.setup (https://github.com/BIONF/fDOG/wiki/Installation#setup-fdog).')
        cfg = general_fn.load_config(pathconfigFile)
        try:
            outPath = cfg['dataPath']
        except:
            try:
                corepath = cfg['corepath']
            except:
                pass
            try:
                searchpath = cfg['searchpath']
            except:
                pass
            try:
                annopath = cfg['annopath']
            except:
                pass

    outPath = os.path.abspath(outPath)
    if not searchpath:
        searchpath = outPath + '/searchTaxa_dir/'
    searchpath = os.path.abspath(searchpath)
    if not corepath:
        corepath = outPath + '/coreTaxa_dir/'
    corepath = os.path.abspath(corepath)
    if not annopath:
        annopath = outPath + '/annotation_dir/'
    annopath = os.path.abspath(annopath)
    return(outPath, searchpath, corepath, annopath)

def create_folders(searchpath, corepath, annopath, spec_name, coreTaxa, noAnno):
    """ Create searchTaxa_dir, coreTaxa_dir and annotation_dir in output folder """
    Path(searchpath).mkdir(parents = True, exist_ok = True)
    genome_path = '%s/%s' %  (searchpath, spec_name)
    Path(genome_path).mkdir(parents = True, exist_ok = True)
    if coreTaxa:
        Path(corepath).mkdir(parents = True, exist_ok = True)
    if not noAnno:
        Path(annopath).mkdir(parents = True, exist_ok = True)
    return(genome_path)


def generate_spec_name(tax_id, name, ver):
    """ Create species name with the format <ABBR_NAME>@<NCBI_TAXID>@<VER> """
    if name == "":
        ncbi_name = tree_fn.check_tax_id(tax_id)
        name = tree_fn.abbr_ncbi_name(ncbi_name)
    return(name+'@'+tax_id+'@'+ver)


def create_genome(args):
    """ Create fa and fai in searchTaxa_dir """
    (faIn, genome_path, spec_name, force, replace, delete) = args
    ### load fasta seq
    in_seq = SeqIO.to_dict((SeqIO.parse(open(faIn), 'fasta')))
    if not os.path.exists(genome_path):
            Path(genome_path).mkdir(parents = True, exist_ok = True)
    genome_file = '%s/%s.fa' % (genome_path, spec_name)
    if (not os.path.exists(os.path.abspath(genome_file))) or (os.stat(genome_file).st_size == 0) or force:
        f = open(genome_file, 'w')
        pipe = 0
        long_id = 0
        mod_id_index = 0
        id_dict = {} # id_dict[ori_id] = mod_id
        for id in in_seq:
            ori_id = id
            seq = str(in_seq[id].seq)
            ### check if ID contains empty char or pipe
            if ' ' in id:
                sys.exit('\033[91mERROR: Sequence IDs (e.g. %s) must not contain space(s)!\033[0m' % id)
            else:
                if '|' in id:
                    tmp = re.split('[_|]', id)
                    tmp = list(OrderedDict.fromkeys(tmp))
                    pipe = 1
                    id = '_'.join(tmp)
                    if not ori_id in id_dict:
                        id_dict[ori_id] = id
            ### check if id longer than 20 character
            if len(id) > 20:
                long_id = 1
                mod_id_index = mod_id_index + 1
                id = '%s_%s' % (spec_name.split('@')[1], mod_id_index)
                if not ori_id in id_dict:
                    id_dict[ori_id] = id
            ### check if seq contains special characters
            if seq[-1] == '*':
                seq = seq[:-1]
            specialChr = 'no'
            if any(c for c in seq if not c.isalpha()):
                specialChr = 'yes'
            if specialChr == 'yes':
                if replace or delete:
                    if replace:
                        seq = re.sub('[^a-zA-Z]', 'X', seq)
                    if delete:
                        seq = re.sub('[^a-zA-Z]', '', seq)
                else:
                    sys.exit('\033[91mERROR: %s sequence contains special character!\033[0m\nYou can use --replace or --delete to solve it.' % (id))
            f.write('>%s\n%s\n' % (id, seq))
        f.close()
        ### create index file
        fasta_fn.read_fasta(genome_file)
        ### write .checked file
        cf = open(genome_file+'.checked', 'w')
        cf.write(str(datetime.now()))
        cf.close()
        ### write ID mapping file and give warning if ID changed
        if len(id_dict) > 0:
            mapping_file = '%s.mapping' % genome_file
            with open(mapping_file, 'w') as mp:
                for o,n in id_dict.items():
                    mp.write('%s\t%s\n' % (o,n))
            if pipe == 1:
                print('\033[94mWARNING: Sequence IDs contain pipe(s). They will be replaced by "_"!\033[0m')
            if long_id == 'yes':
                print('\033[94mWARNING: Some headers longer than 80 characters have been automatically shortened. Please check the %s.mapping file for details!\033[0m' % genome_file)
            if pipe == 1:
                print('\033[94mWARNING: Please check the %s file for details!\033[0m' % mapping_file)
    else:
        print(genome_path + '/' + spec_name + '.fa already exists!')
    return(genome_file)


def create_blastdb(args):
    """ Create blastdb for a given fasta genome_file """
    (searchpath, corepath, outPath, spec_name, genome_file, force, silent) = args
    blast_path = '%s/%s' % (corepath, spec_name)
    if (not os.path.exists(os.path.abspath('%s/%s.phr' % (blast_path, spec_name)))) or force:
        blast_fn.make_blastdb([spec_name, genome_file, outPath, corepath, searchpath, silent])
        ### make symlink to fasta files
        fa_in_genome = "%s/%s/%s.fa" % (searchpath, spec_name, spec_name)
        fai_in_genome = "%s/%s/%s.fa.fai" % (searchpath, spec_name, spec_name)
        fa_in_blast = "%s/%s.fa" % (blast_path, spec_name)
        fai_in_blast = "%s/%s.fa.fai" % (blast_path, spec_name)
        if not os.path.exists(fa_in_blast):
            os.symlink(fa_in_genome, fa_in_blast)
        if not os.path.exists(fai_in_blast):
            os.symlink(fai_in_genome, fai_in_blast)
    else:
        print('Blast DB already exists!')


def create_annoFile(annopath, genome_file, cpus, force):
    """ Create annotation json for a given genome_file """
    annoCmd = 'fas.doAnno -i %s -o %s --cpus %s' % (genome_file, annopath, cpus)
    if force:
        annoCmd = annoCmd + " --force"
    try:
        subprocess.call([annoCmd], shell = True)
    except:
        print('\033[91mERROR: Problem with running fas.doAnno. You can check it with this command:\n%s\033[0m' % annoCmd)
