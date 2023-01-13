# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2022 Vinh Tran
#
#  This script is used to check fdog data which are present in
#  genome_dir, blast_dir and weight_dir
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
from os import listdir
from os.path import isfile, join
from pathlib import Path
import subprocess
import shutil
from Bio import SeqIO
from ete3 import NCBITaxa
import re
import textwrap
from datetime import datetime
import csv
from pkg_resources import get_distribution
from Bio.Blast.Applications import NcbiblastpCommandline


import fdog.libs.zzz as general_fn
import fdog.libs.blast as blast_fn
import fdog.libs.fasta as fasta_fn


def checkOptConflict(concat, replace, delete):
    if concat:
        if not (delete or replace):
            sys.exit('*** ERROR: for rewrite sequences, you need to set either "--delete" or "--replace"!')
    if delete:
        if replace:
            sys.exit('*** ERROR: only one option can be choose between "--replace" and "--delete"')
    if replace:
        if delete:
            sys.exit('*** ERROR: only one option can be choose between "--replace" and "--delete"')


def check_valid_fasta(file):
    """ Check if valid fasta file
    Check if:
    (1) Input file is a fasta file
    (2) If headers are longer than 30 characters
    (3) If headers and sequences contain any space/tab
    (4) If sequences are written in a single line
    """
    spaceChr = (' ', '\t')
    with open(file, 'r') as f:
        f_bkp = f
        # check if input file is a FASTA file
        fasta = SeqIO.parse(f, 'fasta')
        if not any(fasta):
            return({'notFasta': 1})
        else:
            # check for long header
            inSeq = SeqIO.to_dict((SeqIO.parse(open(file), 'fasta')))
            for id in inSeq:
                if len(id) > 30:
                    return({'longHeader': id})
        # check space or tab
        if any(s in f.read() for s in spaceChr):
            return('space')
    # check single line
    nHeader = general_fn.count_line(file, '>', True)
    nSeq = general_fn.count_line(file, '>', False)
    if not nHeader == nSeq:
        return({'multiLine': 1})
    return({'ok': 1})


def check_valid_folder_name(folder):
    """ Check if folder name contains any special characters """
    invalidChr = (' ','|','\t','\'','"','`','´','^','!','$','%','&')
    if any(e in folder for e in invalidChr):
        sys.exit('*** ERROR: Invalid character found in %s' % folder)


def check_valid_seqs(fa_file):
    """ Check if any sequence contains space/tab or special characters """
    spaceChr = (' ', '\t')
    faSeq = SeqIO.parse(open(fa_file),'fasta')
    for fa in faSeq:
        id, seq = fa.description, str(fa.seq)
        c = ''
        if any(e in id for e in spaceChr):
            sys.exit('*** ERROR: Invalid character found in \">%s\" in %s' % (id, fa_file))
        if any(c for c in seq if not c.isalpha()):
            print('*** ERROR: Invalid character "%s" found in the sequence of gene \"%s\" in %s' % (c, id, fa_file))
            sys.exit('You can use "--replace" or "--delete" to solve this issue!')


def rewrite_seqs(fa_file, replace, delete):
    """ Rewrite fasta sequence by replacing or deleting special characters """
    spaceChr = (' ', '\t')
    faSeq = SeqIO.parse(open(fa_file),'fasta')
    with open(fa_file + '.mod', 'w') as tmpOut:
        for fa in faSeq:
            id, seq = fa.description, str(fa.seq)
            if replace:
                seq = re.sub('[^a-zA-Z]', 'X', seq)
            if delete:
                seq = re.sub('[^a-zA-Z]', '', seq)
            tmpOut.write('>%s\n%s\n' % (id, seq))
    os.replace(fa_file + '.mod', fa_file)


def write_faChecked(fa_file):
    """ Add fa.checked file in genome_dir """
    with open(fa_file+'.checked', 'w') as f:
        f.write(str(datetime.now()))


def check_fasta(checkDir, replace, delete, concat):
    """ Check fasta file in genome_dir and blast_dir """
    taxon_list = []
    for fd in general_fn.read_dir(checkDir):
        taxon = fd
        check_valid_folder_name('%s/%s' % (checkDir, taxon))
        for file in listdir('%s/%s' % (checkDir, taxon)):
            if file.endswith('.fa'):
                fa_file = '%s/%s/%s' % (checkDir, taxon, file)
                if os.path.islink(fa_file):
                    fa_file = os.path.realpath(fa_file)
                general_fn.check_file_exist(fa_file)
                checkfa_file = check_valid_fasta(fa_file)
                taxon_list.append(taxon)
                if not os.path.exists('%s.checked' % fa_file):
                    print(taxon)
                    if list(checkfa_file.keys())[0] == 'notFasta':
                        sys.exit('*** ERROR: %s does not look like a fasta file!' % fa_file)
                    elif list(checkfa_file.keys())[0] == 'longHeader':
                        sys.exit('*** ERROR: %s contains long headers! E.g. %s' % (fa_file, list(checkfa_file.values())[0]))
                    elif list(checkfa_file.keys())[0] == 'space':
                        sys.exit('*** ERROR: %s contains spaces/tabs!' % fa_file)
                    elif list(checkfa_file.keys())[0] == 'multiLine':
                        if not concat:
                            print('*** ERROR: %s contains multiple-line sequences!' % fa_file)
                            sys.exit('Please use "--concat" with "--replace" or "--delete" to join them into single lines')
                        else:
                            rewrite_seqs(fa_file, replace, delete)
                    elif list(checkfa_file.keys())[0] == 'ok':
                        if not (delete or replace):
                            check_valid_seqs(fa_file)
                        else:
                            rewrite_seqs(fa_file, replace, delete)
                    write_faChecked(fa_file)
                if not os.path.exists('%s.fai' % fa_file):
                    fasta_fn.read_fasta(fa_file)
    return(taxon_list)


def check_blastdb(blastDir, fdogPath):
    """ Check for outdated blastdb """
    query = '%s/data/infile.fa' % fdogPath
    for fd in general_fn.read_dir(blastDir):
        blast_db = '%s/%s/%s' % (blastDir, fd, fd)
        try:
            blastp_cline = NcbiblastpCommandline(query = query, db = blast_db)
            stdout, stderr = blastp_cline()
        except:
            return([query, blast_db])
        if not os.path.exists('%s/%s/%s.fa.fai' % (blastDir, fd, fd)):
            fai_in_genome = "../../genome_dir/%s/%s.fa.fai" % (fd, fd)
            fai_in_blast = "%s/%s/%s.fa.fai" % (blastDir, fd, fd)
            os.symlink(fai_in_genome, fai_in_blast)
    return([1])


def make_blastdb(blastDir):
    """ Redo (or update) blastdb """
    outPath = blastDir.replace('/blast_dir', '')
    outPath = '/'.join(blastDir.split('/')[0:-1])
    missing = []
    for fd in general_fn.read_dir(blastDir):
        taxon = fd
        ### get genome fasta file
        fa_file = '%s/%s/%s.fa' % (blastDir, taxon, taxon)
        if os.path.islink(fa_file):
            fa_file = os.path.realpath(fa_file)
        if os.path.exists(fa_file):
            ### remove old files
            blast_path = '%s/%s' % (blastDir, taxon)
            shutil.rmtree(blast_path)
            ### Redo blastdb
            Path(blast_path).mkdir(parents = True, exist_ok = True)
            blast_fn.make_blastdb([taxon, fa_file, outPath, True])
            ### make symlink to fasta files
            fa_in_genome = "../../genome_dir/%s/%s.fa" % (taxon, taxon)
            fai_in_genome = "../../genome_dir/%s/%s.fa.fai" % (taxon, taxon)
            fa_in_blast = "%s/%s.fa" % (blast_path, taxon)
            fai_in_blast = "%s/%s.fa.fai" % (blast_path, taxon)
            if not os.path.exists(fa_in_blast):
                os.symlink(fa_in_genome, fa_in_blast)
            if not os.path.exists(fai_in_blast):
                os.symlink(fai_in_genome, fai_in_blast)
            print(taxon)
        else:
            missing.append(taxon)
    return(missing)


def check_missing_json(weightDir, taxon_list):
    """ Check missing annotation for any taxa in blast_dir and genome_dir """
    all_anno = [f for f in listdir(weightDir) if isfile(join(weightDir, f))]
    taxa_anno = [s + '.json' for s in taxon_list]
    s = set(all_anno)
    missing_anno = [x for x in taxa_anno if x not in s]
    return(missing_anno)


def check_complete_anno(weightDir, genomeDir):
    """ Check if an annotation is complete
    I.e. if it contains annotation for all proteins of a species
    """
    all_anno = [f for f in listdir(weightDir) if isfile(join(weightDir, f))]
    for f in all_anno:
        tax = f.replace('.json', '')
        # print('...check annotations for %s' % tax)
        jf = '%s/%s.json' % (weightDir, tax)
        gf = '%s/%s/%s.fa' % (genomeDir, tax, tax)
        cmd = 'fas.checkAnno -s %s -a %s -o %s' % (gf, jf, weightDir)
        try:
            subprocess.call([cmd], shell = True)
        except subprocess.CalledProcessError as e:
            print('*** ERROR: Problem while checking annotation file using fas.checkAnno!')
            print(e.output.decode(sys.stdout.encoding))
            sys.exit()


def check_missing_ncbiID(taxon_list):
    """ Check all taxa in genome_dir and blast_dir
    if they are have valid NCBI taxonomy IDs
    """
    ncbi = NCBITaxa()
    missing_taxa = {}
    present_taxa = {}
    dup_taxa = []
    for t in taxon_list:
        tax_id = t.split('@')[1]
        try:
            taxid2name = ncbi.get_taxid_translator([tax_id])
            if len(taxid2name) < 1:
                if not t+'\t'+str(tax_id) in missing_taxa:
                    missing_taxa[t+'\t'+str(tax_id)] = 1
        except:
            if not t+'\t'+str(tax_id) in missing_taxa:
                missing_taxa[t+'\t'+str(tax_id)] = 1
        if not tax_id in present_taxa:
            present_taxa[tax_id] = t
        else:
            dup_taxa.append('%s\t%s' % (t, present_taxa[tax_id]))
    return(missing_taxa.keys(), dup_taxa)


def main():
    version = get_distribution('fdog').version
    parser = argparse.ArgumentParser(description='You are running fDOG version ' + str(version) + '.')
    parser.add_argument('-g', '--genomeDir', help='Path to search taxa directory (e.g. fdog_dataPath/genome_dir)', action='store', default='')
    parser.add_argument('-b', '--blastDir', help='Path to blastDB directory (e.g. fdog_dataPath/blast_dir)', action='store', default='')
    parser.add_argument('-w', '--weightDir', help='Path to feature annotation directory (e.g. fdog_dataPath/weight_dir)', action='store', default='')
    parser.add_argument('--replace', help='Replace special characters in sequences by "X"', action='store_true', default=False)
    parser.add_argument('--delete', help='Delete special characters in sequences', action='store_true', default=False)
    parser.add_argument('--concat', help='Concatenate multiple-line sequences into single-line', action='store_true', default=False)
    parser.add_argument('--reblast', help='Re-create blast databases', action='store_true', default=False)

    ### get arguments
    args = parser.parse_args()

    genomeDir = args.genomeDir
    blastDir = args.blastDir
    weightDir = args.weightDir
    replace = args.replace
    delete = args.delete
    concat = args.concat
    reblast = args.reblast

    checkOptConflict(concat, replace, delete)
    caution = 0

    ### get fdog dir and assign genomeDir, blastDir, weightDir if not given
    fdogPath = os.path.realpath(__file__).replace('/checkData.py','')
    if not genomeDir or not blastDir or not weightDir:
        pathconfigFile = fdogPath + '/bin/pathconfig.txt'
        if not os.path.exists(pathconfigFile):
            sys.exit('No pathconfig.txt found. Please run fdog.setup (https://github.com/BIONF/fDOG/wiki/Installation#setup-fdog).')
        with open(pathconfigFile) as f:
            dataPath = f.readline().strip()
        if not genomeDir:
            genomeDir = dataPath + "/genome_dir"
        if not blastDir:
            blastDir = dataPath + "/blast_dir"
        if not weightDir:
            weightDir = dataPath + "/weight_dir"

    genomeDir = os.path.abspath(genomeDir)
    blastDir = os.path.abspath(blastDir)
    weightDir = os.path.abspath(weightDir)

    ### check genomeDir and blastDir
    print('=> Checking %s...' % genomeDir)
    genomeTaxa = check_fasta(genomeDir, replace, delete, concat)
    print('=> Checking %s...' % blastDir)
    blastTaxa = check_fasta(blastDir, replace, delete, concat)
    check_blast = check_blastdb(blastDir, fdogPath)
    if not check_blast[0] == 1:
        if not reblast:
            print('*** ERROR: Version incompatible between BlastDB and BLAST program!')
            print('For checking, run: blastp -query %s -db %s' % (check_blast[0], check_blast[1]))
            print('Consider using --reblast option to update old BlastDBs!')
            sys.exit()
        else:
            failed_blast = make_blastdb(blastDir)
            if len(failed_blast) > 0:
                print('*** WARNING: Some BlastDBs cannot be created:\n%s' % ', '.join(failed_blast))
            else:
                print('All old BlastDBs have been updated!')

    ### check weightDir
    print('=> Checking %s...' % weightDir)
    missing_anno = check_missing_json(weightDir, general_fn.join_2lists(genomeTaxa, blastTaxa))
    if len(missing_anno) > 0:
        print('\033[92m*** WARNING: Annotation files not found for:\033[0m')
        print(*missing_anno, sep = "\n")
        print('NOTE: You still can run fdog without FAS using the option "-fasoff"')
        caution = 1
    check_complete_anno(weightDir, genomeDir)

    ### check ncbi IDs
    print('=> Checking NCBI taxonomy IDs...')
    missing_taxa, dup_taxa = check_missing_ncbiID(general_fn.join_2lists(genomeTaxa, blastTaxa))
    if (len(missing_taxa) > 0):
        print('\033[92m*** WARNING: Taxa not found in current local NCBI taxonomy database:\033[0m')
        print(*missing_taxa, sep = "\n")
        print('==> NOTE: You still can run fDOG with those taxa, but they will not be included in the core set compilation!')
        caution = 1
    if (len(dup_taxa) > 0):
        print('\033[92m*** WARNING: These taxa have the same NCBI taxonomy IDs:\033[0m')
        print(*dup_taxa, sep = "\n")
        print('==> NOTE: This could lead to some conflicts!')
        caution = 1

    print('---------------------------------')
    if caution == 1:
        print('==> Done! Data are ready to use WITH CAUTION!')
    else:
        print('==> Done! Data are ready to use!')

if __name__ == '__main__':
    main()
