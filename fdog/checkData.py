# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2022 Vinh Tran
#
#  This script is used to check fdog data which are present in
#  searchTaxa_dir, coreTaxa_dir and annotation_dir
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
import errno
import argparse
from os import listdir
from os.path import isfile, join
from pathlib import Path
import subprocess
import shutil
from Bio import SeqIO
from ete3 import NCBITaxa
import re
from datetime import datetime
import multiprocessing as mp
from tqdm import tqdm
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
    """ Add fa.checked file in searchTaxa_dir """
    with open(fa_file+'.checked', 'w') as f:
        f.write(str(datetime.now()))


def check_fasta(args):
    """ Check fasta file in searchTaxa_dir and coreTaxa_dir """
    (taxon, file, checkDir, replace, delete, concat) = args
    fa_file = '%s/%s/%s' % (checkDir, taxon, file)
    if os.path.islink(fa_file):
        fa_file = os.path.realpath(fa_file)
    general_fn.check_file_exist(fa_file)
    checkfa_file = check_valid_fasta(fa_file)
    if not os.path.exists('%s.checked' % fa_file):
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
    return(taxon)


def run_check_fasta(checkDir, replace, delete, concat):
    """ Run check_fasta fn """
    jobs = []
    for taxon in general_fn.read_dir(checkDir):
        check_valid_folder_name('%s/%s' % (checkDir, taxon))
        for file in listdir('%s/%s' % (checkDir, taxon)):
            if file.endswith('.fa'):
                jobs.append([taxon, file, checkDir, replace, delete, concat])
    cpus = mp.cpu_count()-1
    pool = mp.Pool(cpus)
    taxon_list = []
    for _ in tqdm(pool.imap_unordered(check_fasta, jobs), total=len(jobs)):
        taxon_list.append(_)
    return(taxon_list)


def check_blastdb(args):
    """ Check for outdated blastdb """
    (query, taxon, coreTaxa_dir, searchTaxa_dir) = args
    blast_db = '%s/%s/%s' % (coreTaxa_dir, taxon, taxon)
    try:
        blastp_cline = NcbiblastpCommandline(query = query, db = blast_db)
        stdout, stderr = blastp_cline()
    except:
        return([query, blast_db])
    fai_in_genome = "%s/%s/%s.fa.fai" % (searchTaxa_dir, taxon, taxon)
    fai_in_blast = "%s/%s/%s.fa.fai" % (coreTaxa_dir, taxon, taxon)
    # check if fai_in_blast is a valid symlink
    if os.path.islink(fai_in_blast):
        if not os.path.exists(os.readlink(fai_in_blast)):
            if os.path.exists(fai_in_genome):
                try:
                    os.remove('%s/%s/%s.fa.fai' % (coreTaxa_dir, taxon, taxon))
                except OSError as e:
                    if e.errno != errno.ENOENT:
                        raise
                os.symlink(fai_in_genome, fai_in_blast)
    # or that file doesn't exist
    else:
        if not os.path.exists(fai_in_blast):
            if os.path.exists(fai_in_genome):
                os.symlink(fai_in_genome, fai_in_blast)


def run_check_blastdb(coreTaxa_dir, searchTaxa_dir, fdogPath):
    """ Run check_blastdb fn """
    query = '%s/data/infile.fa' % fdogPath
    jobs = []
    for fd in general_fn.read_dir(coreTaxa_dir):
        jobs.append([query, fd, coreTaxa_dir, searchTaxa_dir])
    cpus = mp.cpu_count()-1
    pool = mp.Pool(cpus)
    out = []
    for _ in tqdm(pool.imap_unordered(check_blastdb, jobs), total=len(jobs)):
        out.append(_)
    return([1])


def create_blastdb(args):
    """ Redo (or update) blastdb """
    (taxon, coreTaxa_dir, searchTaxa_dir, outPath) = args
    fa_file = '%s/%s/%s.fa' % (coreTaxa_dir, taxon, taxon)
    if os.path.islink(fa_file):
        fa_file = os.path.realpath(fa_file)
    if not os.path.exists(fa_file):
        fa_file = '%s/%s/%s.fa' % (searchTaxa_dir, taxon, taxon)
    if os.path.exists(fa_file):
        ### remove old files
        blast_path = '%s/%s' % (coreTaxa_dir, taxon)
        shutil.rmtree(blast_path)
        ### Redo blastdb
        Path(blast_path).mkdir(parents = True, exist_ok = True)
        blast_fn.make_blastdb([taxon, fa_file, outPath, coreTaxa_dir, searchTaxa_dir, True])
        ### make symlink to fasta files
        fa_in_genome = "%s/%s/%s.fa" % (searchTaxa_dir, taxon, taxon)
        fai_in_genome = "%s/%s/%s.fa.fai" % (searchTaxa_dir, taxon, taxon)
        fa_in_blast = "%s/%s.fa" % (blast_path, taxon)
        fai_in_blast = "%s/%s.fa.fai" % (blast_path, taxon)
        if not os.path.exists(fa_in_blast):
            os.symlink(fa_in_genome, fa_in_blast)
        if not os.path.exists(fai_in_blast):
            os.symlink(fai_in_genome, fai_in_blast)
        return None
    else:
        return(taxon)


def run_create_blastdb(coreTaxa_dir, searchTaxa_dir):
    """ Run create_blastdb fn """
    outPath = '/'.join(coreTaxa_dir.split('/')[0:-1])
    jobs = []
    for fd in general_fn.read_dir(coreTaxa_dir):
        jobs.append([fd, coreTaxa_dir, searchTaxa_dir, outPath])
    cpus = mp.cpu_count()-1
    pool = mp.Pool(cpus)
    out = []
    for _ in tqdm(pool.imap_unordered(create_blastdb, jobs), total=len(jobs)):
        out.append(_)
    return([i for i in out if i is not None])


def check_missing_json(annotation_dir, taxon_list):
    """ Check missing annotation for any taxa in coreTaxa_dir and searchTaxa_dir """
    all_anno = [f for f in listdir(annotation_dir) if isfile(join(annotation_dir, f))]
    taxa_anno = [s + '.json' for s in taxon_list]
    s = set(all_anno)
    missing_anno = [x for x in taxa_anno if x not in s]
    return(missing_anno)


def check_complete_anno(args):
    """ Check if an annotation is complete
    I.e. if it contains annotation for all proteins of a species
    """
    (gf,jf, annotation_dir, updateJson) = args
    cmd = 'fas.checkAnno -s %s -a %s -o %s --noAnno' % (gf, jf, annotation_dir)
    if updateJson:
        cmd = '%s --update' % cmd
    try:
        subprocess.call([cmd], shell = True, stdout=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        print('*** ERROR: Problem while checking annotation file using fas.checkAnno!')
        print(e.output.decode(sys.stdout.encoding))
        sys.exit()


def run_check_complete_anno(annotation_dir, searchTaxa_dir, coreTaxa_dir, updateJson):
    """ Run check_complete_anno fn """
    all_anno = [f for f in listdir(annotation_dir) if isfile(join(annotation_dir, f))]
    jobs = []
    for f in all_anno:
        tax = f.replace('.json', '')
        # print('...check annotations for %s' % tax)
        jf = '%s/%s.json' % (annotation_dir, tax)
        gf = '%s/%s/%s.fa' % (searchTaxa_dir, tax, tax)
        if not os.path.exists(gf):
            gf = '%s/%s/%s.fa' % (coreTaxa_dir, tax, tax)
        jobs.append([gf,jf, annotation_dir, updateJson])
    cpus = mp.cpu_count()-1
    pool = mp.Pool(cpus)
    out = []
    for i in jobs:
        check_complete_anno(i)
    # for _ in tqdm(pool.imap_unordered(check_complete_anno, jobs), total=len(jobs)):
    #     out.append(_)
    return None


def check_missing_ncbiID(taxon_list):
    """ Check all taxa in searchTaxa_dir and coreTaxa_dir
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


def run_check(args):
    (searchTaxa_dir, coreTaxa_dir, annotation_dir, replace, delete, concat, reblast, updateJson) = args
    checkOptConflict(concat, replace, delete)
    caution = 0

    ### get fdog dir and assign searchTaxa_dir, coreTaxa_dir, annotation_dir if not given
    fdogPath = os.path.realpath(__file__).replace('/checkData.py','')
    if not searchTaxa_dir or not coreTaxa_dir or not annotation_dir:
        pathconfigFile = fdogPath + '/bin/pathconfig.yml'
        if not os.path.exists(pathconfigFile):
            sys.exit('No pathconfig.yml found. Please run fdog.setup (https://github.com/BIONF/fDOG/wiki/Installation#setup-fdog).')
        cfg = general_fn.load_config(pathconfigFile)
        try:
            dataPath = cfg['dataPath']
        except:
            dataPath = os.getcwd()

        if not searchTaxa_dir:
            try:
                searchTaxa_dir = cfg['searchpath']
            except:
                searchTaxa_dir = dataPath + '/searchTaxa_dir'
        if not coreTaxa_dir:
            try:
                coreTaxa_dir = cfg['corepath']
            except:
                coreTaxa_dir = dataPath + "/coreTaxa_dir"
        if not annotation_dir:
            try:
                annotation_dir = cfg['annopath']
            except:
                annotation_dir = dataPath + "/annotation_dir"

    searchTaxa_dir = os.path.abspath(searchTaxa_dir)
    coreTaxa_dir = os.path.abspath(coreTaxa_dir)
    annotation_dir = os.path.abspath(annotation_dir)

    ### check searchTaxa_dir
    print('=> Checking %s...' % searchTaxa_dir)
    search_taxa = run_check_fasta(searchTaxa_dir, replace, delete, concat)

    ### check coreTaxa_dir
    if reblast:
        print('=> (Re-)Creating blastDBs...')
        failed_blast = run_create_blastdb(coreTaxa_dir, searchTaxa_dir)
        if len(failed_blast) > 0:
            print('*** WARNING: Some BlastDBs cannot be created:\n%s' % ', '.join(failed_blast))
        else:
            print('All old BlastDBs have been updated!')
    print('=> Checking %s...' % coreTaxa_dir)
    core_taxa = run_check_fasta(coreTaxa_dir, replace, delete, concat)
    check_blast = run_check_blastdb(coreTaxa_dir, searchTaxa_dir, fdogPath)

    if not check_blast[0] == 1:
        print('*** ERROR: Version incompatible between BlastDB and BLAST program!')
        print('For checking, run: blastp -query %s -db %s' % (check_blast[0], check_blast[1]))
        print('Consider using --reblast option to update old BlastDBs!')
        sys.exit()

    ### check annotation_dir
    print('=> Checking %s...' % annotation_dir)
    missing_anno = check_missing_json(annotation_dir, general_fn.join_2lists(search_taxa, core_taxa))
    if len(missing_anno) > 0:
        print('\033[92m*** WARNING: Annotation files not found for:\033[0m')
        print(*missing_anno, sep = "\n")
        print('NOTE: You still can run fdog without FAS using the option "-fasoff"')
        caution = 1
    run_check_complete_anno(annotation_dir, searchTaxa_dir, coreTaxa_dir, updateJson)

    ### check ncbi IDs
    print('=> Checking NCBI taxonomy IDs...')
    missing_taxa, dup_taxa = check_missing_ncbiID(general_fn.join_2lists(search_taxa, core_taxa))
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
    return(caution)

def main():
    version = get_distribution('fdog').version
    parser = argparse.ArgumentParser(description='You are running fDOG version ' + str(version) + '.')
    parser.add_argument('-s', '--searchTaxa_dir', help='Path to search taxa directory (e.g. fdog_dataPath/searchTaxa_dir)', action='store', default='')
    parser.add_argument('-c', '--coreTaxa_dir', help='Path to blastDB directory (e.g. fdog_dataPath/coreTaxa_dir)', action='store', default='')
    parser.add_argument('-a', '--annotation_dir', help='Path to feature annotation directory (e.g. fdog_dataPath/annotation_dir)', action='store', default='')
    parser.add_argument('--replace', help='Replace special characters in sequences by "X"', action='store_true', default=False)
    parser.add_argument('--delete', help='Delete special characters in sequences', action='store_true', default=False)
    parser.add_argument('--concat', help='Concatenate multiple-line sequences into single-line', action='store_true', default=False)
    parser.add_argument('--reblast', help='Re-create blast databases', action='store_true', default=False)
    parser.add_argument('--updateJson', help='Update annotation json file to FAS >=1.16', action='store_true', default=False)

    ### get arguments
    args = parser.parse_args()

    searchTaxa_dir = args.searchTaxa_dir
    coreTaxa_dir = args.coreTaxa_dir
    annotation_dir = args.annotation_dir
    replace = args.replace
    delete = args.delete
    concat = args.concat
    reblast = args.reblast
    updateJson = args.updateJson

    caution = run_check([searchTaxa_dir, coreTaxa_dir, annotation_dir, replace, delete, concat, reblast, updateJson])
    if caution == 1:
        print('==> Done! Data are ready to use WITH CAUTION!')
    else:
        print('==> Done! Data are ready to use!')


if __name__ == '__main__':
    main()
