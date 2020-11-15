# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2020 Vinh Tran
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
from Bio import SeqIO
import re
from datetime import datetime
import csv

def checkFileExist(file):
    if not os.path.exists(os.path.abspath(file)):
        sys.exit('%s not found' % file)

def countLine(file,pattern,contain):
    nline = 0
    with open(file, 'r') as f:
        for line in f:
            if contain:
                if pattern in line:
                    nline = nline + 1
            else:
                if not pattern in line:
                    nline = nline + 1
    return(nline)

def join2Lists(first_list, second_list):
    in_first = set(first_list)
    in_second = set(second_list)
    in_second_but_not_in_first = in_second - in_first
    out = first_list + list(in_second_but_not_in_first)
    return(out)

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

def checkValidFasta(file):
    spaceChr = (' ', '\t')
    with open(file, 'r') as f:
        f_bkp = f
        # check if input file a FASTA file
        fasta = SeqIO.parse(f, 'fasta')
        if not any(fasta):
            return('notFasta')
        # check space or tab
        if any(s in f.read() for s in spaceChr):
            return('space')
    # check single line
    nHeader = countLine(file, '>', True)
    nSeq = countLine(file, '>', False)
    if not nHeader == nSeq:
        return('multiLine')
    return('ok')

def checkValidFolderName(folder):
    invalidChr = (' ','|','\t','\'','"','`','´','^','!','$','%','&')
    if any(e in folder for e in invalidChr):
        sys.exit('*** ERROR: Invalid character found in %s' % folder)

def checkValidSeqs(faFile):
    spaceChr = (' ', '\t')
    faSeq = SeqIO.parse(open(faFile),'fasta')
    for fa in faSeq:
        id, seq = fa.description, str(fa.seq)
        if any(e in id for e in spaceChr):
            sys.exit('*** ERROR: Invalid character found in \">%s\" in %s' % (id, faFile))
        if any(c for c in seq if not c.isalpha()):
            print('*** ERROR: Invalid character "%s" found in the sequence of gene \"%s\" in %s' % (c, id, faFile))
            sys.exit('You can use "--replace" or "--delete" to solve this issue!')

def rewriteSeqs(faFile, replace, delete):
    spaceChr = (' ', '\t')
    faSeq = SeqIO.parse(open(faFile),'fasta')
    with open(faFile + '.mod', 'w') as tmpOut:
        for fa in faSeq:
            id, seq = fa.description, str(fa.seq)
            if replace:
                seq = re.sub('[^a-zA-Z]', 'X', seq)
            if delete:
                seq = re.sub('[^a-zA-Z]', '', seq)
            tmpOut.write('>%s\n%s\n' % (id, seq))
    os.replace(faFile + '.mod', faFile)

def writeCheckedFile(faFile):
    with open(faFile+'.checked', 'w') as f:
        f.write(str(datetime.now()))

def checkDataFolder(checkDir, replace, delete, concat):
    taxaList = []
    for fd in listdir(checkDir):
        if not fd.startswith('.'):
            taxon = fd
            checkValidFolderName(checkDir+'/'+taxon)
            getFaCmd = 'ls %s/%s/%s.fa*' % (checkDir, taxon, taxon)
            try:
                faFiles = subprocess.check_output([getFaCmd], shell=True).decode(sys.stdout.encoding).strip().split('\n')
                for faFile in faFiles:
                    if os.path.islink(faFile):
                        faFile = os.path.realpath(faFile)
                    checkFileExist(faFile)
                    if not '.checked' in faFile:
                        if not os.path.exists(faFile+".checked"):
                            checkFaFile = checkValidFasta(faFile)
                            if checkFaFile == 'notFasta':
                                sys.exit('*** ERROR: %s does not look like a fasta file!' % faFile)
                            elif checkFaFile == 'space':
                                sys.exit('*** ERROR: %s contains spaces/tabs!' % faFile)
                            elif checkFaFile == 'multiLine':
                                if not concat:
                                    print('*** ERROR: %s contains multiple-line sequences!' % faFile)
                                    sys.exit('Please use "--concat" with "--replace" or "--delete" to join them into single lines')
                                else:
                                    rewriteSeqs(faFile, replace, delete)
                            elif checkFaFile == 'ok':
                                if not (delete or replace):
                                    checkValidSeqs(faFile)
                                else:
                                    rewriteSeqs(faFile, replace, delete)
                            writeCheckedFile(faFile)
                            print(fd)
                taxaList.append(fd)
            except subprocess.CalledProcessError as e:
                print('*** ERROR: Problem while searching for fasta file')
                print(e.output.decode(sys.stdout.encoding))
                sys.exit()
    return(taxaList)

def checkCompleteAnno(weightDir, taxaList):
    allAnno = [f for f in listdir(weightDir) if isfile(join(weightDir, f))]
    taxaAnno = [s + '.json' for s in taxaList]
    s = set(allAnno)
    missingAnno = [x for x in taxaAnno if x not in s]
    return(missingAnno)

def checkMissingNcbiID(namesDmp, taxaList):
    ncbiId = {}
    with open(namesDmp, 'r') as f:
        lines = f.readlines()
        for x in lines:
            taxId = x.split('\t')[0]
            if not taxId in ncbiId:
                ncbiId[taxId] = 1
    f.close()
    missingTaxa = {}
    presentTaxa = {}
    dupTaxa = []
    for t in taxaList:
        taxId = t.split('@')[1]
        if not taxId in ncbiId:
            if not t+'\t'+str(taxId) in missingTaxa:
                missingTaxa[t+'\t'+str(taxId)] = 1
        if not taxId in presentTaxa:
            presentTaxa[taxId] = t
        else:
            dupTaxa.append('%s\t%s' % (t, presentTaxa[taxId]))
    return(missingTaxa.keys(), dupTaxa)

def main():
    version = '0.0.2'
    parser = argparse.ArgumentParser(description='You are running fdog.checkData version ' + str(version) + '.')
    parser.add_argument('-g', '--genomeDir', help='Path to search taxa directory (e.g. fdog_dataPath/genome_dir)', action='store', default='')
    parser.add_argument('-b', '--blastDir', help='Path to blastDB directory (e.g. fdog_dataPath/blast_dir)', action='store', default='')
    parser.add_argument('-w', '--weightDir', help='Path to feature annotation directory (e.g. fdog_dataPath/weight_dir)', action='store', default='')
    parser.add_argument('--replace', help='Replace special characters in sequences by "X"', action='store_true', default=False)
    parser.add_argument('--delete', help='Delete special characters in sequences', action='store_true', default=False)
    parser.add_argument('--concat', help='Concatenate multiple-line sequences into single-line', action='store_true', default=False)

    ### get arguments
    args = parser.parse_args()

    genomeDir = args.genomeDir
    blastDir = args.blastDir
    weightDir = args.weightDir
    replace = args.replace
    delete = args.delete
    concat = args.concat

    checkOptConflict(concat, replace, delete)
    caution = 0

    ### get fdog dir and assign genomeDir, blastDir, weightDir if not given
    fdogPath = os.path.realpath(__file__).replace('/checkData.py','')
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

    ### check genomeDir and blastDir
    print('=> Checking %s...' % genomeDir)
    genomeTaxa = checkDataFolder(os.path.abspath(genomeDir), replace, delete, concat)
    print('=> Checking %s...' % blastDir)
    blastTaxa = checkDataFolder(os.path.abspath(blastDir), replace, delete, concat)

    ### check weightDir
    print('=> Checking %s...' % weightDir)
    missingAnno = checkCompleteAnno(weightDir, join2Lists(genomeTaxa, blastTaxa))
    if len(missingAnno) > 0:
        print('\033[92m*** WARNING: Annotations not found for:\033[0m')
        print(*missingAnno, sep = "\n")
        print('NOTE: You still can run fdog without FAS using the option "-fasoff"')
        caution = 1

    ### check ncbi IDs
    print('=> Checking NCBI taxonomy IDs...')
    namesDmp = fdogPath + '/taxonomy/names.dmp'
    checkFileExist(namesDmp)
    missingTaxa, dupTaxa = checkMissingNcbiID(namesDmp, join2Lists(genomeTaxa, blastTaxa))
    if (len(missingTaxa) > 0):
        print('\033[92m*** WARNING: Taxa not found in current fdog\'s NCBI taxonomy database:\033[0m')
        print(*missingTaxa, sep = "\n")
        print('NOTE: You still can run fdog, but they will not be included in the core set compilation!')
        caution = 1
    if (len(dupTaxa) > 0):
        print('\033[92m*** WARNING: These taxa have the same NCBI taxonomy IDs:\033[0m')
        print(*dupTaxa, sep = "\n")
        print('NOTE: This could lead to some conflicts!')
        caution = 1

    print('---------------------------------')
    if caution == 1:
        print('Done! Data are ready to use with caution!')
    else:
        print('Done! Data are ready to use!')

if __name__ == '__main__':
    main()
