# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2020 Vinh Tran
#
#  This script is used to prepare data for fdog.
#  For each given genome FASTA file, It will create a folder within genome_dir
#  with the naming scheme of fdog ([Species acronym]@[NCBI ID]@[Proteome version]
#  e.g HUMAN@9606@3), a annotation file in JSON format in weight_dir and
#  a blast DB in blast_dir folder (optional).
#  For a long header of original FASTA sequence, only the first word
#  will be taken as the ID of new fasta file, everything after the
#  first whitespace will be removed. If this first word is not unique,
#  an automatically increasing index will be added.
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
import multiprocessing as mp
from ete3 import NCBITaxa
import csv
from io import StringIO
import re
import shutil
from tqdm import tqdm

def checkFileExist(file):
    if not os.path.exists(os.path.abspath(file)):
        sys.exit('%s not found' % file)

def getTaxName(taxId):
    ncbi = NCBITaxa()
    try:
        ncbiName = ncbi.get_taxid_translator([taxId])[int(taxId)]
        ncbiName = re.sub('[^a-zA-Z1-9\s]+', '', ncbiName)
        taxName = ncbiName.split()
        name = taxName[0][:3].upper()+taxName[1][:2].upper()
    except:
        name = "UNK" + taxId
    return(name)

def parseMapFile(mappingFile):
    nameDict = {}
    with open(mappingFile) as f:
      for line in f:
        if not '#' in line:
            tmp = line.split('\t')
            fileName = tmp[0]
            taxId = tmp[1].strip()
            try:
                taxName = tmp[2].strip()
            except:
                taxName = getTaxName(taxId)
            try:
                ver = tmp[3].strip()
            except:
                ver = 1
            # print(taxName+"@"+str(taxId)+"@"+str(ver))
            nameDict[fileName] = (taxName, str(taxId), str(ver))
    return(nameDict)

def runAddTaxon(args):
    (f,n,i,o,c,v,a,cpus,replace,delete,oldFAS) = args
    cmd = 'fdog.addTaxon -f %s -n %s -i %s -o %s -v %s --cpus %s' % (f,n,i,o,v,cpus)
    if c == True:
        cmd = cmd + ' -c'
    if a == True:
        cmd = cmd + ' -a'
    if oldFAS == True:
        cmd = cmd + ' --oldFAS'
    if replace == True:
        cmd = cmd + ' --replace'
    if delete == True:
        cmd = cmd + ' --delete'
    # print(cmd)
    logFile = o + '/addTaxa2fDog.log'
    cmd = cmd + ' >> ' + logFile
    try:
        subprocess.call([cmd], shell = True)
    except:
        sys.exit('Problem running\n%s' % (cmd))

def main():
    version = '0.0.1'
    parser = argparse.ArgumentParser(description='You are running fdog.addTaxa version ' + str(version) + '.')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-i', '--input', help='Path to input folder', action='store', default='', required=True)
    required.add_argument('-m', '--mapping',
                            help='Tab-delimited text file containing <fasta_filename>tab<taxonID>tab<taxonName>tab<genome_version>. The last 2 columns are optional.',
                            action='store', default='', required=True)
    optional.add_argument('-o', '--outPath', help='Path to output directory', action='store', default='')
    optional.add_argument('-c', '--coreTaxa', help='Include these taxa to core taxa (i.e. taxa in blast_dir folder)', action='store_true', default=False)
    optional.add_argument('-a', '--noAnno', help='Do NOT annotate these taxa using annoFAS', action='store_true', default=False)
    optional.add_argument('--oldFAS', help='Use old verion of FAS (annoFAS ≤ 1.2.0)', action='store_true', default=False)
    optional.add_argument('--cpus', help='Number of CPUs used for annotation. Default = available cores - 1', action='store', default=0, type=int)
    optional.add_argument('--replace', help='Replace special characters in sequences by "X"', action='store_true', default=False)
    optional.add_argument('--delete', help='Delete special characters in sequences', action='store_true', default=False)
    optional.add_argument('-f', '--force', help='Force overwrite existing data', action='store_true', default=False)

    ### get arguments
    args = parser.parse_args()
    folIn = args.input
    mapping = args.mapping
    checkFileExist(mapping)
    outPath = args.outPath
    if outPath == '':
        fdogPath = os.path.realpath(__file__).replace('/addTaxa.py','')
        pathconfigFile = fdogPath + '/bin/pathconfig.txt'
        if not os.path.exists(pathconfigFile):
            sys.exit('No pathconfig.txt found. Please run fdog.setup (https://github.com/BIONF/fDOG/wiki/Installation#setup-fdog).')
        with open(pathconfigFile) as f:
            outPath = f.readline().strip()
    noAnno = args.noAnno
    coreTaxa = args.coreTaxa
    oldFAS = args.oldFAS
    cpus = args.cpus
    if cpus == 0:
        cpus = mp.cpu_count()-2
    replace = args.replace
    delete = args.delete
    force = args.force


    ### get existing genomes
    Path(outPath + "/genome_dir").mkdir(parents = True, exist_ok = True)
    Path(outPath + "/weight_dir").mkdir(parents = True, exist_ok = True)
    genomeFiles = listdir(outPath + "/genome_dir")

    ### generate taxon names from mapping file
    nameDict = parseMapFile(mapping)

    ### read all input fasta files and create addTaxon jobs
    jobs = []
    dupList = {}
    faFiles = [f for f in listdir(folIn) if isfile(join(folIn, f))]
    for f in faFiles:
        # tmp = f.split('.')
        if f in nameDict:
            # check duplicated taxon name in existing data
            taxName = '@'.join(nameDict[f])
            flag = 1
            if taxName in genomeFiles:
                if force:
                    shutil.rmtree(outPath + "/genome_dir/" + taxName)
                    if not noAnno:
                        shutil.rmtree(outPath + "/weight_dir/" + taxName)
                else:
                    flag = 0
                dupList[f] = taxName

            if flag == 1:
                fasta = folIn + '/' + f
                name = nameDict[f][0]
                taxid = nameDict[f][1]
                verProt = nameDict[f][2]
                jobs.append([
                    folIn + '/' + f, nameDict[f][0], nameDict[f][1],
                    outPath, coreTaxa, nameDict[f][2], noAnno, cpus, replace, delete, oldFAS
                ])

    if len(dupList) > 0:
        print("These taxa are probably already present in %s:" % (outPath + "/genome_dir"))
        for f in dupList:
            print('\t'+f+'\t'+dupList[f])
        if force:
            print('They will be deleted and re-compiled!')
        else:
            sys.exit("Please remove them from the mapping file or use different Name/ID/Version!")

    print('Parsing...')
    for job in tqdm(jobs):
        # print('@'.join([job[1],job[2],job[5]]) + '\t' + job[0])
        runAddTaxon(job)

    print('Output can be found in %s' % outPath)

if __name__ == '__main__':
    main()
