# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2020 Vinh Tran
#
#  This script is used to prepare data for fdog.
#  It will create a folder within genome_dir with the naming scheme of
#  fdog ([Species acronym]@[NCBI ID]@[Proteome version], e.g
#  HUMAN@9606@3) and a annotation file in JSON format in weight_dir
#  (optional).
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
from pathlib import Path
from Bio import SeqIO
import subprocess
import multiprocessing as mp
from ete3 import NCBITaxa
import re
from datetime import datetime

def checkFileExist(file):
    if not os.path.exists(os.path.abspath(file)):
        sys.exit('%s not found' % file)

def checkOptConflict(replace, delete):
    if delete:
        if replace:
            sys.exit('*** ERROR: only one option can be choose between "--replace" and "--delete"')
    if replace:
        if delete:
            sys.exit('*** ERROR: only one option can be choose between "--replace" and "--delete"')

def checkTaxId(taxId):
    ncbi = NCBITaxa()
    tmp = ncbi.get_rank([taxId])
    try:
        tmp = ncbi.get_rank([taxId])
        rank = tmp[int(taxId)]
        if not rank == 'species':
            print('\033[92mWARNING: rank of %s is not SPECIES (%s)\033[0m' % (taxId, rank))
        else:
            print('\033[92mNCBI taxon info: %s %s\033[0m' % (taxId, ncbi.get_taxid_translator([taxId])[int(taxId)]))
    except:
        print('\033[92mWARNING: %s not found in NCBI taxonomy database!\033[0m' % taxId)

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

def runBlast(args):
    (specName, specFile, outPath) = args
    blastCmd = 'makeblastdb -dbtype prot -in %s -out %s/blast_dir/%s/%s' % (specFile, outPath, specName, specName)
    try:
        subprocess.call([blastCmd], shell = True)
    except:
        sys.exit('Problem with running %s' % blastCmd)
    fileInGenome = "%s/genome_dir/%s/%s.fa" % (outPath, specName, specName)
    fileInBlast = "%s/blast_dir/%s/%s.fa" % (outPath, specName, specName)
    if not Path(fileInBlast).exists():
        os.symlink(fileInGenome, fileInBlast)

def main():
    version = '0.0.2'
    parser = argparse.ArgumentParser(description='You are running fdog.addTaxon version ' + str(version) + '.')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-f', '--fasta', help='FASTA file of input taxon', action='store', default='', required=True)
    required.add_argument('-i', '--taxid', help='Taxonomy ID of input taxon', action='store', default='', required=True, type=int)
    optional.add_argument('-o', '--outPath', help='Path to output directory', action='store', default='')
    optional.add_argument('-n', '--name', help='Acronym name of input taxon', action='store', default='', type=str)
    optional.add_argument('-v', '--verProt', help='Proteome version', action='store', default=1, type=str)
    optional.add_argument('-c', '--coreTaxa', help='Include this taxon to core taxa (i.e. taxa in blast_dir folder)', action='store_true', default=False)
    optional.add_argument('-a', '--noAnno', help='Do NOT annotate this taxon using annoFAS', action='store_true', default=False)
    optional.add_argument('--oldFAS', help='Use old verion of FAS (annoFAS ≤ 1.2.0)', action='store_true', default=False)
    optional.add_argument('--cpus', help='Number of CPUs used for annotation. Default = available cores - 1', action='store', default=0, type=int)
    optional.add_argument('--replace', help='Replace special characters in sequences by "X"', action='store_true', default=False)
    optional.add_argument('--delete', help='Delete special characters in sequences', action='store_true', default=False)
    optional.add_argument('--force', help='Force overwrite existing data', action='store_true', default=False)

    args = parser.parse_args()

    checkFileExist(args.fasta)
    faIn = args.fasta
    name = args.name.upper()
    taxId = str(args.taxid)
    # outPath = str(Path(args.outPath).resolve())
    outPath = args.outPath #str(Path(args.outPath).resolve())
    if outPath == '':
        fdogPath = os.path.realpath(__file__).replace('/addTaxon.py','')
        pathconfigFile = fdogPath + '/bin/pathconfig.txt'
        if not os.path.exists(pathconfigFile):
            sys.exit('No pathconfig.txt found. Please run fdog.setup (https://github.com/BIONF/fDOG/wiki/Installation#setup-fdog).')
        with open(pathconfigFile) as f:
            outPath = f.readline().strip()
    noAnno = args.noAnno
    coreTaxa = args.coreTaxa
    ver = str(args.verProt)
    oldFAS = args.oldFAS
    cpus = args.cpus
    if cpus == 0:
        cpus = mp.cpu_count()-2
    replace = args.replace
    delete = args.delete
    checkOptConflict(replace, delete)
    force = args.force

    ### species name after fdog naming scheme
    checkTaxId(taxId)
    if name == "":
        name = getTaxName(taxId)
    specName = name+'@'+taxId+'@'+ver
    print('Species name\t%s' % specName)

    ### create file in genome_dir
    print('Parsing FASTA file...')
    Path(outPath + '/genome_dir').mkdir(parents = True, exist_ok = True)
    genomePath = outPath + '/genome_dir/' + specName
    Path(genomePath).mkdir(parents = True, exist_ok = True)
    # load fasta seq
    inSeq = SeqIO.to_dict((SeqIO.parse(open(faIn), 'fasta')))
    specFile = genomePath + '/' + specName + '.fa'
    if (not os.path.exists(os.path.abspath(specFile))) or (os.stat(specFile).st_size == 0) or force:
        f = open(specFile, 'w')
        index = 0
        modIdIndex = 0
        longId = 'no'
        tmpDict = {}
        for id in inSeq:
            seq = str(inSeq[id].seq)
            # check ID
            id = re.sub('\|', '_', id)
            if len(id) > 80:
                # modIdIndex = modIdIndex + 1
                # id = specName + "_" + str(modIdIndex)
                longId = 'yes'
            if not id in tmpDict:
                tmpDict[id] = 1
            else:
                index = index + 1
                id = str(id) + '_' + str(index)
                tmpDict[id] = 1
            # check seq
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
        # write .checked file
        cf = open(specFile+'.checked', 'w')
        cf.write(str(datetime.now()))
        cf.close()
        # warning about long header
        if longId == 'yes':
            print('\033[91mWARNING: Headers are longer than 80 characters. It could cause some troubles!\033[0m')
    else:
        print(genomePath + '/' + specName + '.fa already exists!')

    ### create blast db
    if coreTaxa:
        print('Creating Blast DB...')
        Path(outPath + '/blast_dir').mkdir(parents = True, exist_ok = True)
        if (not os.path.exists(os.path.abspath(outPath + '/blast_dir/' + specName + '/' + specName + '.phr'))) or force:
            try:
                runBlast([specName, specFile, outPath])
            except:
                print('\033[91mProblem with creating BlastDB.\033[0m')
        else:
            print('Blast DB already exists!')

    ### create annotation
    if not noAnno:
        Path(outPath + '/weight_dir').mkdir(parents = True, exist_ok = True)
        annoCmd = 'annoFAS -i %s/%s.fa -o %s --cpus %s' % (genomePath, specName, outPath+'/weight_dir', cpus)
        if force:
            annoCmd = annoCmd + " --force"
        if oldFAS:
            print("running old version of FAS...")
            annoCmd = 'annoFAS -i %s/%s.fa -o %s -n %s --cores %s' % (genomePath, specName, outPath+'/weight_dir', specName, cpus)
        try:
            subprocess.call([annoCmd], shell = True)
        except:
            print('\033[91mProblem with running annoFAS. You can check it with this command:\n%s\033[0m' % annoCmd)

    print('Output for %s can be found in %s within genome_dir [and blast_dir, weight_dir] folder[s]' % (specName, outPath))

if __name__ == '__main__':
    main()
