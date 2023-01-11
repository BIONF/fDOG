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
import shutil
from datetime import datetime
from pkg_resources import get_distribution

import fdog.libs.zzz as general_fn
import fdog.libs.tree as tree_fn
import fdog.libs.blast as blast_fn
import fdog.libs.fasta as fasta_fn


def check_conflict_opts(replace, delete):
    """ Check if both replace and delete are both specified """
    if delete:
        if replace:
            sys.exit('*** ERROR: only one option can be choose between "--replace" and "--delete"')
    if replace:
        if delete:
            sys.exit('*** ERROR: only one option can be choose between "--replace" and "--delete"')


def check_tax_id(taxId):
    """ Check valid taxon ID """
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


def create_genome(faIn, genome_path, spec_name, force):
    """ Create fa and fai in genome_dir """

    ### load fasta seq
    in_seq = SeqIO.to_dict((SeqIO.parse(open(faIn), 'fasta')))
    if not os.path.exists(genome_path):
            Path(genome_path).mkdir(parents = True, exist_ok = True)
    spec_file = genome_path + '/' + spec_name + '.fa'
    if (not os.path.exists(os.path.abspath(spec_file))) or (os.stat(spec_file).st_size == 0) or force:
        f = open(spec_file, 'w')
        # index = 0
        # modIdIndex = 0
        # longId = 'no'
        # tmpDict = {}
        # with open(spec_file + '.mapping', 'a') as mappingFile:
        for id in in_seq:
            seq = str(in_seq[id].seq)
            # check ID
            # oriId = id
            if ' ' in id:
                sys.exit('\033[91mERROR: Sequence IDs (e.g. %s) must not contain space(s)!\033[0m' % id)
            else:
                if '|' in id:
                    print('\033[91mWARNING: Sequence IDs contain pipe(s). They will be replaced by "_"!\033[0m')
                    id = re.sub('\|', '_', id)
            # if len(id) > 20:
            #     modIdIndex = modIdIndex + 1
            #     id = modIdIndex
            #     longId = 'yes'
            # if not id in tmpDict:
            #     tmpDict[id] = 1
            # else:
            #     index = index + 1
            #     id = str(index)
            #     tmpDict[id] = 1
            # mappingFile.write('%s\t%s\n' % (id, oriId))
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
        ### create index file
        fasta_fn.read_fasta(spec_file)
        ### write .checked file
        cf = open(spec_file+'.checked', 'w')
        cf.write(str(datetime.now()))
        cf.close()
        # ### warning about long header
        # if longId == 'yes':
        #     print('\033[91mWARNING: Some headers longer than 80 characters have been automatically shortened. PLease check the %s.mapping file for details!\033[0m' % spec_file)
    else:
        print(genome_path + '/' + spec_name + '.fa already exists!')
    return(spec_file)


def create_blastdb(outPath, spec_name, spec_file, force):
    blast_path = '%s/blast_dir/%s' % (outPath, spec_name)
    if (not os.path.exists(os.path.abspath('%s/%s.phr' % (blast_path, spec_name)))) or force:
        blast_fn.make_blastdb([spec_name, spec_file, outPath])
        ### make symlink to fasta files
        fa_in_genome = "../../genome_dir/%s/%s.fa" % (spec_name, spec_name)
        fai_in_genome = "../../genome_dir/%s/%s.fa.fai" % (spec_name, spec_name)
        fa_in_blast = "%s/%s.fa" % (blast_path, spec_name)
        fai_in_blast = "%s/%s.fa.fai" % (blast_path, spec_name)
        if not os.path.exists(fa_in_blast):
            os.symlink(fa_in_genome, fa_in_blast)
        if not os.path.exists(fai_in_blast):
            os.symlink(fai_in_genome, fai_in_blast)
    else:
        print('Blast DB already exists!')


def create_annoFile(outPath, spec_name, genome_path, cpus, force):
    annoCmd = 'fas.doAnno -i %s/%s.fa -o %s --cpus %s' % (genome_path, spec_name, outPath+'/weight_dir', cpus)
    if force:
        annoCmd = annoCmd + " --force"
    try:
        subprocess.call([annoCmd], shell = True)
    except:
        print('\033[91mProblem with running fas.doAnno. You can check it with this command:\n%s\033[0m' % annoCmd)


def main():
    version = get_distribution('fdog').version
    parser = argparse.ArgumentParser(description='You are running fDOG version ' + str(version) + '.')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-f', '--fasta', help='FASTA file of input taxon', action='store', default='', required=True)
    required.add_argument('-i', '--taxid', help='Taxonomy ID of input taxon', action='store', default='', required=True, type=int)
    optional.add_argument('-o', '--outPath', help='Path to output directory', action='store', default='')
    optional.add_argument('-n', '--name', help='Acronym name of input taxon', action='store', default='', type=str)
    optional.add_argument('-v', '--verProt', help='Proteome version', action='store', default='', type=str)
    optional.add_argument('-c', '--coreTaxa', help='Include this taxon to core taxa (i.e. taxa in blast_dir folder)', action='store_true', default=False)
    optional.add_argument('-a', '--noAnno', help='Do NOT annotate this taxon using fas.doAnno', action='store_true', default=False)
    optional.add_argument('--cpus', help='Number of CPUs used for annotation. Default = available cores - 1', action='store', default=0, type=int)
    optional.add_argument('--replace', help='Replace special characters in sequences by "X"', action='store_true', default=False)
    optional.add_argument('--delete', help='Delete special characters in sequences', action='store_true', default=False)
    optional.add_argument('--force', help='Force overwrite existing data', action='store_true', default=False)

    args = parser.parse_args()

    general_fn.check_file_exist(args.fasta)
    faIn = args.fasta
    name = args.name.upper()
    taxId = str(args.taxid)
    outPath = args.outPath
    if outPath == '':
        fdogPath = os.path.realpath(__file__).replace('/addTaxon.py','')
        pathconfigFile = fdogPath + '/bin/pathconfig.txt'
        if not os.path.exists(pathconfigFile):
            sys.exit('No pathconfig.txt found. Please run fdog.setup (https://github.com/BIONF/fDOG/wiki/Installation#setup-fdog).')
        with open(pathconfigFile) as f:
            outPath = f.readline().strip()
    outPath = os.path.abspath(outPath)
    noAnno = args.noAnno
    coreTaxa = args.coreTaxa
    ver = str(args.verProt)
    if ver == '':
        ver = datetime.today().strftime('%y%m%d')
    cpus = args.cpus
    if cpus == 0:
        cpus = mp.cpu_count()-2
    replace = args.replace
    delete = args.delete
    check_conflict_opts(replace, delete)
    force = args.force

    ### species name after fdog naming scheme
    check_tax_id(taxId)
    if name == "":
        name = tree_fn.get_tax_name(taxId)
    spec_name = name+'@'+taxId+'@'+ver
    print('Species name\t%s' % spec_name)


    ### initiate paths
    Path(outPath + '/genome_dir').mkdir(parents = True, exist_ok = True)
    genome_path = outPath + '/genome_dir/' + spec_name
    Path(genome_path).mkdir(parents = True, exist_ok = True)

    ### remove old folder if force is set
    if force:
        if os.path.exists(outPath + '/genome_dir/' + spec_name):
            shutil.rmtree(outPath + '/genome_dir/' + spec_name)
        if os.path.exists(outPath + '/blast_dir/' + spec_name):
            shutil.rmtree(outPath + '/blast_dir/' + spec_name)

    ### create file in genome_dir
    print('Parsing FASTA file...')
    spec_file = create_genome(faIn, genome_path, spec_name, force)
    out_msg = 'Output for %s can be found in %s within genome_dir'  % (spec_name, outPath)

    ### create blast db
    if coreTaxa:
        print('Creating Blast DB...')
        Path(outPath + '/blast_dir').mkdir(parents = True, exist_ok = True)
        create_blastdb(outPath, spec_name, spec_file, force)
        out_msg = '%s, blast_dir' % out_msg

    ### create annotation
    if not noAnno:
        Path(outPath + '/weight_dir').mkdir(parents = True, exist_ok = True)
        create_annoFile(outPath, spec_name, genome_path, cpus, force)
        if os.path.exists('%s/weight_dir/tmp' % outPath):
            shutil.rmtree('%s/weight_dir/tmp' % outPath)
        out_msg = '%s, weight_dir' % out_msg

    print('%s folder[s]' % out_msg)

if __name__ == '__main__':
    main()
