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
import subprocess
import math
import re
from Bio import SeqIO
from io import StringIO
import random

import fdog.libs.fasta as fasta_fn
import fdog.libs.output as output_fn

##### FUNCTIONS RELATED TO SEQ ALIGNMENT #####

def check_fasta36_executable(fdogPath):
    """ Check if FASTA36 installed in fdogPath """
    try:
        fasta36_cmd = '%s/bin/aligner/bin/ggsearch36' % fdogPath
        subprocess.check_output(fasta36_cmd, shell = True, stderr = subprocess.STDOUT)
        return('%s/bin/aligner/bin/' % fdogPath)
    except:
        try:
            which_fasta36 = subprocess.run(
                'which fasta36', shell = True, capture_output = True, check = True)
            return(which_fasta36.stdout.decode().strip().replace('fasta36',''))
        except subprocess.CalledProcessError as e:
            sys.exit('\033[91mERROR: FASTA36 not found!\033[0m')


def get_muscle_version(aligner):
    """ Check muscle version (3.8 or 5.1)
    Return v3 for v3.8, otherwise v5
    """
    cmd = 'muscle -version'
    try:
        out = subprocess.run(cmd, shell = True, capture_output = True, check = True)
        if 'v3.8' in out.stdout.decode():
            return('v3')
        else:
            return('v5')
    except subprocess.CalledProcessError as e:
        raise RuntimeError("Command '{}' return with error (code {}): {}".format(e.cmd, e.returncode, e.output))
        sys.exit('\033[91mERROR: Error running command\n%s\033[0m' % cmd)


def do_align(aligner, fa_file):
    """ Do alignment using MUSCLE or MAFFT for a multiple fasta file
    Return a dictionary (SeqIO object) containing seq IDs and aligned sequences
    Note: if any input seq is longer than 12.000 aa/nt, only MAFFT can be used
    """
    input_fa = SeqIO.to_dict((SeqIO.parse(open(fa_file), 'fasta')))
    if len(input_fa) == 1:
        return(input_fa)
    # parse output file name (otherwise cause troubles for muscle_v5)
    out_file = fa_file.split('/')[-1].replace('@', '_')
    # check muscle version
    if aligner == 'muscle':
        if get_muscle_version(aligner) == 'v3':
            if fasta_fn.check_long_seq(fa_file, 12000) == 1:
                aligner = 'mafft-linsi'
            else:
                aligner = 'muscle_v3'
        else:
            if fasta_fn.check_long_seq(fa_file, 15000) == 1:
                aligner = 'mafft-linsi'
            else:
                aligner = 'muscle_v5'
    # create alignment command and run
    align_cline = ''
    if aligner == 'muscle_v3':
        align_cline = 'muscle -in %s' % fa_file
    elif aligner == 'muscle_v5':
        align_cline = 'muscle -align %s -output %s.muscle.out' % (fa_file, out_file)
    else:
        align_cline = 'mafft --localpair --maxiterate 1000 %s' % fa_file
    try:
        aln_out = subprocess.run([align_cline], shell = True, capture_output = True, check = True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError("Command '{}' return with error (code {}): {}".format(e.cmd, e.returncode, e.output))
        sys.exit(
            'ERROR: Error doing alignment with %s for %s' % (aligner, fa_file))

    if aligner == 'muscle_v5':
        aln_seq = SeqIO.to_dict((SeqIO.parse(open('%s.muscle.out' % out_file), 'fasta')))
        os.remove('%s.muscle.out' % out_file)
    else:
        aln_io = StringIO(aln_out.stdout.decode().strip())
        aln_seq = SeqIO.to_dict((SeqIO.parse(aln_io,'fasta')))
    return(aln_seq)


def calc_Kimura_dist(aln_dict, id_1, id_2, debug):
    """ Calculate Kimura distance for a pair of sequences
    Input is a dictionary of MSA (see do_align function).
    The Kimura distance is calculated based on perl module
    https://metacpan.org/pod/Bio::Align::ProteinStatistics#D-distance-methods
    """
    matches = 0
    total = 0
    if id_1 in aln_dict and id_2 in aln_dict:
        for a, b in zip(aln_dict[id_1].seq, aln_dict[id_2].seq):
            if a != '-' and b != '-':
                if a == b:
                    matches +=1
                total += 1
        if not total == 0:
            D = 1 - (matches/total)
        else:
            D = 1
        output_fn.print_debug(
            debug, 'Kimura distance',
            'kimura = round(- (math.log( 1 - %s - (0.2 * (%s ** 2)))), 5)' % (D, D))
        try:
            kimura = round(- (math.log( 1 - D - (0.2 * (D ** 2)))), 5)
        except:
            kimura = 999
        return(kimura)
    else:
        sys.exit('ERROR: %s or %s not found in %s!' % (id_1, id_2, aln_dict))


def calc_aln_score(fa1, fa2, aln_strategy = 'local', debugCore = False):
    """ Calculate alignment score for genes in fa2 vs other genes in fa1
    Return dictionary {gene_id:aln_score}
    """
    fdog_path = os.path.realpath(__file__).replace('/libs/alignment.py','')
    fa1_filename = fa1.split("/")[-1]
    fa2_filename = fa2.split("/")[-1]
    if os.path.exists(fa1_filename):
        if fa2_filename == fa1_filename:
            fa2_filename = f'{fa2_filename}.core'
        fa1_filename = f'{fa1_filename}.core'
    os.symlink(fa1, fa1_filename)
    if not fa2_filename == fa1_filename:
        os.symlink(fa2, fa2_filename)
    fasta36_options = f'{fa1_filename} {fa2_filename} -s BP62 -m 9 -d 0 -z -1 -E 100'
    fdog_path = os.path.realpath(__file__).replace('/libs/alignment.py','')
    fasta36_bin = check_fasta36_executable(fdog_path)
    if aln_strategy == 'global':
        fasta36_cmd = '%s/ggsearch36 %s' \
                                % (fasta36_bin, fasta36_options)
    elif aln_strategy == 'glocal':
        fasta36_cmd = '%s/glsearch36 %s' \
                                % (fasta36_bin, fasta36_options)
    else:
        fasta36_cmd = '%s/ssearch36 %s' \
                                % (fasta36_bin, fasta36_options)
    output_fn.print_debug(
        debugCore, 'ALN SCORE',
        'Calculate aln score using FASTA36: %s' % fasta36_cmd)
    try:
        fasta36_out = subprocess.run(
                [fasta36_cmd], shell = True, capture_output = True, check = True)
    except:
        sys.exit('ERROR: Error running FASTA36\n%s' % fasta36_cmd)
    # returns score for genes in fa2
    aln_score = {}
    cand_dict = SeqIO.to_dict((SeqIO.parse(open(fa2), 'fasta')))
    for id in list(cand_dict.keys()):
        aln_score[id[0:60]] = 0
    results = fasta36_out.stdout.decode().split('\n')
    for l in results:
        if len(l) > 1:
            gene_id = l.split()[0]
            if gene_id in aln_score:
                if re.search('\(\s+\d+\)', l):
                    l = re.sub(r'\(\s+','(', l)
                aln_score[gene_id] = aln_score[gene_id] + int(l.split()[2])
    os.remove(fa1_filename)
    if not fa2_filename == fa1_filename:
        os.remove(fa2_filename)
    return(aln_score)
