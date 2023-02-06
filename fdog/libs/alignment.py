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
from Bio.Align.Applications import MuscleCommandline
from Bio.Align.Applications import MafftCommandline
from io import StringIO

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


def do_align(aligner, fa_file):
    """ Do alignment using MUSCLE or MAFFT for a multiple fasta file
    Return a dictionary (SeqIO object) containing seq IDs and aligned sequences
    Note: if any input seq is longer than 12.000 aa/nt, only MAFFT can be used
    """
    if fasta_fn.check_long_seq(fa_file) == 1:
        aligner = 'mafft-linsi'
    if aligner == 'muscle':
        align_cline = MuscleCommandline(input = fa_file)
    else:
        align_cline = MafftCommandline(
            input = fa_file, localpair = True, maxiterate = 1000)
    try:
        stdout, stderr = align_cline()
        aln_io = StringIO(stdout)
        aln_seq = SeqIO.to_dict((SeqIO.parse(aln_io,'fasta')))
        return(aln_seq)
    except:
        sys.exit(
            'ERROR: Error doing alignment with %s for %s\n%s' % (aligner, fa_file, align_cline))


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
        sys.exit('%s or %s not found in %s!' % (id_1, id_2, aln_dict))


def calc_aln_score(fa1, fa2, aln_strategy = 'local', debugCore = False):
    """ Calculate alignment score for genes in fa2 vs other genes in fa1
    Return dictionary {gene_id:aln_score}
    """
    fdog_path = os.path.realpath(__file__).replace('/libs/alignment.py','')
    fasta36_options = '%s %s -s BP62 -m 9 -d 0 -z -1 -E 100' % (fa1, fa2)
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
        sys.exit('Error running FASTA36\n%s' % fasta36_out)
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
    return(aln_score)
