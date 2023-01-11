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
import pyhmmer


##### FUNCTIONS RELATED TO HMM #####

def create_hmm(aln_file, out_file):
    """ Create hmm file for an alinment file """
    hmmbuild_cmd = 'hmmbuild %s %s' % (out_file, aln_file)
    try:
        subprocess.run(
            [hmmbuild_cmd], shell = True,
            stdout = open(os.devnull, 'wb'), check = True)
    except:
        sys.exit('ERROR: Error running hmmbuild %s' % hmmbuild_cmd)


def do_hmmsearch(
        hmm_file, search_fa, evalHmmer = 0.00001, scoreCutoff = 10,
        hitLimit = 10, cpus = os.cpu_count()):
    """ Perform hmmsearch for a hmm file vs a multiple fasta file
    Return a dictionary of hits and their e-value and bit-score
    Only "top" hits are returned. The cutoff is defined by
    max_score / 100 * (100 - scoreCutoff)
    By default, only hits that have at least 90% of the best bit score
    are considers
    """
    hmm_hits = {}
    with pyhmmer.easel.SequenceFile(search_fa, digital = True, alphabet = pyhmmer.easel.Alphabet.amino()) as seq_file:
        sequences = list(seq_file)
    with pyhmmer.plan7.HMMFile(hmm_file) as hmm_file:
        try:
            for hits in pyhmmer.hmmsearch(
                    hmm_file, sequences, E = evalHmmer, cpus = cpus):
                if len(hits) > 0:
                    n = 0
                    for hit in hits:
                        if hit.score >= hits[0].score/100*(100-scoreCutoff):
                            if n < hitLimit:
                                hmm_hits[hit.name.decode('ASCII')] = (
                                    hit.evalue,hit.score)
                                n += 1
        except:
            sys.exit(
                'ERROR: Error running hmmsearch for %s agains %s'
                % (hmm_file, search_fa))
    return(hmm_hits)
