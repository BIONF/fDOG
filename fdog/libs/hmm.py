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
    hmmbuild_cmd = 'hmmbuild --amino %s %s' % (out_file, aln_file)
    try:
        subprocess.run(
            [hmmbuild_cmd], shell = True,
            stdout = open(os.devnull, 'wb'), check = True)
    except:
        sys.exit('ERROR: Error running hmmbuild %s' % hmmbuild_cmd)


def sort_hmm_hits(hmm_hits, hitLimit = 10, scoreCutoff = 10):
    """ Sort HMM hits
    Keep only n hits (n =< hitLimit), and hits that are not less than
    best_hit_domain_score * (100 - scoreCutoff) / 100
    Input hmm_hits is a pyhmmer.plan7.topHits object
    """
    best_score = 0
    score_dict = {}
    for hit in hmm_hits:
        if len(hit.domains) > 0:
            domain_score = hit.domains[0].score
            hit_id = hit.domains[0].hit.name.decode('ASCII')
            if domain_score > best_score:
                best_score = domain_score
            if domain_score >= best_score/100*(100-scoreCutoff):
                if domain_score not in score_dict:
                    score_dict[domain_score] = [hit_id]
                else:
                    score_dict[domain_score].append(hit_id)

    hmm_cand = {}
    n = 1
    score_dict = {
        key:val for key, val in score_dict.items() \
        if key >= best_score/100*(100-scoreCutoff)
    }
    for score in sorted(score_dict):
        if n < hitLimit:
            for id in score_dict[score]:
                hmm_cand[id] = score
                n += 1
    return(hmm_cand)


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
                    hmm_hits = sort_hmm_hits(hits, hitLimit, scoreCutoff)
                    # n = 0
                    # for hit in hits:
                    #     if hit.score >= hits[0].score/100*(100-scoreCutoff):
                    #         if n < hitLimit:
                    #             hmm_hits[hit.name.decode('ASCII')] = (
                    #                 hit.evalue,hit.score)
                    #             n += 1
        except:
            sys.exit(
                'ERROR: Error running hmmsearch for %s agains %s'
                % (hmm_file, search_fa))
    return(hmm_hits)
