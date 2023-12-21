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

import fdog.libs.output as output_fn

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


def sort_hmm_hits(hmm_hits, hmm_score_type = 'domain', hitLimit = 10, scoreCutoff = 10, debug = False):
    """ Sort HMM hits
    Keep only n hits (n =< hitLimit), and hits that are not less than
    best_hit_domain_score * (100 - scoreCutoff) / 100
    Input hmm_hits is a pyhmmer.plan7.topHits object
    """
    best_score = -9999 # some "best" domains still have negative score!
    cutoff = ''
    score_dict = {}
    ori_hits = {}
    best_hit_score = -9999
    if hmm_score_type == 'domain':
        for hit in hmm_hits:
            ori_hits[hit.name.decode('ASCII')] = len(hit.domains)
            best_domain_score = -9999 #hit.domains[0].score
            best_domain_hit = ''
            if len(hit.domains) > 0:
                # get domain with best score for this hit
                for i in hit.domains:
                    if i.score > best_domain_score:
                        best_domain_score = i.score
                        best_domain_hit = i.hit.name.decode('ASCII')
                # add hit to score_dict with increasing domain score
                if best_domain_score > best_score:
                    best_score = best_domain_score
                cutoff = best_score/100*(100-scoreCutoff)
                if best_score < 0:
                    cutoff = best_score/100*(100+scoreCutoff)
                if best_domain_score >= cutoff:
                    if best_domain_score not in score_dict:
                        score_dict[best_domain_score] = [best_domain_hit]
                    else:
                        score_dict[best_domain_score].append(best_domain_hit)
    else:
        for hit in hmm_hits:
            ori_hits[hit.name.decode('ASCII')] = hit.score
            if hit.score > best_hit_score:
                # get hit with best score
                best_hit_score = hit.score
                best_hit = hit.name.decode('ASCII')
            # add to score_dict
            if best_hit_score > best_score:
                best_score = best_hit_score
            cutoff = best_score/100*(100-scoreCutoff)
            if best_score < 0:
                cutoff = best_score/100*(100+scoreCutoff)
            if hit.score >= cutoff:
                if hit.score not in score_dict:
                    score_dict[hit.score] = [hit.name.decode('ASCII')]
                else:
                    score_dict[hit.score].append(hit.name.decode('ASCII'))

    output_fn.print_debug(debug, 'All HMM hits', ori_hits)
    hmm_cand = {}
    n = 0
    score_dict = {
        key:val for key, val in score_dict.items() \
        if key >= cutoff
    }
    output_fn.print_debug(debug, 'Candidate HMM hits', score_dict)
    for score in sorted(score_dict, reverse = True):
        if n < hitLimit:
            for id in score_dict[score]:
                hmm_cand[id] = score
                n += 1
    return(hmm_cand)


def do_hmmsearch(
        hmm_file, search_fa, evalHmmer = 0.00001, scoreCutoff = 10,
        hitLimit = 10, hmm_score_type = 'domain', cpus = os.cpu_count(), debug = False):
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
                    hmm_hits = sort_hmm_hits(hits, hmm_score_type, hitLimit, scoreCutoff, debug)
        except :
            sys.exit(
                'ERROR: Error running hmmsearch for %s agains %s'
                % (hmm_file, search_fa))
    return(hmm_hits)
