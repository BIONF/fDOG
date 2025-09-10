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
import shutil
from pathlib import Path
from ete3 import NCBITaxa
from Bio import SeqIO
import time

import fdog.libs.zzz as general_fn
import fdog.libs.fasta as fasta_fn
import fdog.libs.hmm as hmm_fn
import fdog.libs.alignment as align_fn
import fdog.libs.tree as tree_fn
import fdog.libs.fas as fas_fn
import fdog.libs.output as output_fn
import fdog.libs.orthosearch as ortho_fn


##### FUNCTIONS RELATED TO CORE COMPILATION #####

def get_core_taxa_ids(coreTaxa, corepath):
    """ Get taxonomy IDs for core taxa
    Either from coreTaxa_dir, or from user input list (--coreTaxa)
    Return dictionary {taxID:<TaxName>@<TaxID>@Ver}
    """
    tax_ids = {}
    if not coreTaxa == '':
        ignored_taxa = []
        if os.path.exists(os.path.abspath(coreTaxa)):
            core_taxa = general_fn.read_file(coreTaxa)
        else:
            core_taxa = coreTaxa.split(',')

        for core_taxon in core_taxa:
            if not os.path.exists(
                    os.path.abspath(
                        '%s/%s/%s.phr' % (corepath,core_taxon,core_taxon))):
                ignored_taxa.append(core_taxon)
            else:
                id = core_taxon.split('@')[1]
                if not id in tax_ids:
                    tax_ids[id] = core_taxon
        if len(ignored_taxa) > 0:
            print(
                'WARNING: %s taxa cannot be found at %s\n%s'
                % (len(ignored_taxa), corepath, ignored_taxa))
    else:
        tax_ids = general_fn.get_ids_from_folder(corepath, 'coreTaxa_dir')
    return(tax_ids)


def initiate_core_files(
        seqFile, seqName, refspec, seed_id, hmmpath, annopath, aligner, fasOff):
    hmm_dir = '%s/%s/hmm_dir' % (hmmpath, seqName)
    Path(hmm_dir).mkdir(parents = True, exist_ok = True)
    aln_file = '%s/%s/hmm_dir/%s.aln' % (hmmpath, seqName, seqName)
    aln_seed = align_fn.do_align(aligner, seqFile)
    fasta_fn.write_fasta(aln_seed, aln_file)
    hmm_file = '%s/%s/hmm_dir/%s.hmm' % (hmmpath, seqName, seqName)
    hmm_seed = hmm_fn.create_hmm(aln_file, hmm_file)

    fa_file = '%s/%s/%s.fa' % (hmmpath, seqName, seqName)
    seed_id_mod = '%s|%s|%s' % (seqName, refspec, seed_id)
    input_seed = SeqIO.parse(seqFile,'fasta')
    with open(fa_file, 'w') as initial_core_fa:
        for fa in input_seed:
            initial_core_fa.write('>%s\n%s\n' % (seed_id_mod, str(fa.seq)))

    seed_json = ''
    if not fasOff == True:
        seed_json = fas_fn.get_anno_fas(
            seqName, refspec, seed_id, str(fa.seq), hmmpath, annopath)
    return(aln_file, fa_file, hmm_file, seed_json)


def store_cand_reults(args):
    """ Save intermediate results for a candidate ortholog
    Including:
    1) Candidate joined score of fas & normalised aln score
    in dictionary {taxID:score}
    2) Candidate fasta sequence in dictionary {taxID:fasta_seq}
    3) Update current (best) candidate, normalised aln score and joined score
    """
    (cand_taxid, cand_score, cand_seq,
        curr_cand, curr_aln_score, aln_score_normalized,
        fas_score, fas_dict, ortho_id, ortho_seq) = args
    cand_score[cand_taxid] = float(fas_score) + float(aln_score_normalized)
    cand_seq[cand_taxid] = {ortho_id:ortho_seq}
    curr_aln_score = aln_score_normalized
    curr_candi_score = float(fas_score) + float(aln_score_normalized)
    curr_cand = cand_taxid
    return(
        cand_score, cand_seq, curr_cand, curr_aln_score,
        curr_candi_score, fas_dict)


def validate_candidate(args):
    """ Validate candidate based on its normalised aln score and fas score """
    (aln_score_normalized, cand_args, calc_fas_args, variable_args, debugCore, distDeviation) = args
    (cand_score, cand_seq, curr_cand, curr_aln_score,
                        curr_candi_score, fas_dict) = variable_args
    (cand_taxid, ortho_id, ortho_seq, next_node, first_cand) = cand_args
    (fasOff, seqName, seed_json, spec, seq_id, seq, hmmpath, annopath) = calc_fas_args

    if first_cand == True:
        threshold = 0
    else:
        if next_node == True:
            threshold = curr_candi_score * (1 + distDeviation)
            type = ' (STRICT) '
        else:
            threshold = curr_candi_score
            type = ''

    if aln_score_normalized > threshold - 1:
        if not '%s_%s' % (spec, seq_id) in fas_dict:
            fas_score = fas_fn.calc_fas_cand(calc_fas_args)
            fas_dict['%s_%s' % (spec, seq_id)] = fas_score
            output_fn.print_debug(
                debugCore, '', '-FAS: %s' % fas_score)
        else:
            fas_score = fas_dict['%s_%s' % (spec, seq_id)]
        if float(fas_score) + float(aln_score_normalized) > threshold:
            variable_args = store_cand_reults(
                        [cand_taxid,
                        cand_score, cand_seq,
                        curr_cand, curr_aln_score,
                        aln_score_normalized, fas_score, fas_dict,
                        ortho_id, ortho_seq])
        else:
            output_fn.print_debug(
                debugCore, '',
                '-Joined score not higher than the prev%s! Skip...' % type)
    else:
        output_fn.print_debug(
            debugCore, '',
            '-Aln score %s not higher%s! Skip...' % (aln_score_normalized, type))
    return(variable_args)


def compile_core(args):
    """ Core compilation """
    (seqFile, seqName, refspec, seed_id, coreArgs, pathArgs, orthoArgs, otherArgs, debug) = args
    (minDist, maxDist, coreSize, coreTaxa, distDeviation, alnStrategy, fasOff) = coreArgs
    (outpath, hmmpath, corepath, searchpath, annopath) = pathArgs
    (cpus, debugCore, silentOff, noCleanup, force, append) = otherArgs
    aligner = orthoArgs[-1]
    otherArgs.insert(0, 'NA')

    ncbi = NCBITaxa()
    ### get taxonomy lineage of refspec
    refspec_id = refspec.split('@')[1]
    refspec_lineage = ncbi.get_lineage(refspec_id)

    ### get rank ID and its index in the refspec lineage
    (min_rank, max_rank) = tree_fn.get_rank_range(refspec_lineage, minDist, maxDist, ncbi)
    output_fn.print_debug(debugCore, 'Min & Max-rank', '%s\t%s' % (min_rank, max_rank))

    ### create taxonomy tree from list of core tax
    tax_ids = get_core_taxa_ids(coreTaxa, corepath)
    tree = ncbi.get_topology(tax_ids.keys(), intermediate_nodes = True)
    if debugCore:
        print(tree)

    ### INITIATE FA, ALN, HMM [and anno FAS] FILE FOR SEED
    (aln_file, fa_file, hmm_file, seed_json) = initiate_core_files(
        seqFile, seqName, refspec, seed_id, hmmpath, annopath, aligner, fasOff)

    ### get list of taxa within min and max rank of refspec
    node_dict = tree_fn.get_leaves_dict(
                    refspec_lineage, tree,
                    list(min_rank.values())[0], list(max_rank.values())[0])
    output_fn.print_debug(debugCore, 'Node dictionary', node_dict)

    ### traverse the core taxa tree
    added_taxa = {}
    ignored_taxa = []
    fas_dict = {}
    previous_added_taxon = refspec_id
    for round in range(coreSize - 1):
        output_fn.print_stdout(silentOff, '---------- ROUND %s ----------' % round)
        output_fn.print_debug(
            debugCore, 'CORE COMPILATION',
            '---------- ROUND %s ----------' % round)
        aln_scores = align_fn.calc_aln_score(fa_file, fa_file, alnStrategy, debugCore)
        max_aln_score = max(aln_scores.values()) #0
        if max_aln_score == 0:
            exit('ERROR: Something went wrong with FASTA36. Please run debugCore to investigate!')
        flag_round = 0 # use to stop current round if an ortholog was added
        output_fn.print_debug(debugCore, '', 'ADDED TAXA: %s' % added_taxa.keys())
        cand_seq = {}
        cand_score = {}
        curr_cand = ''
        curr_aln_score = 0
        curr_candi_score = 0
        next_node = False
        for node_id, leaves in node_dict.items():
            if flag_round == 1:
                break
            output_fn.print_debug(
                debugCore, '',
                'NODE %s - %s' % (node_id, ncbi.get_rank([node_id])))
            output_fn.print_debug(
                debugCore, '',
                '-MAX ALN SCORE: %s' % max_aln_score)
            output_fn.print_debug(
                debugCore, '',
                '-PREVIOUS ADDED: %s' % previous_added_taxon)
            leaves.reverse()
            flag_node = 0
            for leaf in leaves:
                if not leaf in tax_ids:
                    continue
                if flag_node == 1:
                    break
                if not leaf == refspec_id and \
                        not leaf in added_taxa and \
                        not leaf in ignored_taxa:
                    output_fn.print_debug(debugCore, '','')
                    output_fn.print_debug(
                        debugCore, '',
                        'Leaf %s - %s' % (leaf, tax_ids[leaf]))
                    if len(curr_cand) > 0 and \
                            not curr_cand in node_dict[node_id]:
                        next_node = True
                        output_fn.print_debug(
                            debugCore, '',
                            '-Current_candidate from different node: %s' \
                            % curr_cand)
                        if curr_candi_score * (1 + distDeviation) > 2:
                            output_fn.print_debug(
                                debugCore, '',
                                '# Current score cannot be defeater! Stop this node!')
                            break
                    else:
                        next_node = False
                        output_fn.print_debug(
                            debugCore, '',
                            '-Current_candidate: %s' % curr_cand)
                    output_fn.print_debug(
                        debugCore, '',
                        '-Current_aln_score: %s' % curr_aln_score)
                    output_fn.print_debug(
                        debugCore, '',
                        '-Current_candi_score: %s' % curr_candi_score)
                    ### compare taxonomy rank with previous added taxon
                    ### ignore if this leaf closer to the refspec than to
                    ### previous added taxon
                    ancestor_to_ref = tree_fn.get_ancestor(refspec_id, leaf, ncbi)
                    check_ancestor_to_ref = tree_fn.check_common_ancestor(
                            previous_added_taxon, list(ancestor_to_ref.keys())[0],
                            minDist, maxDist, ncbi)
                    if check_ancestor_to_ref == 0:
                        ignored_taxa.append(leaf)
                        output_fn.print_debug(
                            debugCore, '',
                            '-Closer to refspec than previous added taxon!')
                        continue
                    ### continue process only if this leaf is within min and max rank
                    ### of the previous added taxon
                    ancestor = tree_fn.get_ancestor(previous_added_taxon, leaf, ncbi)
                    check_ancestor = tree_fn.check_common_ancestor(
                            previous_added_taxon, list(ancestor.keys())[0],
                            minDist, maxDist, ncbi)
                    if check_ancestor == 1:
                        output_fn.print_debug(
                            debugCore, '',
                            '-Ancestor %s with %s accepted' \
                            % (ancestor, previous_added_taxon))
                        ### run ortholog search
                        otherArgs[2] = debug
                        otherArgs[0] = tax_ids[leaf]
                        hamstr_out = ortho_fn.run_hamstr([seqName, refspec, pathArgs,
                            orthoArgs, otherArgs])
                        if len(hamstr_out) > 1:
                            ### calculate alignment score
                            ortho = list(hamstr_out.items())[-1]
                            tmp_fa = '%s/%s/%s_%s.fa' \
                                % (hmmpath, seqName, seqName, leaf)
                            with open(tmp_fa, 'w') as tmp_fa_out:
                                tmp_fa_out.write('>%s\n%s\n' \
                                    % (ortho[0][0:len(ortho[0])-2], ortho[1]))
                            aln_score = align_fn.calc_aln_score(fa_file, tmp_fa, alnStrategy, debugCore)
                            output_fn.print_debug(
                                debugCore, '',
                                '-Max: %s - Aln: %s' % (max_aln_score, aln_score))
                            aln_score_normalized = \
                                    list(aln_score.values())[0] / max_aln_score
                            output_fn.print_debug(
                                debugCore, '',
                                '-Normalized_aln_score: %s' % aln_score_normalized)
                            os.remove(tmp_fa)

                            ### validate candidate
                            if len(cand_score) == 0 \
                                    and len(curr_cand) == 0:
                                first_cand = True
                            else:
                                first_cand = False
                            calc_fas_args = (fasOff, seqName, seed_json,
                                tax_ids[leaf], ortho[0].split('|')[-2],
                                ortho[1] ,hmmpath, annopath)
                            cand_args = (leaf, ortho[0][0:len(ortho[0])-2],
                                ortho[1], next_node, first_cand)
                            variable_args = (cand_score, cand_seq,
                                curr_cand, curr_aln_score,
                                curr_candi_score, fas_dict)
                            (cand_score, cand_seq, curr_cand,
                                curr_aln_score, curr_candi_score,
                                fas_dict) = validate_candidate([
                                        aln_score_normalized, cand_args,
                                        calc_fas_args, variable_args, debugCore,
                                        distDeviation])
                            if curr_candi_score == 2:
                                flag_node = 1
                                output_fn.print_debug(
                                    debugCore, '',
                                    '-Max score achieved! Stop this node!')
                        elif not len(hamstr_out) > 1 and len(cand_seq) == 0:
                            ignored_taxa.append(leaf)
                            output_fn.print_debug(
                                debugCore, '',
                                '-No ortholog found!')
                    else:
                        ignored_taxa.append(leaf)
                        output_fn.print_debug(
                            debugCore, '',
                            '-Not considered due to ancestor %s with %s\n' \
                                % (ancestor, previous_added_taxon))
                else:
                    output_fn.print_debug(
                        debugCore, '',
                        '%s - %s skipped' % (leaf, tax_ids[leaf]))
            output_fn.print_debug(
                debugCore, '', 'Node candidates: %s' % cand_score)
            if len(cand_score) > 0 \
                    and cand_score[curr_cand] == 2:
                output_fn.print_debug(
                    debugCore, '',
                    '# MAX SCORE ACCHIEVED! Stop this round!')
                flag_round = 1
            if next_node == True \
                    and cand_score[curr_cand] * (1 + distDeviation) > 2:
                output_fn.print_debug(
                    debugCore, '',
                    '# CURRENT SCORE CANNOT BE DEFEATED! Stop this round!')
                flag_round = 1

        if len(cand_seq) > 0:
            output_fn.print_debug(
                debugCore, '',
                '# ADD THIS TAXON TO CORE GROUP\t%s - %s\n' \
                % (curr_cand, tax_ids[curr_cand]))
            previous_added_taxon = curr_cand
            added_taxa[curr_cand] = {tax_ids[curr_cand]:cand_score[curr_cand]}
            ### update seqName.fa and hmm_dir/seqName.hmm
            fasta_fn.append_to_fasta_file(fa_file, cand_seq[curr_cand])
            aln_seed = align_fn.do_align(aligner, fa_file)
            fasta_fn.write_fasta(aln_seed, aln_file)
            hmm_seed = hmm_fn.create_hmm(aln_file, hmm_file)
            os.remove(aln_file)
    ### remove temp json files
    for file in os.listdir('%s/%s' % (hmmpath, seqName)):
        if file.endswith('.json'):
            os.remove('%s/%s/%s' % (hmmpath, seqName, file))
    output_fn.print_debug(
        debugCore, 'CORE COMPILATION',
        'All added taxa %s' % added_taxa)
    if len(added_taxa) < coreSize - 1:
        output_fn.print_stdout(
            silentOff,
            'WARNING: Only %s/%s orthologs in the core group' \
            % (len(added_taxa) + 1, coreSize))


def run_compile_core(args):
    (seqFile, seqName, refspec, seed_id, reuseCore, forceCore, coreArgs,
            pathArgs, orthoCoreArgs, otherCoreArgs, debug) = args
    (outpath, hmmpath, corepath, searchpath, annopath) = pathArgs
    (cpus, debugCore, silentOff, noCleanup, force, append) = otherCoreArgs[-6:]
    begin = time.time()
    fdogPath = os.path.realpath(__file__).replace('/libs/corecompile.py','')
    align_fn.check_fasta36_executable(fdogPath)

    coreHmmfile = '%s/%s/hmm_dir/%s.hmm' % (hmmpath, seqName, seqName)
    coreHmmfile = os.path.abspath(coreHmmfile)
    compile_core_check = 1
    ncbi = ''
    if reuseCore == True:
        general_fn.check_file_exist(coreHmmfile)
        compile_core_check = 0
    else:
        if os.path.exists(coreHmmfile):
            if forceCore == True:
                print('WARNING: Existing %s core group will be deleted!' % seqName)
                shutil.rmtree('%s/%s' % (hmmpath, seqName))
            else:
                sys.exit(
                    'WARNING: Core group %s exists in %s! ' % (seqName, hmmpath)
                    + 'You still can run with --forceCore or --reuseCore option')
    if compile_core_check == 1:
        compile_core([seqFile, seqName, refspec, seed_id, coreArgs, pathArgs,
                    orthoCoreArgs, otherCoreArgs[-6:], debug])
    end = time.time()
    return([seqName, '{:5.3f}s'.format(end - begin)])
