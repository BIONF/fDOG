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

import os
from Bio import SeqIO
import multiprocessing as mp
from tqdm import tqdm
import time

import fdog.libs.zzz as general_fn
import fdog.libs.fasta as fasta_fn
import fdog.libs.blast as blast_fn
import fdog.libs.hmm as hmm_fn
import fdog.libs.alignment as align_fn
import fdog.libs.preparation as prepare_fn
import fdog.libs.output as output_fn


##### REVERSE SEARCH OF ONE pHMM HIT (step 4 of HaMStR) #####
def reverse_check(hmm_hit, seed_ids, seqName, search_taxon, outpath,
        refspec_db, refspec_seqs, search_seqs,
        evalBlast, lowComplexityFilter, checkCoorthologsRefOff,
        aligner, cpus, debug, silentOff, noCleanup):
    """ Reverse BLAST search of one pHMM hit against the reference proteome
    seed_ids: reference sequences that may serve as seed, in order of preference.
              Default (one seed): [first reference sequence of the core group].
              With --matchedSeed: all reference sequences of the core group.
    Return the seed ID the hit was accepted with, or None if it was rejected:
    (4a) the best BLAST hit is one of the seeds -> accepted with that seed
    (4b) otherwise (co-ortholog test): accepted with the first seed for which the
         best BLAST hit is closer to the seed than to the hmm hit
    """
    hmm_hit_fa = '%s/hmm_%s_%s_%s.fa' % (
                        outpath, seqName, search_taxon, hmm_hit)
    with open(hmm_hit_fa, 'w') as hmm_fa_out:
        hmm_fa_out.write('>%s_%s\n%s' % (search_taxon, hmm_hit, search_seqs.fetch(hmm_hit)))
    blast_xml = blast_fn.do_blastsearch(
            hmm_hit_fa, refspec_db, evalBlast = evalBlast, lowComplexityFilter = lowComplexityFilter, cpus = cpus)
    blast_out = blast_fn.parse_blast_xml(blast_xml)
    output_fn.print_debug(debug, 'BLAST hits', blast_out)
    if noCleanup == False:
        os.remove(hmm_hit_fa)
    if not len(blast_out['hits'].keys()) > 0:
        return(None)
    best_blast_hit = list(blast_out['hits'].keys())[0]
    if best_blast_hit == hmm_hit and len(blast_out['hits'].keys()) > 1:
        best_blast_hit = list(blast_out['hits'].keys())[0]
    ### (4a) if refspec_seq_id == best blast hit
    if best_blast_hit in seed_ids:
        output_fn.print_stdout(
            silentOff,
            '%s accepted (best blast hit is ref)' % (blast_out['query']))
        return(best_blast_hit)
    ### (4b) else, check for co-ortholog ref
    if checkCoorthologsRefOff == False:
        for seed_id in seed_ids:
            aln_fa = '%s/blast_%s_%s_%s_%s_%s.fa' % (
                        outpath, seqName, seed_id, search_taxon,
                        hmm_hit, best_blast_hit)
            with open(aln_fa, 'w') as aln_fa_out:
                aln_fa_out.write(
                    '>%s\n%s\n>%s_%s\n%s\n>%s\n%s' % (
                        seed_id, refspec_seqs.fetch(seed_id),
                        search_taxon, hmm_hit, search_seqs.fetch(hmm_hit),
                        best_blast_hit, refspec_seqs.fetch(best_blast_hit)
                    )
                )
            fasta_fn.remove_dup(aln_fa)
            aln_seq = align_fn.do_align(aligner, aln_fa)
            output_fn.print_debug(
                debug, 'Alignment for checking co-ortholog ref', aln_seq)
            br_dist = align_fn.calc_Kimura_dist(aln_seq, best_blast_hit, seed_id, debug)
            bh_dist = align_fn.calc_Kimura_dist(aln_seq, best_blast_hit, f'{search_taxon}_{hmm_hit}', debug)
            output_fn.print_debug(
                debug, 'Check if distance blast_vs_ref < blast_vs_hmm',
                'd_br = %s; d_bh = %s' % (br_dist, bh_dist))
            if noCleanup == False:
                os.remove(aln_fa)
            if br_dist == bh_dist == 0 or br_dist < bh_dist:
                output_fn.print_stdout(
                    silentOff,
                    '%s accepted (best blast hit is co-ortholog to ref)'
                    % (blast_out['query'])
                )
                return(seed_id)
    return(None)


##### FUNCTION FOR HMM-BASED ORTHOLOG SEARCH (HaMStR) #####
def hamstr(args):
    (seqName, hmmpath, corepath, searchpath, outpath,
        refspec, seed_id, search_taxon,
        evalHmmer, hitLimit, hmmScoreType, scoreCutoff,
        evalBlast, lowComplexityFilter,
        checkCoorthologsRefOff, rbh, rep,
        aligner, cpus, debug, silentOff, noCleanup,
        seed_ids, windowFallback) = args
    """ Ortholog search algorithm for a hmm core group agains a search taxon
    Implemented based on HaMStR https://doi.org/10.1186/1471-2148-9-157
    seed_ids: reference sequences that may serve as seed ([seed_id] by default,
              all reference sequences of the core group with --matchedSeed)
    windowFallback: if no hit within the score window is accepted, check the
              other hits among the hitLimit best hits (--windowFallback)
    Return (ortho_final, seed_info). seed_info records which seed was used for
    the reported ortholog(s), its position in seed_ids, and whether it was
    found only by the fallback.
    """
    ### (0) Dict for storing candidate and final orthologs (key=id, value=seq)
    ortho_candi = {}
    candi_seed = {} # seed each accepted candidate was accepted with (key=id, value=seed id)
    ortho_final = {}
    seed_info = {'search_taxon': search_taxon, 'seed': None, 'seed_rank': None, 'fallback': False}
    ### (00) Parse input files
    hmm_file = '%s/%s/hmm_dir/%s.hmm' % (hmmpath, seqName, seqName)
    refspec_db = '%s/%s/%s' % (corepath, refspec, refspec)
    refspec_fa = '%s/%s/%s.fa' % (corepath, refspec, refspec)
    search_fa = '%s/%s/%s.fa' % (searchpath, search_taxon, search_taxon)
    ### (000) Adapt parameters
    if rbh == True:
        checkCoorthologsRefOff = True
        rep = True

    ### PRINT JOB PARAMETERS
    output_fn.print_stdout(
        silentOff,
        '\n### Ortholog search ###'
        + '\nSeed: %s\nRefspec: %s\n' % (seqName, refspec)
        + 'Ref_seqID: %s\n' % ','.join(seed_ids)
        + 'Search taxon: %s' % search_taxon)
    output_fn.print_debug(
        debug, 'Parameters',
        'HMM evalue cutoff: %s\nHMM hit limit: %s\n' % (evalHmmer, hitLimit)
        + 'HMM hit score cutoff: %s\n' % scoreCutoff
        + 'BLAST evalue cutoff: %s\n' % evalBlast
        + 'Blast low complexity filter: %s\n' % lowComplexityFilter
        + 'Turn off check for co-orthologs ref: %s\n' % checkCoorthologsRefOff
        + 'Window fallback: %s\n' % windowFallback
        + 'Aligner: %s' % aligner)

    ### (1) Do hmmsearch for query hmm against search taxon fasta
    if windowFallback == True:
        (hmm_hits, reserve_hits) = hmm_fn.do_hmmsearch_with_reserve(
            hmm_file, search_fa, evalHmmer, scoreCutoff, hitLimit, hmmScoreType, cpus, debug)
    else:
        hmm_hits = hmm_fn.do_hmmsearch(
            hmm_file, search_fa, evalHmmer, scoreCutoff, hitLimit, hmmScoreType, cpus, debug)
        reserve_hits = {}
    output_fn.print_debug(debug, 'Sorted HMM hits', hmm_hits)
    ### (2) Read fasta file of refspec and search taxon
    refspec_seqs = fasta_fn.read_fasta(refspec_fa)
    search_seqs = fasta_fn.read_fasta(search_fa)
    ### (3) Do re-blast search for each hmm hit against refspec
    if len(hmm_hits) == 0:
        output_fn.print_stdout(
            silentOff, 'WARNING: No HMM hit found!')
    else:
        ### (4) check reciprocity; first for the hits within the score window,
        ### then (only with --windowFallback, and only if nothing was accepted)
        ### for the other hits among the hitLimit best hits
        for phase, phase_hits in enumerate([hmm_hits, reserve_hits]):
            if phase == 1:
                if len(ortho_candi) > 0 or len(phase_hits) == 0:
                    break
                seed_info['fallback'] = True
                output_fn.print_stdout(
                    silentOff,
                    'No hit within the score window accepted. Checking %s further hit(s)...'
                    % len(phase_hits))
                output_fn.print_debug(debug, 'Further HMM hits (fallback)', phase_hits)
            for hmm_hit in phase_hits:
                accepted_seed = reverse_check(
                    hmm_hit, seed_ids, seqName, search_taxon, outpath,
                    refspec_db, refspec_seqs, search_seqs,
                    evalBlast, lowComplexityFilter, checkCoorthologsRefOff,
                    aligner, cpus, debug, silentOff, noCleanup)
                if accepted_seed:
                    ortho_candi[hmm_hit] = search_seqs.fetch(hmm_hit)
                    candi_seed[hmm_hit] = accepted_seed
            # remove seed protein from candidate list
            # (a reference sequence that was accepted with itself as seed)
            if search_taxon == refspec:
                for cand in list(ortho_candi.keys()):
                    if cand == candi_seed[cand]:
                        ortho_candi.pop(cand)
                        candi_seed.pop(cand)

        ### (5) check co-ortholog if more than 1 HMM hits are accepted
        if len(ortho_candi) == 0:
            seed_info['fallback'] = False
            output_fn.print_stdout(
                silentOff, 'WARNING: Reciprocity not fulfulled! No ortholog found!')
        else:
            output_fn.print_debug(
                debug, 'Candidates for checking co-orthologs', ortho_candi.keys())
            best_ortho = list(ortho_candi.keys())[0]
            # the seed for step 5 (and for the output) is the seed that the
            # best ortholog was accepted with
            seed_used = candi_seed[best_ortho]
            seed_info['seed'] = seed_used
            seed_info['seed_rank'] = seed_ids.index(seed_used) + 1
            ortho_final = fasta_fn.add_seq_to_dict(
                ortho_final, '%s|%s|%s|1' % (seqName, search_taxon, best_ortho),
                ortho_candi[best_ortho])

            if rep == False:
                if len(ortho_candi) > 1:
                    aln_co_fa = '%s/coortho_%s_%s.fa' % (
                                outpath, seqName, search_taxon)
                    with open(aln_co_fa, 'w') as aln_co_fa_out:
                        aln_co_fa_out.write(('>%s\n%s\n') %
                            (seed_used, refspec_seqs.fetch(seed_used)))
                        for cand in ortho_candi:
                            aln_co_fa_out.write(('>%s\n%s\n') %
                                (f'{search_taxon}_{cand}', ortho_candi[cand]))
                    aln_co_seq = align_fn.do_align(aligner, aln_co_fa)
                    output_fn.print_debug(
                        debug, 'Alignment for checking co-orthologs', aln_co_seq)
                    if noCleanup == False:
                        os.remove(aln_co_fa)
                    best_dist = align_fn.calc_Kimura_dist(
                                    aln_co_seq, seed_used, f'{search_taxon}_{best_ortho}', debug)
                    for cand in ortho_candi:
                        if not cand == best_ortho:
                            candi_dist = align_fn.calc_Kimura_dist(
                                aln_co_seq, f'{search_taxon}_{best_ortho}', f'{search_taxon}_{cand}', debug)
                            output_fn.print_debug(
                                debug,
                                'Check if distance bestHmm_vs_ref > '
                                + 'other_vs_bestHmm',
                                'd_best = %s; d_other = %s'
                                % (best_dist, candi_dist))
                            if candi_dist < best_dist:
                                if not cand == seed_used:
                                    ortho_final = fasta_fn.add_seq_to_dict(
                                        ortho_final,
                                        '%s|%s|%s|0' \
                                            % (seqName, search_taxon, cand),
                                        ortho_candi[cand])
            output_fn.print_stdout(
                silentOff,
                '=> %s orthologs found: %s'
                % (len(ortho_final), list(ortho_final.keys())))
    return(ortho_final, seed_info)


def run_hamstr(args, matchedSeed = False, windowFallback = False, return_seed_info = False):
    """ Perform ortholog search based on hamstr approach
    matchedSeed: accept a candidate whose best reverse BLAST hit is any of the
        reference sequences in the core group; that sequence is then the seed
        (only relevant if the reference species has several sequences in the
        core group, i.e. with --reuseCore)
    windowFallback: if no hit within the score window is accepted, check the
        other hits among the hitLimit best hits
    return_seed_info: also return a list with one seed_info dict per search taxon
    """

    (seqName, refspec, pathArgs, orthoArgs, otherArgs) = args
    (outpath, hmmpath, corepath, searchpath, annopath) = pathArgs
    (checkCoorthologsRefOff, rbh, rep, evalBlast, lowComplexityFilter,
                        evalHmmer, hitLimit, hmmScoreType, scoreCutoff, aligner) = orthoArgs
    (searchTaxa, cpus, debug, silentOff, noCleanup, force, append) = otherArgs

    hamstr_jobs = []
    ### get ref seqID
    core_fa = '%s/%s/%s.fa' % (hmmpath, seqName, seqName)
    seed_id = prepare_fn.get_seed_id_from_fa(core_fa, refspec)
    if matchedSeed == True:
        seed_ids = prepare_fn.get_seed_ids_from_fa(core_fa, refspec)
    else:
        seed_ids = [seed_id]

    ### get search taxa from user defined list (as a file or directly a list)
    if not searchTaxa == '':
        ignored_taxa = []
        if os.path.exists(os.path.abspath(searchTaxa)):
            search_taxa = general_fn.read_file(searchTaxa)
        else:
            search_taxa = searchTaxa.split(',')

        for search_taxon in search_taxa:
            if os.path.exists(
                    os.path.abspath(
                        '%s/%s/%s.fa' % (searchpath,search_taxon,search_taxon))):
                hamstr_jobs.append([
                    seqName, hmmpath, corepath, searchpath, outpath,
                    refspec, seed_id, search_taxon,
                    evalHmmer, hitLimit, hmmScoreType, scoreCutoff,
                    evalBlast, lowComplexityFilter,
                    checkCoorthologsRefOff, rbh, rep,
                    aligner, cpus, debug, silentOff, noCleanup,
                    seed_ids, windowFallback
                ])
            else:
                ignored_taxa.append(search_taxon)
        if len(ignored_taxa) > 0:
            print(
                'WARNING: %s taxa cannot be found at %s\n%s'
                % (len(ignored_taxa), searchpath, ignored_taxa))
    ### get search taxa from searchpath (searchTaxa_dir)
    else:
        for search_taxon in general_fn.read_dir(searchpath):
            if os.path.exists(
                    os.path.abspath(
                        '%s/%s/%s.fa' % (searchpath,search_taxon,search_taxon))):
                hamstr_jobs.append([
                    seqName, hmmpath, corepath, searchpath, outpath,
                    refspec, seed_id, search_taxon,
                    evalHmmer, hitLimit, hmmScoreType, scoreCutoff,
                    evalBlast, lowComplexityFilter,
                    checkCoorthologsRefOff, rbh, rep,
                    aligner, cpus, debug, silentOff, noCleanup,
                    seed_ids, windowFallback
                ])

    ### do ortholog search
    hamstr_out = {}
    seed_infos = []
    if len(hamstr_jobs) > 0:
        if cpus > len(hamstr_jobs):
            cpus = len(hamstr_jobs)
        output_fn.print_stdout(
                silentOff, 'Ortholog search for %s taxa...' % len(hamstr_jobs))
        if debug == True or silentOff == True or len(hamstr_jobs) == 1:
            for job in hamstr_jobs:
                (tmp_out, tmp_info) = hamstr(job)
                hamstr_out = {**hamstr_out, **tmp_out}
                seed_infos.append(tmp_info)
        else:
            pool = mp.Pool(cpus)
            for _ in tqdm(
                    pool.imap_unordered(hamstr, hamstr_jobs),
                    total=len(hamstr_jobs)):
                (tmp_out, tmp_info) = _
                if len(tmp_out) > 0:
                    hamstr_out = {**hamstr_out, **tmp_out}
                seed_infos.append(tmp_info)

    ### Get seed seq
    # the seed written into the output (and used by FAS as the reference) is the
    # seed that the orthologs were accepted with. By default this is always the
    # first reference sequence. If different search taxa were accepted with
    # different seeds (only possible with --matchedSeed), the first reference
    # sequence is written and a warning is printed.
    used_seeds = []
    for info in seed_infos:
        if info['seed'] and not info['seed'] in used_seeds:
            used_seeds.append(info['seed'])
    out_seed = seed_id
    if len(used_seeds) == 1:
        out_seed = used_seeds[0]
    elif len(used_seeds) > 1:
        print(
            'WARNING: Search taxa of %s were accepted with different reference sequences (%s). '
            'The output lists %s as seed sequence. See %s.seed_info.tsv for the seed of each taxon.'
            % (seqName, ','.join(used_seeds), seed_id, seqName))
    refspec_fa = '%s/%s/%s.fa' % (corepath, refspec, refspec)
    refspec_seqs = fasta_fn.read_fasta(refspec_fa)
    seed_id_mod = '%s|%s|%s|1' % (seqName, refspec, out_seed)
    seed_seq = refspec_seqs.fetch(out_seed)

    ### return
    if return_seed_info == True:
        return({**{seed_id_mod:seed_seq}, **hamstr_out}, seed_infos)
    return({**{seed_id_mod:seed_seq}, **hamstr_out})


def get_search_taxa_ids(searchTaxa, searchpath):
    """ Get taxonomy IDs for search taxa
    Either from searchTaxa_dir, or from user input list (--searchTaxa)
    Return dictionary {taxID:<TaxName>@<TaxID>@Ver}
    """
    tax_ids = {}
    if not searchTaxa == '' and not len(searchTaxa) == 0:
        ignored_taxa = []
        if os.path.exists(os.path.abspath(searchTaxa)):
            search_taxa = general_fn.read_file(searchTaxa)
        else:
            search_taxa = searchTaxa.split(',')

        for search_taxon in search_taxa:
            if not os.path.exists(
                    os.path.abspath(
                        '%s/%s/%s.fa' % (searchpath,search_taxon,search_taxon))):
                ignored_taxa.append(search_taxon)
            else:
                id = search_taxon.split('@')[1]
                if not id in tax_ids:
                    tax_ids[id] = search_taxon
        if len(ignored_taxa) > 0:
            print(
                'WARNING: %s taxa cannot be found at %s\n%s'
                % (len(ignored_taxa), searchpath, ignored_taxa))
    else:
        tax_ids = general_fn.get_ids_from_folder(searchpath, 'searchTaxa_dir')
    return(tax_ids)
