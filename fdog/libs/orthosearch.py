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


##### FUNCTION FOR HMM-BASED ORTHOLOG SEARCH (HaMStR) #####
def hamstr(args):
    (seqName, hmmpath, corepath, searchpath, outpath,
        refspec, seed_id, search_taxon,
        evalHmmer, hitLimit, scoreCutoff,
        evalBlast, lowComplexityFilter,
        checkCoorthologsRefOff, rbh, rep,
        aligner, cpus, debug, silentOff, noCleanup) = args
    """ Ortholog search algorithm for a hmm core group agains a search taxon
    Implemented based on HaMStR https://doi.org/10.1186/1471-2148-9-157
    """
    ### (0) Dict for storing candidate and final orthologs (key=id, value=seq)
    ortho_candi = {}
    ortho_final = {}
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
        + 'Ref_seqID: %s\n' % seed_id
        + 'Search taxon: %s' % search_taxon)
    output_fn.print_debug(
        debug, 'Parameters',
        'HMM evalue cutoff: %s\nHMM hit limit: %s\n' % (evalHmmer, hitLimit)
        + 'HMM hit score cutoff: %s\n' % scoreCutoff
        + 'BLAST evalue cutoff: %s\n' % evalBlast
        + 'Blast low complexity filter: %s\n' % lowComplexityFilter
        + 'Turn off check for co-orthologs ref: %s\n' % checkCoorthologsRefOff
        + 'Aligner: %s' % aligner)

    ### (1) Do hmmsearch for query hmm against search taxon fasta
    hmm_hits = hmm_fn.do_hmmsearch(
            hmm_file, search_fa, evalHmmer, scoreCutoff, hitLimit, cpus, debug)
    output_fn.print_debug(debug, 'Sorted HMM hits', hmm_hits)
    ### (2) Read fasta file of refspec and search taxon
    refspec_seqs = fasta_fn.read_fasta(refspec_fa)
    search_seqs = fasta_fn.read_fasta(search_fa)
    ### (3) Do re-blast search for each hmm hit against refspec
    for hmm_hit in hmm_hits:
        if not hmm_hit == seed_id: # only if search taxon == refspec
            hmm_hit_fa = '%s/hmm_%s_%s_%s.fa' % (
                                outpath, seqName, search_taxon, hmm_hit)
            with open(hmm_hit_fa, 'w') as hmm_fa_out:
                hmm_fa_out.write('>%s\n%s' % (hmm_hit, search_seqs.fetch(hmm_hit)))
            blast_xml = blast_fn.do_blastsearch(
                    hmm_hit_fa, refspec_db, evalBlast = evalBlast, lowComplexityFilter = lowComplexityFilter)
            blast_out = blast_fn.parse_blast_xml(blast_xml)
            output_fn.print_debug(debug, 'BLAST hits', blast_out)
            if noCleanup == False:
                os.remove(hmm_hit_fa)
            ### (4) check reciprocity
            ### (4a) if refspec_seq_id == best blast hit
            if len(blast_out['hits'].keys()) > 0:
                best_blast_hit = list(blast_out['hits'].keys())[0]
                if best_blast_hit == hmm_hit and len(blast_out['hits'].keys()) > 1:
                    best_blast_hit = list(blast_out['hits'].keys())[1]
                if seed_id == best_blast_hit:
                    output_fn.print_stdout(
                        silentOff,
                        '%s accepted (best blast hit is ref)' % (blast_out['query']))
                    ortho_candi[hmm_hit] = search_seqs.fetch(hmm_hit)
                    continue
                else:
                    ### (4b) else, check for co-ortholog ref
                    if checkCoorthologsRefOff == False:
                        aln_fa = '%s/blast_%s_%s_%s_%s_%s.fa' % (
                                    outpath, seqName, seed_id, search_taxon,
                                    hmm_hit, best_blast_hit)
                        with open(aln_fa, 'w') as aln_fa_out:
                            aln_fa_out.write(
                                '>%s\n%s\n>%s\n%s\n>%s\n%s' % (
                                    seed_id, refspec_seqs.fetch(seed_id),
                                    hmm_hit, search_seqs.fetch(hmm_hit),
                                    best_blast_hit, refspec_seqs.fetch(best_blast_hit)
                                )
                            )
                        fasta_fn.remove_dup(aln_fa)
                        aln_seq = align_fn.do_align(aligner, aln_fa)
                        output_fn.print_debug(
                            debug, 'Alignment for checking co-ortholog ref', aln_seq)
                        br_dist = align_fn.calc_Kimura_dist(aln_seq, best_blast_hit, seed_id, debug)
                        bh_dist = align_fn.calc_Kimura_dist(aln_seq, best_blast_hit, hmm_hit, debug)
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
                            ortho_candi[hmm_hit] = search_seqs.fetch(hmm_hit)
                            continue
    ### (5) check co-ortholog if more than 1 HMM hits are accepted
    if len(ortho_candi) == 0:
        output_fn.print_stdout(
            silentOff, 'WARNING: Reciprocity not fulfulled! No ortholog found!')
    else:
        best_ortho = list(ortho_candi.keys())[0]
        if not best_ortho == seed_id:
            ortho_final = fasta_fn.add_seq_to_dict(
                ortho_final, '%s|%s|%s|1' % (seqName, search_taxon, best_ortho),
                ortho_candi[best_ortho])
        if rep == False:
            if len(ortho_candi) > 1:
                aln_co_fa = '%s/coortho_%s_%s.fa' % (
                            outpath, seqName, search_taxon)
                with open(aln_co_fa, 'w') as aln_co_fa_out:
                    aln_co_fa_out.write(('>%s\n%s\n') %
                        (seed_id, refspec_seqs.fetch(seed_id)))
                    for cand in ortho_candi:
                        aln_co_fa_out.write(('>%s\n%s\n') %
                            (cand, ortho_candi[cand]))
                aln_co_seq = align_fn.do_align(aligner, aln_co_fa)
                output_fn.print_debug(
                    debug, 'Alignment for checking co-orthologs', aln_co_seq)
                if noCleanup == False:
                    os.remove(aln_co_fa)
                best_dist = align_fn.calc_Kimura_dist(
                                aln_co_seq, seed_id, best_ortho, debug)
                for cand in ortho_candi:
                    if not cand == best_ortho:
                        candi_dist = align_fn.calc_Kimura_dist(
                            aln_co_seq, best_ortho, cand, debug)
                        output_fn.print_debug(
                            debug,
                            'Check if distance bestHmm_vs_ref > '
                            + 'other_vs_bestHmm',
                            'd_best = %s; d_other = %s'
                            % (best_dist, candi_dist))
                        if candi_dist < best_dist:
                            if not cand == seed_id:
                                ortho_final = fasta_fn.add_seq_to_dict(
                                    ortho_final,
                                    '%s|%s|%s|0' \
                                        % (seqName, search_taxon, cand),
                                    ortho_candi[cand])
        output_fn.print_stdout(
            silentOff,
            '=> %s orthologs found: %s'
            % (len(ortho_final), list(ortho_final.keys())))
    return(ortho_final)


def run_hamstr(args):
    """ Perform ortholog search based on hamstr approach """

    (seqName, refspec, pathArgs, orthoArgs, otherArgs) = args
    (outpath, hmmpath, corepath, searchpath, annopath) = pathArgs
    (checkCoorthologsRefOff, rbh, rep, evalBlast, lowComplexityFilter,
                        evalHmmer, hitLimit, scoreCutoff, aligner) = orthoArgs
    (searchTaxa, cpus, debug, silentOff, noCleanup, force, append) = otherArgs

    hamstr_jobs = []
    ### get ref seqID
    core_fa = '%s/%s/%s.fa' % (hmmpath, seqName, seqName)
    seed_id = prepare_fn.get_seed_id_from_fa(core_fa, refspec)

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
                    evalHmmer, hitLimit, scoreCutoff,
                    evalBlast, lowComplexityFilter,
                    checkCoorthologsRefOff, rbh, rep,
                    aligner, cpus, debug, silentOff, noCleanup
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
                    evalHmmer, hitLimit, scoreCutoff,
                    evalBlast, lowComplexityFilter,
                    checkCoorthologsRefOff, rbh, rep,
                    aligner, cpus, debug, silentOff, noCleanup
                ])

    ### do ortholog search
    hamstr_out = {}
    if len(hamstr_jobs) > 0:
        output_fn.print_stdout(
                silentOff, 'Ortholog search for %s taxa...' % len(hamstr_jobs))
        if debug == True or silentOff == True or len(hamstr_jobs) == 1:
            for job in hamstr_jobs:
                tmp_out = hamstr(job)
                hamstr_out = {**hamstr_out, **tmp_out}
        else:
            pool = mp.Pool(cpus)
            for _ in tqdm(
                    pool.imap_unordered(hamstr, hamstr_jobs),
                    total=len(hamstr_jobs)):
                if len(_) > 0:
                    hamstr_out = {**hamstr_out, **_}

    ### Get seed seq
    refspec_fa = '%s/%s/%s.fa' % (corepath, refspec, refspec)
    refspec_seqs = fasta_fn.read_fasta(refspec_fa)
    seed_id_mod = '%s|%s|%s|1' % (seqName, refspec, seed_id)
    seed_seq = refspec_seqs.fetch(seed_id)

    ### return
    return({**{seed_id_mod:seed_seq}, **hamstr_out})
