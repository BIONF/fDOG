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
from pathlib import Path
from os.path import isfile, join
import argparse
import subprocess
import re
import shutil
import multiprocessing as mp
from tqdm import tqdm
from ete3 import NCBITaxa
from pkg_resources import get_distribution
import time

import fdog.libs.zzz as general_fn
import fdog.libs.preparation as prepare_fn
import fdog.libs.orthosearch as ortho_fn
import fdog.libs.corecompile as core_fn
import fdog.libs.fas as fas_fn
import fdog.libs.tree as tree_fn
import fdog.libs.output as output_fn



def get_sorted_files(directory):
    list = os.listdir(directory)
    pairs = []
    for file in list:
        location = os.path.join(directory, file)
        if isfile(location):
            size = os.path.getsize(location)
            pairs.append((size, file))
    pairs.sort(key=lambda s: s[0], reverse=True)
    return([x[1] for x in pairs])


def get_seed_name(seedFile):
    seqName = seedFile.rsplit('.', 1)[0]
    seqName = re.sub('[\|\.]', '_', seqName)
    return(seqName)


def create_core_jobs(args):
    (seed, core_options, other_options, inFol, outpath, silentOff) = args
    (coreArgs, orthoCoreArgs, otherCoreArgs) = core_options
    (refspec, reuseCore, forceCore, pathArgs, debug) = other_options
    (outpath, hmmpath, corepath, searchpath, annopath) = pathArgs
    seqFile = ('%s/%s' % (inFol, seed))
    seqName = get_seed_name(seed)
    if not os.path.exists('%s/core_orthologs/%s/hmm_dir/%s.hmm' % (outpath, seqName, seqName)) or forceCore == True:
        seed_id = prepare_fn.identify_seed_id(seqFile, refspec, corepath, debug, silentOff)
        if not seed_id == 'None':
            return([seqFile, seqName, refspec, seed_id,
                        reuseCore, forceCore, coreArgs, pathArgs, orthoCoreArgs,
                        otherCoreArgs, debug])
        else:
            print(f'WARNING: Cannot identify seed ID for {seqFile}!')


def compile_core(core_options, other_options, seeds, inFol, cpus, outpath, silentOff, jobName):
    core_compilation_jobs = []
    (coreArgs, orthoCoreArgs, otherCoreArgs) = core_options
    (cpus, debugCore, silentOff, noCleanup, force, append) = otherCoreArgs
    (refspec, reuseCore, forceCore, pathArgs, debug) = other_options
    (outpath, hmmpath, corepath, searchpath, annopath) = pathArgs
    pool = mp.Pool(cpus)
    begin = time.time()
    print('Preparing core compilation jobs...')
    core_job_file = '%s/%s_core_jobs.list' % (outpath, jobName)
    if os.path.exists(core_job_file) and os.stat(core_job_file).st_size > 0:
        print('... file contains jobs found (%s)' % core_job_file)
        core_compilation_jobs = general_fn.read_pyobj_file(core_job_file)
    else:
        prepare_jobs = []
        for seed in seeds:
            prepare_jobs.append([seed, core_options, other_options, inFol, outpath, silentOff])
        for _ in tqdm(pool.imap_unordered(create_core_jobs, prepare_jobs), total=len(prepare_jobs)):
            core_compilation_jobs.append(_)
        general_fn.save_pyobj(core_compilation_jobs, core_job_file)
    end = time.time()
    print('==> %s jobs will be run. Preparing finished in %s' % (len(core_compilation_jobs), '{:5.3f}s'.format(end - begin)))
    if len(core_compilation_jobs) > 0:
        core_runtime = []
        if debugCore == True or silentOff == True or len(core_compilation_jobs) == 1:
            for job in core_compilation_jobs:
                tmp_out = core_fn.run_compile_core(job)
                core_runtime.append(tmp_out)
        else:
            for _ in tqdm(pool.imap_unordered(core_fn.run_compile_core, core_compilation_jobs), total=len(core_compilation_jobs)):
                core_runtime.append(_)
        pool.close()
        pool.join()
        out = []
        for r in core_runtime:
            out.append('\t'.join(r))
        return(out)


def search_ortholog(options, seeds, inFol, outpath):
    (orthoArgs, otherArgs, pathArgs, refspec) = options
    (searchTaxa, cpus, debug, silentOff, noCleanup, force, append) = otherArgs
    ortho_runtime = []
    n = 1
    for seed in seeds:
        begin = time.time()
        seqFile = [inFol + '/' + seed]
        seqName = get_seed_name(seed)
        print('... %s (%s / %s)' % (seqName, n, len(seeds)))
        if not os.path.exists('%s/%s.extended.fa' % (outpath, seqName)) or force == True:
            hamstr_out = ortho_fn.run_hamstr([seqName, refspec, pathArgs, orthoArgs, otherArgs])
            output_fn.write_hamstr(hamstr_out, outpath, seqName, force, append)
            end = time.time()
            ortho_runtime.append('%s\t%s' % (seqName, '{:5.3f}s'.format(end - begin)))
        n += 1
    return(ortho_runtime)


def join_outputs(outpath, jobName, seeds, keep, silentOff):
    finalFa = '%s/%s.extended.fa' % (outpath, jobName)
    single_output_fol = '%s/%s_singleOutput' % (outpath, jobName)
    Path('%s/%s_singleOutput' % (outpath, jobName)).mkdir(parents=True, exist_ok=True)
    with open(finalFa,'wb') as wfd:
        for seed in seeds:
            seqName = get_seed_name(seed)
            resultFile = '%s/%s.extended.fa'  % (outpath, seqName)
            if silentOff == True:
                print(resultFile)
            if os.path.exists(resultFile):
                with open(resultFile,'rb') as fd:
                    shutil.copyfileobj(fd, wfd)
                if not os.path.exists('%s/%s.extended.fa' % (single_output_fol, seqName)):
                    shutil.move(resultFile, single_output_fol)
            if os.path.exists(outpath + '/' + seqName + '.fa'):
                os.remove(outpath + '/' + seqName + '.fa')
            if os.path.exists(os.getcwd() + '/' + seqName + '.fa'):
                os.remove(os.getcwd() + '/' + seqName + '.fa')
    if keep == True:
        try:
            print('Compressing single outputs...')
            shutil.make_archive(single_output_fol, 'gztar', single_output_fol)
        except:
            shutil.make_archive(single_output_fol, 'tar', single_output_fol)
    shutil.rmtree(single_output_fol)


def main():
    version = get_distribution('fdog').version
    parser = argparse.ArgumentParser(description='You are running fDOG version ' + str(version) + '.',
                                     epilog="For more information on certain options, please refer to the wiki pages "
                                            "on github: https://github.com/BIONF/fDOG/wiki")
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--seqFolder', help='Input folder containing the seed sequences (protein only) in fasta format',
                            action='store', default='', required=True)
    required.add_argument('--jobName', help='Job name. This will also be file name for the output',
                            action='store', default='', required=True)
    required.add_argument('--refspec', help='Reference taxon. It should be the species the seed sequence was derived from',
                            action='store', default='', required=True)

    optional_paths = parser.add_argument_group('Non-default directory options')
    optional_paths.add_argument('--outpath', help='Output directory', action='store', default='')
    optional_paths.add_argument('--hmmpath', help='Path for the core ortholog directory', action='store', default='')
    optional_paths.add_argument('--corepath', help='Path for the core taxa directory', action='store', default='')
    optional_paths.add_argument('--searchpath', help='Path for the search taxa directory', action='store', default='')
    optional_paths.add_argument('--annopath', help='Path for the pre-calculated feature annotion directory', action='store', default='')
    optional_paths.add_argument('--pathFile', help='Config file contains paths to data folder (in yaml format)', action='store', default='')

    core_options = parser.add_argument_group('Core compilation options')
    core_options.add_argument('--coreOnly', help='Compile only the core orthologs', action='store_true', default=False)
    core_options.add_argument('--reuseCore', help='Reuse existing core set of your sequence', action='store_true', default=False)
    core_options.add_argument('--minDist', help='Minimum systematic distance of primer taxa for the core set compilation. Default: genus',
                            choices=['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'superkingdom'],
                            action='store', default='genus')
    core_options.add_argument('--maxDist', help='Maximum systematic distance of primer taxa for the core set compilation. Default: kingdom',
                            choices=['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'superkingdom'],
                            action='store', default='kingdom')
    core_options.add_argument('--coreSize', help='Maximul number of orthologs in core set. Default: 6', action='store', default=6, type=int)
    core_options.add_argument('--coreTaxa', help='List of primer taxa that should exclusively be used for the core set compilation', action='store', default='')
    core_options.add_argument('--CorecheckCoorthologsOff', help='Turn off checking for co-ortholog of the reverse search during the core compilation',
                                action='store_true', default=False)
    core_options.add_argument('--coreRep', help='Obtain only the sequence being most similar to the corresponding sequence in the core set rather than all putative co-orthologs',
                                action='store_true', default=False)
    core_options.add_argument('--coreHitLimit', help='Number of hits of the initial pHMM based search that should be evaluated via a reverse search. Default: 3',
                                action='store', default=3, type=int)
    core_options.add_argument('--distDeviation', help='The deviation in score in percent (0 = 0 percent, 1 = 100 percent) allowed for two taxa to be considered similar. Default: 0.05',
                                action='store', default=0.05, type=float)
    core_options.add_argument('--alnStrategy', help='Specify the alignment strategy during core ortholog compilation. Default: local',
                                choices=['local', 'glocal', 'global'],
                                action='store', default='local')

    ortho_options = parser.add_argument_group('Ortholog search strategy options')
    ortho_options.add_argument('--searchTaxa', help='Specify file contains list of search taxa', action='store', default='')
    ortho_options.add_argument('--group', help='Allows to limit the search to a certain systematic group', action='store', default='')
    ortho_options.add_argument('--checkCoorthologsRefOff', help='Turn off checking for co-ortholog of the reverse search during the final ortholog search',
                                action='store_true', default=False)
    ortho_options.add_argument('--rbh', help='Requires a reciprocal best hit during the ortholog search to accept a new ortholog',
                                action='store_true', default=False)
    ortho_options.add_argument('--rep', help='Obtain only the sequence being most similar to the corresponding sequence in the core set rather than all putative co-orthologs',
                                action='store_true', default=False)
    ortho_options.add_argument('--lowComplexityFilter', help='Switch the low complexity filter for the blast search on. Default: False',
                                action='store_true', default=False)
    ortho_options.add_argument('--evalBlast', help='E-value cut-off for the Blast search. Default: 0.00001',
                                action='store', default=0.00005, type=float)
    ortho_options.add_argument('--evalHmmer', help='E-value cut-off for the HMM search. Default: 0.00001',
                                action='store', default=0.00005, type=float)
    ortho_options.add_argument('--hitLimit', help='number of hits of the initial pHMM based search that should be evaluated via a reverse search. Default: 10',
                                action='store', default=10, type=int)
    ortho_options.add_argument('--scoreCutoff', help='Define the percent range of the hmms core of the best hit up to which a candidate of the hmmsearch will be subjected for further evaluation. Default: 10',
                                action='store', default=10, type=int)

    fas_options = parser.add_argument_group('FAS options')
    fas_options.add_argument('--coreFilter',
                                help='Specifiy mode for filtering core orthologs by FAS score. In \'relaxed\' mode candidates with insufficient FAS score will be disadvantaged. In \'strict\' mode candidates with insufficient FAS score will be deleted from the candidates list. The option \'--minScore\' specifies the cut-off of the FAS score.',
                                choices=['relaxed', 'strict'], action='store', default='')
    fas_options.add_argument('--minScore', help='Specify the threshold for coreFilter. Default: 0.75', action='store', default=0.75, type=float)

    addtionalIO = parser.add_argument_group('Other I/O options')
    addtionalIO.add_argument('--append', help='Append the output to existing output files', action='store_true', default=False)
    addtionalIO.add_argument('--force', help='Overwrite existing ortholog search output files', action='store_true', default=False)
    addtionalIO.add_argument('--forceCore', help='Overwrite existing core set of your sequence', action='store_true', default=False)
    addtionalIO.add_argument('--noCleanup', help='Temporary output will NOT be deleted. Default: False', action='store_true', default=False)
    addtionalIO.add_argument('--keep', help='Keep output of individual seed sequence. Default: False', action='store_true', default=False)
    addtionalIO.add_argument('--debug', help='Set this flag to obtain more detailed information about the ortholog search progress', action='store_true', default=False)
    addtionalIO.add_argument('--debugCore', help='Set this flag to obtain more detailed information about the core compilation actions', action='store_true', default=False)
    addtionalIO.add_argument('--silentOff', help='Show more output to terminal', action='store_true', default=False)

    optional = parser.add_argument_group('Other options')
    optional.add_argument('--fasOff', help='Turn OFF FAS support', action='store_true', default=False)
    optional.add_argument('--aligner', help='Choose between mafft-linsi or muscle for the multiple sequence alignment. DEFAULT: muscle',
        choices=['mafft-linsi', 'muscle'], action='store', default='muscle')
    optional.add_argument('--cpus', help='Determine the number of threads to be run in parallel. Default: 4', action='store', default=4, type=int)

    ### get arguments
    args = parser.parse_args()

    # required arguments
    inFol = os.path.abspath(args.seqFolder)
    jobName = args.jobName
    refspec = args.refspec

    # path arguments
    outpath = os.path.abspath(args.outpath)
    hmmpath = args.hmmpath
    corepath = args.corepath
    searchpath = args.searchpath
    annopath = args.annopath
    pathFile = args.pathFile

    # core compilation arguments
    coreOnly = args.coreOnly
    reuseCore = args.reuseCore
    minDist = args.minDist
    maxDist = args.maxDist
    coreSize = args.coreSize
    coreTaxa = args.coreTaxa
    if not coreTaxa == '':
        if os.path.exists(os.path.abspath(coreTaxa)):
            coreTaxa = os.path.abspath(coreTaxa)
    CorecheckCoorthologsOff = args.CorecheckCoorthologsOff
    coreRep = args.coreRep
    coreHitLimit = args.coreHitLimit
    distDeviation = args.distDeviation
    alnStrategy = args.alnStrategy

    # ortholog search arguments
    searchTaxa = args.searchTaxa
    if not searchTaxa == '':
        if os.path.exists(os.path.abspath(searchTaxa)):
            searchTaxa = os.path.abspath(searchTaxa)
    group = args.group
    if not group == '' and not searchTaxa == '':
        print('WARNING: Both --group and --searchTaxa are specified. Search taxa will be obtained only from %s!' % searchTaxa)
        group = ''
    checkCoorthologsRefOff = args.checkCoorthologsRefOff
    rbh = args.rbh
    rep = args.rep
    lowComplexityFilter = args.lowComplexityFilter
    evalBlast = args.evalBlast
    evalHmmer = args.evalHmmer
    hitLimit = args.hitLimit
    scoreCutoff = args.scoreCutoff

    # fas arguments
    fasOff = args.fasOff
    coreFilter = args.coreFilter
    minScore = args.minScore

    # other I/O arguments
    append = args.append
    force = args.force
    forceCore = args.forceCore
    noCleanup = args.noCleanup
    keep = args.keep
    debug = args.debug
    debugCore = args.debugCore
    silentOff = args.silentOff

    # others
    aligner = args.aligner
    cpus = args.cpus
    if cpus > os.cpu_count():
        cpus = os.cpu_count()


    begin = time.time()
    ##### Check and group parameters
    (inFol, hmmpath, corepath, searchpath, annopath) = prepare_fn.check_input(
                    [inFol, refspec, outpath, hmmpath,
                    corepath, searchpath, annopath, pathFile])
    pathArgs = [outpath, hmmpath, corepath, searchpath, annopath]

    if not fasOff:
        check_fas = fas_fn.check_fas_executable()
        if check_fas == 0:
            sys.exit('ERROR: FAS is not executable! You still can use fDOG with --fasOff!')

    ### START
    Path(outpath).mkdir(parents=True, exist_ok=True)
    multiLog = open(outpath + '/' + jobName + '_log.txt', "w")

    print('PID %s - Jobname %s'% (str(os.getpid()), jobName))
    multiLog.write('PID %s - Jobname %s\n'% (str(os.getpid()), jobName))
    seeds = get_sorted_files(inFol)
    end = time.time()
    print('==> Sort seed files finished in ' + '{:5.3f}s'.format(end - begin))

    ##### DO CORE COMPILATION
    if reuseCore == False:
        print('Starting compiling core orthologs...')
        start = time.time()
        coreArgs = [minDist, maxDist, coreSize, coreTaxa, distDeviation,
                    alnStrategy, fasOff]
        orthoCoreArgs = [CorecheckCoorthologsOff, rbh, True, evalBlast,
                        lowComplexityFilter, evalHmmer, coreHitLimit,
                        scoreCutoff, aligner] # rep = True
        otherCoreArgs = [cpus, debugCore, silentOff, noCleanup, force, append]
        core_options = [coreArgs, orthoCoreArgs, otherCoreArgs]
        other_options = [refspec, reuseCore, forceCore, pathArgs, debug]
        core_runtime = compile_core(core_options, other_options, seeds, inFol, cpus, outpath, silentOff, jobName)
        end = time.time()
        multi_core_time = '{:5.3f}'.format(end-start)
        print('==> Core compilation finished in %ss\n' % multi_core_time)
        if len(core_runtime) > 1:
            multiLog.write('==> Core compilation finished in %ss\n%s\n' % (multi_core_time, '\n'.join(core_runtime)))
        else:
            multiLog.write('==> Core compilation finished in %ss\n' % multi_core_time)
    else:
        if not os.path.exists(hmmpath):
            sys.exit('--reuseCore was set, but no core orthologs found in %s! You could use --hmmpath to manually specify the core ortholog directory.' % outpath)


    ##### DO ORTHOLOG SEARCH USING HMM (HAMSTR)
    finalFa = '%s/%s.extended.fa' % (outpath, jobName)
    if not coreOnly:
        print('Searching orthologs...')
        start = time.time()
        if not os.path.exists(finalFa) or force == True:
            ### get list of search taxa
            if not group == '':
                ### Check valid taxonomy group
                ncbi = NCBITaxa()
                group_id = ncbi.get_name_translator([group])
                if len(group_id) == 0:
                    exit('ERROR: Taxon group "%s" invalid!' % group)
                ### create taxonomy tree from list of search taxa
                searchTaxa = []
                tax_ids = core_fn.get_core_taxa_ids(coreTaxa, corepath)

                for tax_id in tax_ids.keys():
                    check = tree_fn.check_taxon_group(group_id[group][0], tax_id, ncbi)
                    if check == True:
                        searchTaxa.append(tax_ids[tax_id])
                if debugCore:
                    print(searchTaxa)
                if len(searchTaxa) == 0:
                    exit('ERROR: No taxon found within %s taxonomy group!' % group)
                else:
                    searchTaxa = ','.join(searchTaxa)

            if len(searchTaxa) == '':
                searchTaxa = general_fn.read_dir(searchpath)
                searchTaxa = ','.join(searchTaxa)

            ### do ortholog search
            orthoArgs = [checkCoorthologsRefOff, rbh, rep, evalBlast,
                        lowComplexityFilter, evalHmmer, hitLimit, scoreCutoff, aligner]
            otherArgs = [searchTaxa, cpus, debug, silentOff, noCleanup, force, append]
            ortho_options = [orthoArgs, otherArgs, pathArgs, refspec]
            ortho_runtime = search_ortholog(ortho_options, seeds, inFol, outpath)
            end = time.time()
            multi_ortho_time = '{:5.3f}'.format(end-start)
            print('==> Ortholog search finished in %ss\n' % multi_ortho_time)
            multiLog.write('==> Ortholog search finished in %ss\n%s\n' % (multi_ortho_time, '\n'.join(ortho_runtime)))
            ### join output
            print('Joining single outputs...')
            start = time.time()
            join_outputs(outpath, jobName, seeds, keep, silentOff)
            end = time.time()
            print('==> Joining outputs finished in %ss\n' % '{:5.3f}'.format(end-start))

        ##### DO FINAL FAS CALCULATION
        if not fasOff:
            try:
                fasVersion = subprocess.run(['fas.run --version'], shell = True, capture_output = True, check = True)
            except:
                sys.exit('Problem with FAS! Please check https://github.com/BIONF/FAS or turn it off if not needed!')
            if os.path.exists(finalFa):
                start = time.time()
                fas_fn.calc_fas_multi(finalFa, outpath, annopath, cpus)
                end = time.time()
                print('==> FAS calculation finished in ' + '{:5.3f}s'.format(end - start))
                multiLog.write('==> FAS calculation finished in ' + '{:5.3f}s'.format(end - start))
        else:
            output_fn.hamstr_2_profile(finalFa)

        end = time.time()
        print('==> fdogs.run finished in ' + '{:5.3f}s'.format(end - begin))

if __name__ == '__main__':
    main()
