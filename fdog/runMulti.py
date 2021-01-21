# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2020 Vinh Tran
#
#  This script is used to run fdog with multiple seed sequences.
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
from os import listdir
from os.path import isfile, join
import time
import argparse
import subprocess
from pathlib import Path
import multiprocessing as mp
import re
from tqdm import tqdm
import fdog.runSingle as fdogFn
import shutil

def getSortedFiles(directory):
    list = os.listdir(directory)
    pairs = []
    for file in list:
        location = os.path.join(directory, file)
        if isfile(location):
            size = os.path.getsize(location)
            pairs.append((size, file))
    pairs.sort(key=lambda s: s[0], reverse=True)
    return([x[1] for x in pairs])

def prepare(args, step):
    (seqFile, seqName, fdogPath, refspec, minDist, maxDist, coreOrth,
    append, force, cleanup, group, blast, db,
    outpath, hmmpath, blastpath, searchpath, weightpath,
    coreOnly, reuseCore, coreTaxa, coreStrict, CorecheckCoorthologsRef, coreRep, coreHitLimit, distDeviation,
    fasoff, countercheck, coreFilter, minScore,
    strict, checkCoorthologsRef, rbh, rep, ignoreDistance, lowComplexityFilterOff, evalBlast, evalHmmer, evalRelaxfac, hitLimit, autoLimit, scoreThreshold, scoreCutoff, aligner, local, glocal, searchTaxa,
    cpu, hyperthread, debug, silent) = args

    mute = False
    if step == 'core':
        coreOnly = True
        silent = True
        mute = True
    else:
        reuseCore = True
        fasoff = True
        if silent == True:
            mute = True
    ### check input arguments
    seqFile, hmmpath, blastpath, searchpath, weightpath = fdogFn.checkInput([fdogPath, seqFile, refspec, outpath, hmmpath, blastpath, searchpath, weightpath])
    # group arguments
    basicArgs = [fdogPath, seqFile, seqName, refspec, minDist, maxDist, coreOrth]
    ioArgs = [append, force, cleanup, group, blast, db]
    pathArgs = [outpath, hmmpath, blastpath, searchpath, weightpath]
    coreArgs = [coreOnly, reuseCore, coreTaxa, coreStrict, CorecheckCoorthologsRef, coreRep, coreHitLimit, distDeviation]
    fasArgs = [fasoff, countercheck, coreFilter, minScore]
    orthoArgs = [strict, checkCoorthologsRef, rbh, rep, ignoreDistance, lowComplexityFilterOff, evalBlast, evalHmmer, evalRelaxfac, hitLimit, autoLimit, scoreThreshold, scoreCutoff, aligner, local, glocal, searchTaxa]
    otherArgs = [cpu, hyperthread, debug, True]
    return(basicArgs, ioArgs, pathArgs, coreArgs, orthoArgs, fasArgs, otherArgs, mute)

def getSeedName(seedFile):
    seqName = seedFile.split('.')[0]
    seqName = re.sub('[\|\.]', '_', seqName)
    return(seqName)

def getIndividualRuntime(step, outpath, seeds):
    logFile = outpath + '/runtime_core.txt'
    searchTerm = 'Core set compilation finished in'
    if step == 'ortho':
        logFile = outpath + '/runtime_ortho.txt'
        searchTerm = 'Ortholog search completed in'
    log = open(logFile, "w")
    for seed in seeds:
        seqName = getSeedName(seed)
        logFile = outpath + '/' + seqName + '/fdog.log'
        if os.path.exists(logFile):
            with open(logFile, 'r') as f:
                for line in f:
                    if searchTerm in line:
                        runtime = line.split()[-2]
                        log.write('%s\t%s\n' % (seqName, runtime))
        else:
            missing = open(outpath + '/missing.txt', 'a+')
            missing.write(step + '\t' + seqName + '\n')
    log.close()

def compileCore(options, seeds, inFol, cpu, outpath):
    print('Starting compiling core orthologs...')
    start = time.time()
    coreCompilationJobs = []
    for seed in seeds:
        seqFile = [inFol + '/' + seed]
        seqName = getSeedName(seed)
        (basicArgs, ioArgs, pathArgs, coreArgs, orthoArgs, fasArgs, otherArgs, mute) = prepare(seqFile + [seqName] + options, 'core')
        coreCompilationJobs.append([basicArgs, ioArgs, pathArgs, coreArgs, orthoArgs, fasArgs, otherArgs, mute])
    pool = mp.Pool(cpu)
    coreOut = []
    for _ in tqdm(pool.imap_unordered(fdogFn.runSingle, coreCompilationJobs), total=len(coreCompilationJobs)):
        coreOut.append(_)
    pool.close()
    pool.join()
    end = time.time()
    # read logs file to get runtime for individual seeds
    getIndividualRuntime('core', outpath, seeds)
    multiCoreTime = '{:5.3f}'.format(end-start)
    print('==> Core compiling finished in %s sec' % multiCoreTime) #'{:5.3f}s'.format(end-start))
    return(multiCoreTime)

def searchOrtho(options, seeds, inFol, cpu, outpath):
    print('Searching orthologs for...')
    start = time.time()
    coreCompilationJobs = []
    for seed in seeds:
        seqFile = [inFol + '/' + seed]
        seqName = getSeedName(seed)
        (basicArgs, ioArgs, pathArgs, coreArgs, orthoArgs, fasArgs, otherArgs, mute) = prepare(seqFile + [seqName] + options, 'ortholog')
        if mute == True:
            print(seed)
        else:
            print('\n##### ' + seed)
        fdogFn.runSingle([basicArgs, ioArgs, pathArgs, coreArgs, orthoArgs, fasArgs, otherArgs, mute])
    end = time.time()
    # read logs file to get runtime for individual seeds
    getIndividualRuntime('ortho', outpath, seeds)
    multiOrthoTime = '{:5.3f}'.format(end-start)
    print('==> Ortholog search finished in %s sec' % multiOrthoTime)
    return(multiOrthoTime)

def joinOutputs(outpath, jobName, seeds, keep):
    print('Joining single outputs...')
    finalFa = '%s/%s.extended.fa' % (outpath, jobName)
    Path(outpath+'/singleOutput').mkdir(parents=True, exist_ok=True)
    with open(finalFa,'wb') as wfd:
        for seed in seeds:
            seqName = getSeedName(seed)
            resultFile = '%s/%s/%s.extended.fa'  % (outpath, seqName, seqName)
            if os.path.exists(resultFile):
                with open(resultFile,'rb') as fd:
                    shutil.copyfileobj(fd, wfd)
                shutil.move(outpath + '/' + seqName, outpath + '/singleOutput')
            else:
                Path(outpath+'/missingOutput').mkdir(parents=True, exist_ok=True)
                shutil.move(outpath + '/' + seqName, outpath + '/missingOutput')
                os.remove(outpath + '/' + seqName + '.fa')
    if keep == True:
        try:
            print('Compressing single outputs...')
            shutil.make_archive(outpath + '/' + jobName + '_singleOutput', 'gztar', outpath+'/singleOutput')
        except:
            shutil.make_archive(outpath + '/' + jobName + '_singleOutput', 'tar', outpath+'/singleOutput')
    shutil.rmtree(outpath + '/singleOutput')
    return(finalFa)

def calcFAS (outpath, extendedFa, weightpath, cpu):
    print('Starting calculating FAS scores...')
    start = time.time()
    fasCmd = 'fdogFAS -i %s -w %s --cores %s' % (extendedFa, weightpath, cpu)
    try:
        subprocess.call([fasCmd], shell = True)
        end = time.time()
        if os.path.exists(outpath + '/tmp'):
            shutil.rmtree(outpath + '/tmp')
        fasTime = '{:5.3f}s'.format(end-start)
        print('==> FAS calculation finished in %s sec' % fasTime)
        return(fasTime)
    except:
        sys.exit('Problem running\n%s' % (fasCmd))

def main():
    version = '0.0.16'
    parser = argparse.ArgumentParser(description='You are running fdogs.run version ' + str(version) + '.')
    parser.add_argument('--version', action='version', version=str(version))
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--input', help='Input folder containing the seed sequences (protein only) in fasta format',
                            action='store', default='', required=True)
    required.add_argument('--jobName', help='Job name. This will also be file name for the output',
                            action='store', default='', required=True)
    required.add_argument('--refspec', help='Reference taxon. It should be the species the seed sequence was derived from',
                            action='store', default='', required=True)

    optional_paths = parser.add_argument_group('Non-default directory options')
    optional_paths.add_argument('--outpath', help='Output directory', action='store', default='')
    optional_paths.add_argument('--hmmpath', help='Path for the core ortholog directory', action='store', default='')
    optional_paths.add_argument('--blastpath', help='Path for the blastDB directory', action='store', default='')
    optional_paths.add_argument('--searchpath', help='Path for the search taxa directory', action='store', default='')
    optional_paths.add_argument('--weightpath', help='Path for the pre-calculated feature annotion directory', action='store', default='')

    addtionalIO = parser.add_argument_group('Other I/O options')
    addtionalIO.add_argument('--append', help='Append the output to existing output files', action='store_true', default=False)
    addtionalIO.add_argument('--force', help='Overwrite existing output files', action='store_true', default=False)
    addtionalIO.add_argument('--cleanup', help='Temporary output will be deleted. Default: True', action='store_true', default=True)
    addtionalIO.add_argument('--keep', help='Keep output of individual seed sequence. Default: False', action='store_true', default=False)
    addtionalIO.add_argument('--group', help='Allows to limit the search to a certain systematic group', action='store', default='')
    addtionalIO.add_argument('--blast', help='Determine sequence id and refspec automatically. Note, the chosen sequence id and reference species does not necessarily reflect the species the sequence was derived from.',
                                action='store_true', default=False)
    addtionalIO.add_argument('--db', help='Run fdog in database mode. Requires a mySql database. Only for internal use.', action='store_true', default=False)

    core_options = parser.add_argument_group('Core compilation options')
    core_options.add_argument('--coreOnly', help='Compile only the core orthologs', action='store_true', default=False)
    core_options.add_argument('--reuseCore', help='Reuse existing core set of your sequence', action='store_true', default=False)
    core_options.add_argument('--minDist', help='Minimum systematic distance of primer taxa for the core set compilation. Default: genus',
                            choices=['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom'],
                            action='store', default='genus')
    core_options.add_argument('--maxDist', help='Maximum systematic distance of primer taxa for the core set compilation. Default: kingdom',
                            choices=['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom'],
                            action='store', default='kingdom')
    core_options.add_argument('--coreOrth', help='Number of orthologs added to the core set. Default: 5', action='store', default=5, type=int)
    core_options.add_argument('--coreTaxa', help='List of primer taxa that should exclusively be used for the core set compilation', action='store', default='')
    core_options.add_argument('--coreStrict', help='An ortholog is only then accepted when the reciprocity is fulfilled for each sequence in the core set',
                                action='store_true', default=False)
    core_options.add_argument('--CorecheckCoorthologsRef', help='During the core compilation, an ortholog also be accepted when its best hit in the reverse search is not the core ortholog itself, but a co-ortholog of it',
                                action='store_true', default=False)
    core_options.add_argument('--coreRep', help='Obtain only the sequence being most similar to the corresponding sequence in the core set rather than all putative co-orthologs',
                                action='store_true', default=False)
    core_options.add_argument('--coreHitLimit', help='Number of hits of the initial pHMM based search that should be evaluated via a reverse search. Default: 3',
                                action='store', default=3, type=int)
    core_options.add_argument('--distDeviation', help='The deviation in score in percent (0 = 0 percent, 1 = 100 percent) allowed for two taxa to be considered similar. Default: 0.05',
                                action='store', default=0.05, type=float)

    ortho_options = parser.add_argument_group('Search strategy options')
    ortho_options.add_argument('--strict', help='An ortholog is only then accepted when the reciprocity is fulfilled for each sequence in the core set',
                                action='store_true', default=False)
    ortho_options.add_argument('--checkCoorthologsRef', help='During the final ortholog search, accept an ortholog also when its best hit in the reverse search is not the core ortholog itself, but a co-ortholog of it',
                                action='store_true', default=False)
    ortho_options.add_argument('--rbh', help='Requires a reciprocal best hit during the ortholog search to accept a new ortholog',
                                action='store_true', default=False)
    ortho_options.add_argument('--rep', help='Obtain only the sequence being most similar to the corresponding sequence in the core set rather than all putative co-orthologs',
                                action='store_true', default=False)
    ortho_options.add_argument('--ignoreDistance', help='Ignore the distance between Taxa and to choose orthologs only based on score',
                                action='store_true', default=False)
    ortho_options.add_argument('--lowComplexityFilterOff', help='Switch on or off the low complexity filter for the blast search. Default: False',
                                action='store_true', default=False)
    ortho_options.add_argument('--evalBlast', help='E-value cut-off for the Blast search. Default: 0.00005',
                                action='store', default=0.00005, type=float)
    ortho_options.add_argument('--evalHmmer', help='E-value cut-off for the HMM search. Default: 0.00005',
                                action='store', default=0.00005, type=float)
    ortho_options.add_argument('--evalRelaxfac', help='The factor to relax the e-value cut-off (Blast search and HMM search). Default: 10',
                                action='store', default=10, type=int)
    ortho_options.add_argument('--hitLimit', help='number of hits of the initial pHMM based search that should be evaluated via a reverse search. Default: 10',
                                action='store', default=10, type=int)
    ortho_options.add_argument('--autoLimit', help='Invoke a lagPhase analysis on the score distribution from the hmmer search. This will determine automatically a hit limit for each query. Note, it will be effective for both the core compilation and the final ortholog search',
                                action='store_true', default=False)
    ortho_options.add_argument('--scoreThreshold', help='Instead of setting an automatic hit limit, you can specify with this flag that only candidates with an hmm score no less than x percent of the hmm score of the best hit are further evaluated. Default: x = 10. You can change this cutoff with the option -scoreCutoff. Note, it will be effective for both the core compilation and the final ortholog search',
                                action='store_true', default=False)
    ortho_options.add_argument('--scoreCutoff', help='In combination with -scoreThreshold you can define the percent range of the hmms core of the best hit up to which a candidate of the hmmsearch will be subjected for further evaluation. Default: 10',
                                action='store', default=10, type=int)
    ortho_options.add_argument('--aligner', help='Choose between mafft-linsi or muscle for the multiple sequence alignment. DEFAULT: muscle',
                                choices=['mafft-linsi', 'muscle'], action='store', default='muscle')
    ortho_options.add_argument('--local', help='Specify the alignment strategy during core ortholog compilation. Default: True',
                                action='store_true', default=True)
    ortho_options.add_argument('--glocal', help='Specify the alignment strategy during core ortholog compilation. Default: False',
                                action='store_true', default=False)
    ortho_options.add_argument('--searchTaxa', help='Specify list of search taxa', action='store', default='')

    fas_options = parser.add_argument_group('FAS options')
    fas_options.add_argument('--fasoff', help='Turn OFF FAS support', action='store_true', default=False)
    fas_options.add_argument('--countercheck', help='The FAS score will be computed in two ways', action='store_true', default=True)
    fas_options.add_argument('--coreFilter',
                                help='Specifiy mode for filtering core orthologs by FAS score. In \'relaxed\' mode candidates with insufficient FAS score will be disadvantaged. In \'strict\' mode candidates with insufficient FAS score will be deleted from the candidates list. The option \'--minScore\' specifies the cut-off of the FAS score.',
                                choices=['relaxed', 'strict'], action='store', default='')
    fas_options.add_argument('--minScore', help='Specify the threshold for coreFilter. Default: 0.75', action='store', default=0.75, type=float)

    optional = parser.add_argument_group('Other options')
    optional.add_argument('--cpu', help='Determine the number of threads to be run in parallel. Default: 4', action='store', default=4, type=int)
    optional.add_argument('--hyperthread', help='Set this flag to use hyper threading. Default: False', action='store_true', default=False)
    optional.add_argument('--debug', help='Set this flag to obtain more detailed information about the programs actions', action='store_true', default=False)
    optional.add_argument('--silentOff', help='Show more output to terminal', action='store_true', default=False)

    ### get arguments
    args = parser.parse_args()

    # required arguments
    inFol = os.path.abspath(args.input)
    jobName = args.jobName
    refspec = args.refspec

    minDist = args.minDist
    maxDist = args.maxDist
    coreOrth = args.coreOrth

    # path arguments
    outpath = os.path.abspath(args.outpath)
    hmmpath = args.hmmpath
    blastpath = args.blastpath
    searchpath = args.searchpath
    weightpath = args.weightpath

    # other I/O arguments
    append = args.append
    force = args.force
    cleanup = args.cleanup
    keep = args.keep
    group = args.group
    blast = args.blast
    db = args.db

    # core compilation arguments
    coreOnly = args.coreOnly
    reuseCore = args.reuseCore
    coreTaxa = args.coreTaxa
    coreStrict = args.coreStrict
    CorecheckCoorthologsRef = args.CorecheckCoorthologsRef
    coreRep = args.coreRep
    coreHitLimit = args.coreHitLimit
    distDeviation = args.distDeviation

    # ortholog search arguments
    strict = args.strict
    checkCoorthologsRef = args.checkCoorthologsRef
    rbh = args.rbh
    rep = args.rep
    ignoreDistance = args.ignoreDistance
    lowComplexityFilterOff = args.lowComplexityFilterOff
    evalBlast = args.evalBlast
    evalHmmer = args.evalHmmer
    evalRelaxfac = args.evalRelaxfac
    hitLimit = args.hitLimit
    autoLimit = args.autoLimit
    scoreThreshold = args.scoreThreshold
    scoreCutoff = args.scoreCutoff
    aligner = args.aligner
    local = args.local
    glocal = args.glocal
    searchTaxa = args.searchTaxa

    # fas arguments
    fasoff = args.fasoff
    countercheck = args.countercheck
    coreFilter = args.coreFilter
    minScore = args.minScore

    # others
    cpu = args.cpu
    hyperthread = args.hyperthread
    debug = args.debug
    silentOff = args.silentOff
    if silentOff == True:
        silent = False
    else:
        silent = True

    ### get fdog and data path
    fdogPath = os.path.realpath(__file__).replace('/runMulti.py','')
    pathconfigFile = fdogPath + '/bin/pathconfig.txt'
    if not os.path.exists(pathconfigFile):
        sys.exit('No pathconfig.txt found. Please run fdog.setup (https://github.com/BIONF/fDOG/wiki/Installation#setup-fdog).')
    with open(pathconfigFile) as f:
        dataPath = f.readline().strip()
    if hmmpath == '':
        hmmpath = dataPath + '/core_orthologs'
    if blastpath == '':
        blastpath = dataPath + '/blast_dir'
    if searchpath == '':
        searchpath = dataPath + '/genome_dir'
    if weightpath == '':
        weightpath = dataPath + '/weight_dir'

    ### join options
    options = [fdogPath, refspec, minDist, maxDist, coreOrth,
                append, force, cleanup, group, blast, db,
                outpath, hmmpath, blastpath, searchpath, weightpath,
                coreOnly, reuseCore, coreTaxa, coreStrict, CorecheckCoorthologsRef, coreRep, coreHitLimit, distDeviation,
                fasoff, countercheck, coreFilter, minScore,
                strict, checkCoorthologsRef, rbh, rep, ignoreDistance, lowComplexityFilterOff, evalBlast, evalHmmer, evalRelaxfac, hitLimit, autoLimit, scoreThreshold, scoreCutoff, aligner, local, glocal, searchTaxa,
                cpu, hyperthread, debug, silent]

    ### START
    multiLog = open(outpath + '/' + jobName + '_log.txt', "w")
    fdogStart = time.time()
    seeds = getSortedFiles(inFol)
    print('PID ' + str(os.getpid()))
    multiLog.write('PID ' + str(os.getpid()) + '\n')

    ### run core compilation
    if reuseCore == False:
        multiCoreTime = compileCore(options, seeds, inFol, cpu, outpath)
        multiLog.write('==> Core compilation finished in %s sec\n' % multiCoreTime)

    ### do ortholog search
    if coreOnly == False:
        ### create list of search taxa
        searchTaxa = ''
        searchGroup = 'all'
        if not group == '':
            print('Creating list for search taxa...')
            searchTaxa = '%s/searchTaxa.txt' % (outpath)
            searchGroup = group
            cmd = 'perl %s/bin/getSearchTaxa.pl -i %s -b %s -h %s -r %s -n %s -t %s/taxonomy -o %s' % (fdogPath, searchpath, evalBlast, evalHmmer, evalRelaxfac, searchGroup, fdogPath, searchTaxa)
            try:
                subprocess.call([cmd], shell = True)
            except:
                sys.exit('Problem running\n%s' % (cmd))
        ### run ortholog search
        multiOrthoTime = searchOrtho(options, seeds, inFol, cpu, outpath)
        multiLog.write('==> Ortholog search finished in %s sec\n' % multiOrthoTime)
        ### join output
        finalFa = joinOutputs(outpath, jobName, seeds, keep)
        ### calculate FAS scores
        if fasoff == False:
            fasTime = calcFAS(outpath, finalFa, weightpath, cpu)
            multiLog.write('==> FAS calculation finished in %s sec\n' % fasTime)

    fdogEnd = time.time()
    print('==> fdogs.run finished in ' + '{:5.3f}s'.format(fdogEnd-fdogStart))
    multiLog.write('==> fdogs.run finished in ' + '{:5.3f}s'.format(fdogEnd-fdogStart))
    multiLog.close()

if __name__ == '__main__':
    main()
