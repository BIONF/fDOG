# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2020 Vinh Tran
#
#  This script is used to run fdog for one seed sequence.
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
import argparse
import subprocess
from pathlib import Path
import yaml

def checkFileExist(file):
    if not os.path.exists(os.path.abspath(file)):
        sys.exit('%s not found' % file)

def load_config(config_file):
    with open(config_file, 'r') as stream:
        try:
            return yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

def checkInput(args):
    (fdogPath, seqFile, refspec, outpath, hmmpath, blastpath, searchpath, weightpath) = args
    # create output directory
    Path(outpath).mkdir(parents=True, exist_ok=True)
    Path(hmmpath).mkdir(parents=True, exist_ok=True)
    # check path existing
    for path in [hmmpath, blastpath, searchpath, weightpath]:
        checkFileExist(path)
    # check for seqFile
    if not os.path.exists(os.path.abspath(seqFile)):
        if not os.path.exists(fdogPath + '/data/' + seqFile):
            sys.exit('%s not found in %s or %s' % (seqFile, os.getcwd(), fdogPath + '/data/'))
        else:
            seqFile = fdogPath + '/data/' + seqFile
    else:
        seqFile = os.path.abspath(seqFile)
    # check refspec
    if not os.path.exists(os.path.abspath(blastpath+'/'+refspec)):
        exit('Reference taxon %s not found in %s' % (refspec, blastpath))
    return (seqFile, hmmpath, blastpath, searchpath, weightpath)

def getfdogInfo(fdogPath, infoType):
    if os.path.exists(fdogPath + '/bin/oneSeq.pl'):
        cmd = subprocess.Popen([fdogPath + '/bin/oneSeq.pl', infoType], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        msg, err = cmd.communicate()
        print(msg.decode('UTF-8').strip())
        print(err.decode('UTF-8').strip())
        exit()
    else:
        exit('%s not found' % (fdogPath + '/bin/oneSeq.pl'))

def runSingle(args):
    (basicArgs, ioArgs, pathArgs, coreArgs, orthoArgs, fasArgs, otherArgs, mute) = args
    # basic command
    (fdogPath, seqFile, seqName, refspec, minDist, maxDist, coreOrth) = basicArgs
    cmd = 'perl %s/bin/oneSeq.pl -seqFile=%s -seqName=%s -refspec=%s' % (fdogPath, seqFile, seqName, refspec)
    # add paths
    (outpath, hmmpath, blastpath, searchpath, weightpath) = pathArgs
    cmd = cmd + ' -outpath=%s -hmmpath=%s -blastpath=%s -searchpath=%s -weightpath=%s' % (outpath, hmmpath, blastpath, searchpath, weightpath)
    # add other I/O options
    (append, force, noCleanup, group, blast, db) = ioArgs
    if append == True:
        cmd = cmd + ' -append'
    if force == True:
        cmd = cmd + ' -force'
    if noCleanup == False:
        cmd = cmd + ' -cleanup'
    if blast == True:
        cmd = cmd + ' -blast'
    if db == True:
        cmd = cmd + ' -db'
    if not group == '':
        cmd = cmd + ' -group=%s' % group
    # add core compilation options
    (coreOnly, reuseCore, coreTaxa, coreStrict, CorecheckCoorthologsRef, coreRep, coreHitLimit, distDeviation) = coreArgs
    if coreOnly == True:
        cmd = cmd + ' -coreOnly'
    if reuseCore == True:
        cmd = cmd + ' -reuseCore'
    else:
        cmd = cmd + ' -minDist=%s -maxDist=%s -coreOrth=%s' % (minDist, maxDist, coreOrth)
    if not coreTaxa == '':
        cmd = cmd + ' -coreTaxa=%s' % coreTaxa
    if coreStrict == True:
        cmd = cmd + ' -coreStrict'
    if CorecheckCoorthologsRef == True:
        cmd = cmd + ' -CorecheckCoorthologsRef'
    if coreRep == True:
        cmd = cmd + ' -coreRep'
    if not coreHitLimit == 3:
        cmd = cmd + ' -coreHitLimit=%s' % coreHitLimit
    if not distDeviation == 0.05:
        cmd = cmd + ' -distDeviation=%s' % distDeviation
    # add ortholo search options
    (strict, checkCoorthologsRef, rbh, rep, ignoreDistance, lowComplexityFilterOff, evalBlast, evalHmmer, evalRelaxfac, hitLimit, autoLimit, scoreThreshold, scoreCutoff, aligner, local, glocal, searchTaxa) = orthoArgs
    if strict == True:
        cmd = cmd + ' -strict'
    if checkCoorthologsRef == True:
        cmd = cmd + ' -checkCoorthologsRef'
    if rbh == True:
        cmd = cmd + ' -rbh'
    if rep == True:
        cmd = cmd + ' -rep'
    if ignoreDistance == True:
        cmd = cmd + ' -ignoreDistance'
    if lowComplexityFilterOff == True:
        cmd = cmd + ' -filter=F'
    if not evalBlast == 0.00005:
        cmd = cmd + ' -evalBlast=%s' % evalBlast
    if not evalHmmer == 0.00005:
        cmd = cmd + ' -evalHmmer=%s' % evalHmmer
    if not evalRelaxfac == 10:
        cmd = cmd + ' -evalRelaxfac=%s' % evalRelaxfac
    if not hitLimit == 10:
        cmd = cmd + ' -hitLimit=%s' % hitLimit
    if autoLimit == True:
        cmd = cmd + ' -autoLimit'
    if scoreThreshold:
        cmd = cmd + ' -scoreThreshold'
    if not scoreCutoff == 10:
        cmd = cmd + ' -scoreCutoff=%s' % scoreCutoff
    if not aligner == 'muscle':
        cmd = cmd + ' -aligner=%s' % aligner
    if glocal == True:
        cmd = cmd + ' -glocal'
    if not searchTaxa == '':
        checkFileExist(searchTaxa)
        searchTaxa = os.path.abspath(searchTaxa)
        cmd = cmd + ' -searchTaxa=%s' % searchTaxa
    # add fas options
    (fasoff, countercheck, coreFilter, minScore) = fasArgs
    if fasoff == True:
        cmd = cmd + ' -fasoff'
    else:
        if countercheck == True:
            cmd = cmd + ' -countercheck'
        if not coreFilter == '':
            if minScore > 0:
                cmd = cmd + ' -coreFilter=%s -minScore=%s' % (coreFilter, minScore)
    # add other options
    (cpu, hyperthread, debug, silent) = otherArgs
    cmd = cmd + ' -cpu=%s' % cpu
    if hyperthread == True:
        cmd = cmd + ' -hyperthread'
    if debug == True:
        cmd = cmd + ' -debug'
    if silent == True:
        cmd = cmd + ' -silent'
    # print(cmd)
    if mute == True:
        cmd = cmd + ' > /dev/null 2>&1'
    try:
        subprocess.call([cmd], shell = True)
    except:
        sys.exit('Problem running\n%s' % (cmd))

def main():
    version = '0.0.23'
    parser = argparse.ArgumentParser(description='You are running fdog.run version ' + str(version) + '.')
    parser.add_argument('--version', action='version', version=str(version))
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--seqFile', help='Input file containing the seed sequence (protein only) in fasta format',
                            action='store', default='', required=True)
    required.add_argument('--seqName', help='Job name. This will also be file name for the output',
                            action='store', default='', required=True)
    required.add_argument('--refspec', help='Reference taxon. It should be the species the seed sequence was derived from',
                            action='store', default='', required=True)

    optional_paths = parser.add_argument_group('Non-default directory options')
    optional_paths.add_argument('--outpath', help='Output directory', action='store', default='')
    optional_paths.add_argument('--hmmpath', help='Path for the core ortholog directory', action='store', default='')
    optional_paths.add_argument('--blastpath', help='Path for the blastDB directory', action='store', default='')
    optional_paths.add_argument('--searchpath', help='Path for the search taxa directory', action='store', default='')
    optional_paths.add_argument('--weightpath', help='Path for the pre-calculated feature annotion directory', action='store', default='')
    optional_paths.add_argument('--pathFile', help='Config file contains paths to data folder (in yaml format)', action='store', default='')

    addtionalIO = parser.add_argument_group('Other I/O options')
    addtionalIO.add_argument('--append', help='Append the output to existing output files', action='store_true', default=False)
    addtionalIO.add_argument('--force', help='Overwrite existing output files', action='store_true', default=False)
    addtionalIO.add_argument('--noCleanup', help='Temporary output will NOT be deleted. Default: False', action='store_true', default=False)
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

    ortho_options = parser.add_argument_group('Ortholog search strategy options')
    ortho_options.add_argument('--searchTaxa', help='Specify file contains list of search taxa', action='store', default='')
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
    seqFile = args.seqFile
    seqName = args.seqName
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
    pathFile = args.pathFile

    # other I/O arguments
    append = args.append
    force = args.force
    noCleanup = args.noCleanup
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
    dataPath = ''
    fdogPath = os.path.realpath(__file__).replace('/runSingle.py','')
    pathconfigFile = fdogPath + '/bin/pathconfig.txt'
    if not os.path.exists(pathconfigFile):
        sys.exit('No pathconfig.txt found. Please run fdog.setup (https://github.com/BIONF/fDOG/wiki/Installation#setup-fdog).')
    if pathFile == '':
        with open(pathconfigFile) as f:
            dataPath = f.readline().strip()
    else:
        cfg = load_config(pathFile)
        try:
            dataPath = cfg['dataPath']
        except:
            dataPath = 'config'

    if hmmpath == '':
        hmmpath = outpath + '/core_orthologs'
    #     hmmpath = dataPath + '/core_orthologs'
    #     if dataPath == 'config':
    #         try:
    #             hmmpath = cfg['hmmpath']
    #         except:
    #             sys.exit('hmmpath not found in %s' % pathFile)

    if blastpath == '':
        blastpath = dataPath + '/blast_dir'
        if dataPath == 'config':
            try:
                blastpath = cfg['blastpath']
            except:
                sys.exit('blastpath not found in %s' % pathFile)
    if searchpath == '':
        searchpath = dataPath + '/genome_dir'
        if dataPath == 'config':
            try:
                searchpath = cfg['searchpath']
            except:
                sys.exit('searchpath not found in %s' % pathFile)
    if weightpath == '':
        weightpath = dataPath + '/weight_dir'
        if dataPath == 'config':
            try:
                weightpath = cfg['weightpath']
            except:
                sys.exit('weightpath not found in %s' % pathFile)

    ### check input arguments
    seqFile, hmmpath, blastpath, searchpath, weightpath = checkInput([fdogPath, seqFile, refspec, outpath, hmmpath, blastpath, searchpath, weightpath])
    # group arguments
    basicArgs = [fdogPath, seqFile, seqName, refspec, minDist, maxDist, coreOrth]
    ioArgs = [append, force, noCleanup, group, blast, db]
    pathArgs = [outpath, hmmpath, blastpath, searchpath, weightpath]
    coreArgs = [coreOnly, reuseCore, coreTaxa, coreStrict, CorecheckCoorthologsRef, coreRep, coreHitLimit, distDeviation]
    fasArgs = [fasoff, countercheck, coreFilter, minScore]
    orthoArgs = [strict, checkCoorthologsRef, rbh, rep, ignoreDistance, lowComplexityFilterOff, evalBlast, evalHmmer, evalRelaxfac, hitLimit, autoLimit, scoreThreshold, scoreCutoff, aligner, local, glocal, searchTaxa]
    otherArgs = [cpu, hyperthread, debug, silent]

    ### run fdog
    runSingle([basicArgs, ioArgs, pathArgs, coreArgs, orthoArgs, fasArgs, otherArgs, False])

if __name__ == '__main__':
    main()
