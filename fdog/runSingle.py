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
import argparse
import subprocess
from ete3 import NCBITaxa
from pkg_resources import get_distribution
import time

import fdog.libs.preparation as prepare_fn
import fdog.libs.orthosearch as ortho_fn
import fdog.libs.corecompile as core_fn
import fdog.libs.fas as fas_fn
import fdog.libs.tree as tree_fn
import fdog.libs.output as output_fn


def main():
    version = get_distribution('fdog').version
    parser = argparse.ArgumentParser(description='You are running fDOG version ' + str(version) + '.',
                                     epilog="For more information on certain options, please refer to the wiki pages "
                                            "on github: https://github.com/BIONF/fDOG/wiki")
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--seqFile', help='Input file containing the seed sequence (protein only) in fasta format',
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
    seqFile = args.seqFile
    seqName = args.jobName
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
    if seqFile == 'infile.fa':
        fdogPath = os.path.realpath(__file__).replace('/runSingle.py','')
        seqFile = '%s/data/infile.fa' % fdogPath

    (seqFile, hmmpath, corepath, searchpath, annopath) = prepare_fn.check_input(
                    [seqFile, refspec, outpath, hmmpath,
                    corepath, searchpath, annopath, pathFile])
    pathArgs = [outpath, hmmpath, corepath, searchpath, annopath]

    if not fasOff:
        check_fas = fas_fn.check_fas_executable()
        if check_fas == 0:
            sys.exit('ERROR: FAS is not executable! You still can use fDOG with --fasOff!')

    ##### Identify seed ID from refspec genome
    if reuseCore:
        core_fa = '%s/%s/%s.fa' % (hmmpath, seqName, seqName)
        seed_id = prepare_fn.get_seed_id_from_fa(core_fa, refspec)
    else:
        seed_id = prepare_fn.identify_seed_id(seqFile, refspec, corepath, debug, silentOff)
    print('Identified seed ID: %s' % seed_id)

    ##### DO CORE COMPILATION
    # start = time.time()
    coreArgs = [minDist, maxDist, coreSize, coreTaxa, distDeviation,
                alnStrategy, fasOff]
    orthoCoreArgs = [CorecheckCoorthologsOff, rbh, True, evalBlast,
                    lowComplexityFilter, evalHmmer, coreHitLimit,
                    scoreCutoff, aligner] # rep = True
    otherCoreArgs = [cpus, debugCore, silentOff, noCleanup, force, append]
    print('Compiling core set for %s' % seqName)
    core_runtime = core_fn.run_compile_core([seqFile, seqName, refspec, seed_id, reuseCore,
        forceCore, coreArgs, pathArgs, orthoCoreArgs, otherCoreArgs, debug])
    print('==> Core compilation finished in %s' % core_runtime[1])


    ##### DO ORTHOLOG SEARCH USING CORE HMM (HAMSTR)
    if not coreOnly:
        start = time.time()
        # check existing output
        finalOutfile = '%s/%s.extended.fa' % (outpath, seqName)
        finalOutfile = os.path.abspath(finalOutfile)
        output_fn.check_output_exist(finalOutfile, force, append)
        # get list of search taxa
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
                output_fn.print_debug(debugCore, 'Search taxa', searchTaxa)
            if len(searchTaxa) == 0:
                exit('ERROR: No taxon found within %s taxonomy group!' % group)
            else:
                searchTaxa = ','.join(searchTaxa)
        # do ortholog search
        orthoArgs = [checkCoorthologsRefOff, rbh, rep, evalBlast,
                    lowComplexityFilter, evalHmmer, hitLimit, scoreCutoff, aligner]
        otherArgs = [searchTaxa, cpus, debug, silentOff, noCleanup, force, append]
        hamstr_out = ortho_fn.run_hamstr([seqName, refspec, pathArgs, orthoArgs, otherArgs])
        output_fn.write_hamstr(hamstr_out, outpath, seqName, force, append)
        end = time.time()
        print('==> Ortholog search finished in ' + '{:5.3f}s'.format(end - start))

        ##### DO FINAL FAS CALCULATION
        if not fasOff:
            try:
                fasVersion = subprocess.run(['fas.run --version'], shell = True, capture_output = True, check = True)
            except:
                sys.exit('Problem with FAS! Please check https://github.com/BIONF/FAS or turn it off if not needed!')
            if os.path.exists(finalOutfile):
                start = time.time()
                fas_fn.calc_fas_multi(finalOutfile, outpath, annopath, cpus)
                end = time.time()
                print('==> FAS calculation finished in ' + '{:5.3f}s'.format(end - start))
        else:
            output_fn.hamstr_2_profile(finalOutfile)

        end = time.time()
        print('==> fdog.run finished in ' + '{:5.3f}s'.format(end - begin))

if __name__ == '__main__':
    main()
