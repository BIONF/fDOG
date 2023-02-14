# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2022 Vinh Tran
#
#  This script is used to setup fdog: install dependencies and
#  download pre-computed data
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
import platform
import argparse
import subprocess
import shutil
from ete3 import NCBITaxa
from pathlib import Path
from pkg_resources import get_distribution

import fdog.libs.zzz as general_fn
import fdog.libs.fas as fas_fn
import fdog.libs.alignment as align_fn


def check_conda_env():
    """ Return if a conda env is currently using """
    if 'CONDA_DEFAULT_ENV' in os.environ:
        if not os.environ['CONDA_DEFAULT_ENV'] == 'base':
            return(True)
    return(False)


def get_source_path():
    """ Get path of installed fDOG library """
    fdogPath = os.path.realpath(__file__).replace('/setupfDog.py','')
    return(fdogPath)


def get_data_path(fdogPath):
    """ Get path of fDOG data """
    pathconfigFile = fdogPath + '/bin/pathconfig.yml'
    if not os.path.exists(pathconfigFile):
        sys.exit('No pathconfig.yml found. Please run fdog.setup (https://github.com/BIONF/fDOG/wiki/Installation#setup-fdog).')
    else:
        cfg = general_fn.load_config(pathconfigFile)
        try:
            dataPath = cfg['datapath']
        except:
            try:
                corepath = cfg['corepath']
            except:
                pass
            try:
                searchpath = cfg['searchpath']
            except:
                pass
            try:
                annopath = cfg['annopath']
            except:
                pass
            dataPath = 'Core taxa: %s\nSearch taxa: %s\nAnnotations: %s' % (corepath, searchpath, annopath)
        return(dataPath)


def install_fas(woFAS):
    """ Install greedyFAS """
    if not woFAS:
        ### check if fas already installed
        try:
            fasVersion = subprocess.run(['fas.run --version'], shell = True, capture_output = True, check = True)
        except:
            print('=> greedyFAS (https://github.com/BIONF/FAS)')
            install_fas_cmd = 'pip install greedyFAS'
            try:
                subprocess.check_output([install_fas_cmd], shell = True, stderr = subprocess.STDOUT)
            except subprocess.CalledProcessError as e:
                sys.exit('\033[91mERROR: Problem with installing FAS! Please do it manually. See: https://github.com/BIONF/FAS!\033[0m')
        ### check if fas installed but not yet configured
        check_fas = fas_fn.check_fas_executable()


def install_fasta36(fdogPath, cwd):
    """ Install FASTA36 from source """
    try:
        subprocess.check_output(['which fasta36'], shell = True, stderr = subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        print('=> FASTA36 (https://github.com/wrpearson/fasta36)')
        fasta36v = '36.3.8h_04-May-2020'
        fasta36url = 'https://github.com/wrpearson/fasta36/archive/refs/tags'
        fasta36file = 'v%s.tar.gz' % fasta36v
        if not os.path.exists('%s/bin/aligner/bin/fasta36' % fdogPath):
            if os.path.exists('%s/bin/aligner' % fdogPath):
                shutil.rmtree('%s/bin/aligner' % fdogPath)
            general_fn.download_file(fasta36url, fasta36file)
            shutil.unpack_archive(fasta36file, '%s/bin/' % fdogPath, 'gztar')
            os.remove(fasta36file)
            shutil.move('%s/bin/fasta36-%s' % (fdogPath, fasta36v), '%s/bin/aligner' % fdogPath)
            if 'Darwin' in platform.uname():
                make_cmd = 'make -f %s/bin/aligner/make/Makefile.os_x86_64 all' % fdogPath
            elif 'Linux' in platform.uname():
                make_cmd = 'make -f %s/bin/aligner/make/Makefile.linux64_sse2 all' % fdogPath
            else:
                sys.exit('\033[91mERROR: Cannot identify type of system (neither Linux nor Darwin/MacOS)\033[0m')
            try:
                print('Compiling fasta36. Please wait...')
                os.chdir('%s/bin/aligner/src' % fdogPath)
                subprocess.run(make_cmd, shell = True, check = True)
            except:
                sys.exit('\033[91mERROR: Cannot install FASTA36!\033[0m')
            os.chdir(cwd)
            if not os.path.exists('%s/bin/aligner/bin/fasta36' % fdogPath):
                sys.exit('\033[91mERROR: fasta36 not found! Please install it manually!\033[0m')
            else:
                print('FASTA36 installed at %s/bin/aligner/' % fdogPath)
        else:
            fasta36_path = align_fn.check_fasta36_executable(fdogPath)
            print('FASTA36 found at %s' % fasta36_path)


def check_dependencies(fdogPath):
    """ Check for missing dependencies
    Dependencies are specified in fdog/data/dependencies.txt file
    """
    missing = []
    dependencies = '%s/data/dependencies.txt' % fdogPath
    for tool in general_fn.read_file(dependencies):
        function = tool
        if tool == 'hmmer':
            function = 'hmmsearch'
        if tool == 'ncbi-blast+':
            function = 'blastp'
        try:
            subprocess.check_output(['which %s' % function], shell = True, stderr = subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            missing.append(tool)
    return(missing)


def download_data(dataPath, resetData):
    """ Downloade pre-calculated fDOG data """
    data_fdog_file = "data_HaMStR-2019c.tar.gz"
    checksum_data = "1748371655 621731824 $data_fdog_file"

    genome_path = '%s/searchTaxa_dir' % dataPath
    Path(genome_path).mkdir(parents = True, exist_ok = True)

    if len(general_fn.read_dir(genome_path)) < 1 or resetData:
        data_url = 'https://applbio.biologie.uni-frankfurt.de/download/hamstr_qfo'
        if os.path.exists(data_fdog_file) and resetData:
            os.remove(data_fdog_file)
        general_fn.download_file(data_url, data_fdog_file)
        try:
            print('Extracting %s...' % data_fdog_file)
            shutil.unpack_archive(data_fdog_file, dataPath, 'gztar')
        except:
            sys.exit('\033[91mERROR: Cannot extract %s to %s!\033[0m' % (data_fdog_file, dataPath))
        if 'genome_dir' in general_fn.read_dir(dataPath):
            os.rename('%s/genome_dir' % dataPath, '%s/searchTaxa_dir' % dataPath)
            os.rename('%s/blast_dir' % dataPath, '%s/coreTaxa_dir' % dataPath)
            os.rename('%s/weight_dir' % dataPath, '%s/annotation_dir' % dataPath)
        check_cmd = 'fdog.checkData -s %s/searchTaxa_dir -c %s/coreTaxa_dir -a %s/annotation_dir --reblast' % (dataPath, dataPath, dataPath)
        try:
            print('Checking downloaded data...')
            subprocess.run([check_cmd], stdout = subprocess.DEVNULL, check = True, shell = True)
        except:
            print('\033[96mWARNING: Problem with validating downloaded data. Please run fdog.checkData manually!\033[0m')
        os.remove(data_fdog_file)
        print('fDOG data downloaded and saved at %s' % dataPath)
    else:
        print('fDOG data found at %s' % dataPath)


def write_pathconfig(fdogPath, dataPath):
    """ Write data directories to pathconfig file """
    Path('%s/bin' % fdogPath).mkdir(parents = True, exist_ok = True)
    pathconfigFile = '%s/bin/pathconfig.yml' % fdogPath
    if os.path.exists(pathconfigFile):
        os.remove(pathconfigFile)
    with open(pathconfigFile, 'w') as cf:
        cf.write('datapath: \'%s\'\n' % dataPath)
        cf.write('corepath: \'%s/coreTaxa_dir\'\n' % dataPath)
        cf.write('searchpath: \'%s/searchTaxa_dir\'\n' % dataPath)
        cf.write('annopath: \'%s/annotation_dir\'\n' % dataPath)


def main():
    version = get_distribution('fdog').version
    parser = argparse.ArgumentParser(description='You are running fDOG version ' + str(version) + '.')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    required.add_argument('-d', '--dataPath', help='Output path for fDOG data', action='store', default='', required=True)
    optional.add_argument('--getSourcepath', help='Get path to installed fdog package', action='store_true', default=False)
    optional.add_argument('--getDatapath', help='Get fDOG default data path', action='store_true', default=False)
    optional.add_argument('--woFAS', help='Do not install FAS (https://github.com/BIONF/FAS)', action='store_true', default=False)
    optional.add_argument('--force', help='Force installing', action='store_true', default=False)
    optional.add_argument('--resetData', help='Re-download precalculated fDOG data', action='store_true', default=False)

    ### parse arguments
    args = parser.parse_args()
    dataPath = os.path.abspath(args.dataPath)
    woFAS = args.woFAS
    force = args.force
    resetData = args.resetData


    ### get install path
    fdogPath = get_source_path()
    if args.getSourcepath:
        print(fdogPath)
        sys.exit()

    ### get data path
    if args.getDatapath:
        dataPath = get_data_path(fdogPath)
        print(dataPath)
        sys.exit()

    ### check if pathconfig file exists
    pathconfigFile = '%s/bin/pathconfig.yml' % fdogPath
    demo_cmd = 'fdog.run --seqFile infile.fa --jobName test --refspec HUMAN@9606@3'
    if os.path.exists(pathconfigFile) and not force:
        check_fas = 1
        if not woFAS:
            check_fas = fas_fn.check_fas_executable()
        if check_fas == 1:
            print('fDOG seems to be ready to use!')
            print('You can test fDOG using the following command:\n%s' % demo_cmd)
        else:
            print('fDOG seems to be ready to use without FAS!')
            print('You can test fDOG using the following command:\n%s --fasOff' % demo_cmd)
        sys.exit()

    ### get ncbi taxonomy database for ete3
    print('*** Creating local NCBI taxonomy database...')
    ncbi = NCBITaxa()

    ### install dependencies
    print('*** Installing dependencies...')
    ## FAS
    if not woFAS:
        install_fas(woFAS)
    ## hmmer, blast+, clustalw, mafft, muscle
    missing_tools = check_dependencies(fdogPath)
    if len(missing_tools) > 0:
        if check_conda_env() == True:
            req_file = '%s/data/conda_requirements.yml' % fdogPath
            print('=> Dependencies in %s' % req_file)
            conda_install_cmd = 'conda install -c bioconda --file %s -y' % (req_file)
            try:
                subprocess.call([conda_install_cmd], shell = True)
            except:
                sys.exit('\033[91mERROR: Cannot install conda packages in %s!\033[0m' % req_file)
        else:
            install_cmd = 'sudo apt-get install -y -qq <tool>'
            sys.exit('\033[91mERROR: Please install these tools manually:\n%s\nusing the command: %s!\033[0m' % (', '.join(missing_tools), install_cmd))
    else:
        print('=> Dependencies in %s/data/dependencies.txt already installed!' % fdogPath)
    ## fasta36
    install_fasta36(fdogPath, os.getcwd())

    ### download pre-calculated data
    print('*** Downloading precalculated data...')
    ### Remove data if resetData is used
    if resetData:
        if os.path.exists(dataPath):
            print('fDOG data found in %s will be deleted! Enter to continue.' % dataPath)
            if general_fn.query_yes_no(''):
                shutil.rmtree(dataPath)

    Path(dataPath).mkdir(parents = True, exist_ok = True)
    download_data(dataPath, resetData)

    ### create pathconfig file
    write_pathconfig(fdogPath, dataPath)

    print('\033[96m==> FINISHED! fDOG data can be found at %s\033[0m' % dataPath)
    print('You can test fDOG using the following command:\n%s' % demo_cmd)

if __name__ == '__main__':
    main()
