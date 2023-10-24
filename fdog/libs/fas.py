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
import shutil
import greedyFAS.annoFAS.annoModules as annoFas

import fdog.libs.zzz as general_fn


##### FUNCTIONS RELATED TO FAS #####

def check_fas_executable():
    try:
        subprocess.check_output(['fas.setup -t ./ --check'], shell = True, stderr = subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        print('\033[96m%s\033[0m' % e.output.decode(sys.stdout.encoding).strip())
        print('FAS installed but fas.setup still need to be run if you want to use it!')
        return(0)
    return(1)


def get_tool_fas_path():
    """ Get path to FAS annotation tools """
    cmd = 'fas.setup -t ~/ -c'
    try:
        out = subprocess.run(
                [cmd], shell = True, capture_output = True, check = True)
        tool_path = out.stdout.decode().split('\n')[0].split()[6].replace('.','')
        return(tool_path)
    except:
        sys.exit('ERROR: fas.setup cannot be called!')


def get_anno_fas(seqName, spec, seq_id, seq, hmmpath, annopath):
    """ Get annotation for a seq_id from existing json file in annopath """
    out_json = '%s/%s/%s_%s.json' % (hmmpath, seqName, spec, seq_id)
    if not os.path.exists(out_json):
        tmp_seed_fa = '%s/%s/%s_%s.fa' % (hmmpath, seqName, seqName, spec)
        with open(tmp_seed_fa, 'w') as tmp_seed_fa_out:
            tmp_seed_fa_out.write('>%s\n%s\n' % (seq_id, seq))
        spec_anno = '%s/%s.json' % (annopath, spec)
        try:
            anno_dict = annoFas.extractAnno(tmp_seed_fa, spec_anno)
            anno_dict['clan'] = annoFas.getClans(
                get_tool_fas_path(), anno_dict['feature'])
            annoFas.save2json(
                anno_dict, '%s_%s' % (spec, seq_id), '%s/%s' % (hmmpath, seqName))
            os.remove(tmp_seed_fa)
        except:
            sys.exit(
                'ERROR: Annotation for %s cannot be found in %s'
                % (seq_id, spec_anno))
    return(out_json)


def calc_pairwise_fas(seed_json, query_json, seqName, hmmpath):
    """ Calculate FAS score for a pair seed and query protein
    Input are two anno json files for seed and query
    Return a value between 0 and 1
    """
    general_fn.check_file_exist(seed_json)
    general_fn.check_file_exist(query_json)

    fas_cmd = 'fas.run -s %s -q %s --no_config' % (seed_json, query_json)
    fas_cmd = '%s -a %s/%s --raw --tsv --domain --cpus 1 -o %s/%s' \
                % (fas_cmd, hmmpath, seqName, hmmpath, seqName)
    try:
        fas_out = subprocess.run(
                [fas_cmd], shell = True, capture_output = True, check = True)
    except:
        sys.exit('ERROR: Error running FAS\n%s' % fas_cmd)
    results = fas_out.stdout.decode().split('\n')
    for l in results:
        if l.startswith('#') and len(l.split('\t')) > 1:
            return(l.split('\t')[-1])
    return('')


def calc_fas_cand(args):
    """ Calculate FAS score for a ortholog candidate against seed
    Ortholog candidate defined by spec, seq_id and seq
    """
    (fasOff, seqName, seed_json, spec, seq_id, seq, hmmpath, annopath) = args
    if not fasOff == True:
        query_json = get_anno_fas(seqName, spec, seq_id, seq, hmmpath, annopath)
        fas_score = calc_pairwise_fas(seed_json, query_json, seqName, hmmpath)
    else:
        fas_score = 1
    return(fas_score)


def calc_fas_multi (input_fa, outpath, annopath, cpus):
    """ Calculate pairwise FAS scores for all orthologs vs seed protein
    input_fa is the default <seqName>.extended.fa output file of fDOG
    Output will be <seqName>_forward.
    """
    fasCmd = 'fas.runFdogFas -i %s -w %s --cores %s --redo_anno' % (input_fa, annopath, cpus)
    try:
        subprocess.call([fasCmd], shell = True)
        if os.path.exists(outpath + '/tmp'):
            shutil.rmtree(outpath + '/tmp')
    except:
        sys.exit('Problem running\n%s' % (fasCmd))
