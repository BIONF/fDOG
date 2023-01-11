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


##### GENERAL FUNCTIONS FOR FILES, FOLDERS AND GENERAL VARIABLES #####

def check_file_exist(file):
    """ Exit if a file does not exist"""
    if not os.path.exists(os.path.abspath(file)):
        sys.exit('%s not found' % file)


def read_file(file):
    """ Read a file and return list of lines"""
    if os.path.exists(file):
        with open(file, 'r') as f:
            lines = f.read().splitlines()
            f.close()
            return(lines)
    else:
        sys.exit('%s not found' % file)


def read_dir(dir):
    """ Return list of directories from a given path """
    check_file_exist(dir)
    out_dirs = []
    p = os.listdir(dir)
    for i in p:
        if os.path.isdir('%s/%s' % (dir, i)):
            out_dirs.append(i)
    return(out_dirs)


def get_ids_from_folder(folder, type):
    """ Get taxonomy IDs for from blast_dir, genome_dir or weight_dir
    Return dictionary {taxID:<TaxName>@<TaxID>@Ver}
    """
    tax_ids = {}

    for name in read_dir(folder):
        if type == 'weight_dir':
            if not name.endswith('.json'):
                continue
            else:
                name = name.replace('.json','')
        else:
            if not os.path.isdir('%s/%s' % (folder, name)):
                continue
        id = name.split('@')[1]
        if not id in tax_ids:
            tax_ids[id] = name
    return(tax_ids)


def load_config(config_file):
    """ Load a YAML file and return as a dictionary """
    with open(config_file, 'r') as stream:
        try:
            return yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)


def matching_elements(dictionary, search_string):
    """ Search for a string in dictionary's values
    Return {key:val} where string was found in val
    """
    return {key:val for key,val in dictionary.items() if search_string == val}


def remove_dup_in_dict(dictionary):
    """ Find and remove duplicated or empty values of a dictionary """
    tmp_dict = {'_'.join(val) : key for key, val in dictionary.items()}
    res = {val : key.split('_') for key, val in tmp_dict.items()}
    res = {key : val for key, val in res.items() if len(val[0]) > 0}
    return(res)
