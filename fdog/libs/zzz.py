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
import ssl
import urllib.request
import yaml
import time
import pickle


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


def load_config(config_file):
    """ Load a YAML file and return as a dictionary """
    with open(config_file, 'r') as stream:
        try:
            return yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)


def download_progress(count, block_size, total_size):
    global start_time
    if count == 0:
        start_time = time.time()
        return
    duration = time.time() - start_time
    progress_size = int(count * block_size)
    speed = int(progress_size / (1024 * duration))
    percent = int(count * block_size * 100 / total_size)
    if percent > 100:
        percent = 100
    sys.stdout.write("\r...%d%%, %d MB, %d KB/s, %d seconds passed" %
                     (percent, progress_size / (1024 * 1024), speed, duration))
    sys.stdout.flush()


def download_file(url, file):
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    download_file = urllib.request.URLopener(context=ctx)
    print('Downloading %s' % (url + '/' + file))
    urllib.request.urlretrieve(url + '/' + file, file, download_progress)
    print(' ... done!')


def count_line(file, pattern, contain):
    """ Count lines in file that contain (or not) a pattern """
    nline = 0
    with open(file, 'r') as f:
        for line in f:
            if contain:
                if pattern in line:
                    nline = nline + 1
            else:
                if not pattern in line:
                    nline = nline + 1
    return(nline)


def get_ids_from_folder(folder, type):
    """ Get taxonomy IDs for from coreTaxa_dir, searchTaxa_dir or annotation_dir
    Return dictionary {taxID:<TaxName>@<TaxID>@Ver}
    """
    tax_ids = {}

    for name in read_dir(folder):
        if type == 'annotation_dir':
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


def join_2lists(first_list, second_list):
    """ Join 2 lists """
    in_first = set(first_list)
    in_second = set(second_list)
    in_second_but_not_in_first = in_second - in_first
    out = first_list + list(in_second_but_not_in_first)
    return(out)


def save_pyobj(obj, out_file):
    """ Save a python object to out_file """
    with open(out_file, 'wb') as obj_out:
        pickle.dump(obj, obj_out)


def read_pyobj_file(in_file):
    """ Read a python object from an in_file """
    with open(in_file, 'rb') as obj_file:
        return(pickle.load(obj_file))


def query_yes_no(question, default='yes'):
    valid = {'yes': True, 'y': True, 'ye': True,
             'no': False, 'n': False}
    if default is None:
        prompt = ' [y/n] '
    elif default == 'yes':
        prompt = ' [Y/n] '
    elif default == 'no':
        prompt = ' [y/N] '
    else:
        raise ValueError('invalid default answer: "%s"' % default)
    while True:
        choice = sys.stdin.readline().rstrip().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write('Please respond with "yes" or "no" '
                             '(or "y" or "n").\n')


def search_string_in_file(file, string):
    """ Search for a string in file
    Return 0 if not found, 1 if found
    """
    flag = 0
    with open(file, 'r') as fp:
        for l_no, line in enumerate(fp):
            if string in line:
                flag = 1
                break
    return(flag)
