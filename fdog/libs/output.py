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
from Bio import SeqIO

import fdog.libs.zzz as general_fn


##### FUNCTIONS FOR OUTPUT #####

def print_debug(debug, cat, msg):
    """ Print msg of a category in debug mode """
    if debug == True:
        if cat == '':
            print('#DEBUG#\t%s' % msg)
        else:
            print('#DEBUG#\t%s\n#DEBUG#\t%s' % (cat, msg))


def print_stdout(silentOff, msg):
    """ Print stdout """
    if silentOff == True:
        print(msg)


def check_output_exist(outfile, force, append):
    """ Check if outfile exists
    And decide depends on the choice of force or append option
    """
    if os.path.exists(outfile):
        if force == True:
            print('WARNING: %s will be deleted!' % outfile)
            os.remove(outfile)
        elif append == True:
            general_fn.check_file_exist(outfile)
            print('Result will be appended to %s!' % outfile)
        else:
            sys.exit(
                'WARNING: %s exists! ' % outfile
                + 'You still can run with --force or --append option')


def write_hamstr(hamstr_result, outpath, seqName, force, append):
    """ Write result of ortholog search into seqName.extended.fa """
    outfile = '%s/%s.extended.fa' % (outpath, seqName)
    outfile = os.path.abspath(outfile)
    check_output_exist(outfile, force, append)

    ### Write to output.extended.fa
    ortho_count = len(hamstr_result) - 1
    if append == True:
        if os.path.exists(outfile):
            old_result_tmp = SeqIO.to_dict((SeqIO.parse(open(outfile),'fasta')))
            old_result = {}
            for old_id in old_result_tmp:
                if not old_id in hamstr_result:
                    old_result[old_id] = str(old_result_tmp[old_id].seq)
            hamstr_result = {**hamstr_result, **old_result}

    with open(outfile, 'w') as out_file:
        for id, seq in hamstr_result.items():
            out_file.write('>%s\n%s\n' % (id, seq))
    return(
        'Found %s ortholog(s)!\nOutput file: %s' % (ortho_count, outfile))


def hamstr_2_profile(fa_file):
    """ Convert extended.fa file into phyloprofile file """
    if os.path.exists(fa_file):
        pp_file = fa_file.replace('.extended.fa', '.phyloprofile')
        fa = SeqIO.to_dict((SeqIO.parse(open(fa_file),'fasta')))
        with open(pp_file, 'w') as pp:
            pp.write('geneID\tncbiID\torthoID\n')
            for id in list(fa.keys()):
                tmp = id.split('|')
                pp.write('%s\tncbi%s\t%s\n' % (tmp[0], tmp[1].split('@')[1], id))


def add_all_taxa(pp_file, searchTaxa):
    """ Add all "missing" search taxa into phyloprofile file """
    missing_taxa = [] # missing_taxa = [ncbi_id]
    for taxon in searchTaxa.split(','):
        flag = general_fn.search_string_in_file(pp_file, taxon)
        if flag == 0:
            missing_taxa.append(taxon.split('@')[1])
    first_gene = ''
    if os.path.exists(pp_file):
        with open(pp_file, 'a') as pp:
            for line in general_fn.read_file(pp_file):
                if not line.startswith('geneID'):
                    if not first_gene:
                        first_gene = line.split('\t')[0]
                        for i in missing_taxa:
                            if len(line.split('\t')) == 5:
                                pp.write(f'{first_gene}\tncbi{i}\tfdogMA\tNA\tNA\n')
                            else:
                                pp.write(f'{first_gene}\tncbi{i}\tfdogMA\n')
                        break
