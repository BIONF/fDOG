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

from pysam import FastaFile
from Bio import SeqIO

import fdog.libs.zzz as general_fn


##### FUNCTIONS RELATED TO FASTA SEQ #####

def add_seq_to_dict(dict, id, seq):
    """ Add fasta sequence to a dictionary """
    if not id in dict:
        dict[id] = seq
    return(dict)


def read_fasta(fa_file):
    """ Read LARGE fasta file and return fasta object
    Sequence can be get using fasta_object.fetch(seq_id)
    """
    fasta_object = FastaFile(fa_file)
    return(fasta_object)


def write_fasta(fa_dict, out_file):
    """ Write sequences in SeqIO dict into output file """
    with open(out_file, 'w') as out:
        for seq in fa_dict:
            out.write('>%s\n%s\n' % (seq, fa_dict[seq].seq))


def append_to_fasta_file(fa_file, new_fa_dict):
    """ Append a dict of fasta seq to an existing fasta file """
    general_fn.check_file_exist(fa_file)
    existing_seq = SeqIO.to_dict(SeqIO.parse(open(fa_file),'fasta'))
    with open(fa_file, 'a') as fa_out:
        for id, seq in new_fa_dict.items():
            if not id in existing_seq:
                fa_out.write('>%s\n%s\n' % (id, seq))


def check_long_seq(fa_file, max_len):
    """ Check if any sequence longer than max_len
    (12.000 aa/nt for muscle v3; 20.000 for muscle v5)"""
    fa_seq = SeqIO.parse(open(fa_file),'fasta')
    for fa in fa_seq:
        if len(fa.seq) > max_len:
            return(1)
    return(0)


def remove_dup(fa_file):
    """ Remove duplicated sequences (filter by headers) """
    tmp = {}
    fa_seq = SeqIO.parse(open(fa_file),'fasta')
    for fa in fa_seq:
        if not fa.id in tmp:
            tmp[fa.id] = fa.seq
    with open(fa_file, 'w') as out:
        for id, seq in tmp.items():
            out.write('>%s\n%s\n' % (id, seq))
