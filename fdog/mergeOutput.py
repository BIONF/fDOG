# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2022 Vinh Tran
#
#  This script is used to merge all output files (.extended.fa, .phyloprofile,
#  _forward.domains, _reverse.domains) in a given directory into one file each.
#
#  This script is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License <http://www.gnu.org/licenses/> for
#  more details
#
#  Contact: dosch@bio.uni-frankfurt.de
#
#######################################################################

import sys
import os
from os import listdir as ldir
import argparse
import yaml
from importlib.metadata import version, PackageNotFoundError
from Bio import SeqIO

def createConfigPP(phyloprofile, domains_0, ex_fasta, directory, out):
    settings = dict(
        mainInput = '%s/%s' % (directory, phyloprofile),
        fastaInput = '%s/%s' % (directory, ex_fasta)
    )
    if not domains_0 == None:
        settings['domainInput'] = '%s/%s' % (directory, domains_0)
    settings['clusterProfile'] = 'TRUE'
    with open('%s.config.yml' % (out), 'w') as outfile:
        yaml.dump(settings, outfile, default_flow_style = False)


def main():
    fdog_version = version("fdog")
    parser = argparse.ArgumentParser(description='You are running fDOG version ' + str(fdog_version) + '.')
    parser.add_argument('-i', '--input',
                        help='Input directory, where all single output (.extended.fa, .phyloprofile, _forward.domains, _reverse.domains) can be found',
                        action='store', default='', required=True)
    parser.add_argument('-o', '--output', help='Output name', action='store', default='', required=True)
    args = parser.parse_args()

    directory = args.input
    out = args.output
    if not os.path.exists(os.path.abspath(directory)):
        sys.exit('%s not found' % directory)
    else:
        directory = os.path.abspath(directory)

    phyloprofile = None
    domains_0 = None
    domains_1 = None
    domains_2 = None
    domains_3 = None
    ex_fasta = None
    lines_seen = set()
    lines_seen_2 = set()
    lines_seen_3 = set()
    lines_seen_4 = set()
    fa_seq_id = set()
    for infile in ldir(directory):
        if infile.endswith('.phyloprofile') and not infile == out + '.phyloprofile':
            if not phyloprofile:
                phyloprofile = out + '.phyloprofile'
                phyloprofile_out = open(phyloprofile, 'w')
            with open(directory + '/' + infile, 'r') as reader:
                lines = reader.readlines()
                for line in lines:
                    if line not in lines_seen: # not a duplicate
                        phyloprofile_out.write(line)
                        lines_seen.add(line)
        elif infile.endswith('_forward.domains') and not infile == out + '_forward.domains':
            if not domains_0:
                domains_0 = out + '_forward.domains'
                domains_0_out = open(domains_0, 'w')
            with open(directory + '/' + infile, 'r') as reader:
                lines = reader.readlines()
                for line in lines:
                    if line not in lines_seen: # not a duplicate
                        domains_0_out.write(line)
                        lines_seen.add(line)
        elif infile.endswith('_reverse.domains') and not infile == out + '_reverse.domains':
            if not domains_1:
                domains_1 = open(out + '_reverse.domains', 'w')
            with open(directory + '/' + infile, 'r') as reader:
                lines = reader.readlines()
                for line in lines:
                    if line not in lines_seen_2: # not a duplicate
                        domains_1.write(line)
                        lines_seen_2.add(line)
        elif infile.endswith('.forward.domains') and not infile == out + '.forward.domains':
            if not domains_2:
                domains_2 = out + '.forward.domains'
                domains_2_out = open(domains_2, 'w')
            with open(directory + '/' + infile, 'r') as reader:
                lines = reader.readlines()
                for line in lines:
                    if line not in lines_seen_3: # not a duplicate
                        domains_2_out.write(line)
                        lines_seen_3.add(line)
        elif infile.endswith('.reverse.domains') and not infile == out + '.reverse.domains':
            if not domains_3:
                domains_3 = open(out + '.reverse.domains', 'w')
            with open(directory + '/' + infile, 'r') as reader:
                lines = reader.readlines()
                for line in lines:
                    if line not in lines_seen_4: # not a duplicate
                        domains_3.write(line)
                        lines_seen_4.add(line)
        elif infile.endswith('.extended.fa') and not infile == out + '.extended.fa':
            if not ex_fasta:
                ex_fasta = out + '.extended.fa'
                ex_fasta_out = open(ex_fasta, 'w')
            inSeq = SeqIO.to_dict((SeqIO.parse(open(directory + '/' + infile), 'fasta')))
            for seq in inSeq:
                if not seq in fa_seq_id:
                    ex_fasta_out.write('>%s\n%s\n' % (seq, inSeq[seq].seq))
                    fa_seq_id.add(seq)

    if phyloprofile:
        phyloprofile_out.close()
    if domains_0:
        domains_0_out.close()
    if domains_1:
        domains_1.close()
    if ex_fasta:
        ex_fasta_out.close()

    createConfigPP(phyloprofile, domains_0, ex_fasta, directory, out)
    print('Done! Output files:\n%s/%s.*' % (directory,out))


if __name__ == "__main__":
    sys.exit(main())
