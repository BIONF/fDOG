# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2021 Hannah Muelbaier
#
#  This script is used to prepare the core group used as input for fDOG-Assembly from a fasta file of an ortholog group.
#
#  This script is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License <http://www.gnu.org/licenses/> for
#  more details
#
#  Contact: hannah.muelbaier@gmail.com
#
#######################################################################

############################ imports ##################################
import os
import os.path
import sys
import argparse
import fdog.libs.alignment as align_fn
import fdog.libs.zzz as general_fn

def check_fasta(file):
    nHeader = general_fn.count_line(file, '>', True)
    nSeq = general_fn.count_line(file, '>', False)
    if not nHeader == nSeq:
        return(1)
    return(0)

def make_single_line_fasta(input, gene, out_folder):
    print(out_folder)
    output = out_folder + gene + ".fa"
    print(output)
    with open(input, 'r') as f_input, open(output, 'w') as f_output:
        block = []
        for line in f_input:
            if line.startswith('>'):
                if block:
                    f_output.write(''.join(block) + '\n')
                    block = []
                f_output.write(line)
            else:
                block.append(line.strip())

        if block:
            f_output.write(''.join(block) + '\n')
    return (output)

def makeMSA(out_folder, gene, fasta_file):
    aln_file = out_folder + gene + '.aln'
    if align_fn.get_muscle_version('muscle') == 'v3':
        os.system('muscle -quiet -in %s -out %s' % (fasta_file, aln_file))
        #print("muscle -quiet -in " + output_file + " -out " + aln_file)
    else:
        os.system('muscle -quiet -align %s -out %s' % (fasta_file, aln_file))
    return aln_file

def makeHMM(out_folder, gene, aln_file):
    hmm_dir = out_folder + 'hmm_dir'
    os.system('mkdir %s >/dev/null 2>&1' % (hmm_dir))
    out_file = '%s/%s.hmm' % (hmm_dir, gene)
    hmmbuild_cmd = 'hmmbuild --amino %s %s' % (out_file, aln_file)
    os.system(hmmbuild_cmd)
    return out_file


def main():

    #################### handle user input #####################################
    version = '0.0.1'
    ################### initialize parser ######################################
    parser = argparse.ArgumentParser(description='You are running fdog.addCoreGroup version ' + str(version) + '.')
    ################## required arguments ######################################
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--fasta', help='Path to fasta file of ortholog group.', action='store', default='', required=True)
    required.add_argument('--out', help='Path to output folder.', action='store', default='', required=True)
    required.add_argument('--geneName', help='Path to output folder.', action='store', default='', required=True)
    args = parser.parse_args()

    fasta_file_input = args.fasta
    out_folder = args.out
    gene = args.geneName


    out_folder = out_folder + '/' + gene + '/'
    os.system('mkdir %s >/dev/null 2>&1' % (out_folder))

    if check_fasta(fasta_file_input) == 1:
        fasta_file = make_single_line_fasta(fasta_file_input, gene, out_folder)
    else:
        fasta_file = out_folder + gene + '.fa'
        os.system('cp ' + fasta_file_input + ' ' + fasta_file)

    aln_file = makeMSA(out_folder, gene, fasta_file)
    hmm_file = makeHMM(out_folder, gene, aln_file)

    print('Core group located at %s. Fasta file: %s; MSA: %s; HMM: %s' % (out_folder, fasta_file, aln_file, hmm_file))

main()
