import fdog.libs.addtaxon as addTaxon_fn
import fdog.libs.tree as tree_fn
import fdog.libs.zzz as general_fn
import sys
import os
import argparse

def check_fasta(file):
    nHeader = general_fn.count_line(file, '>', True)
    nSeq = general_fn.count_line(file, '>', False)
    if not nHeader == nSeq:
        return(1)
    nPipe = general_fn.count_line(file, '|', True)
    if nPipe > 0:
        return(1)
    return(0)

def check_path(path):
    if not os.path.exists(path):
        return False
    else:
        if os.path.isfile(path):
            return "File"
        else:
            return "Path"

def parse_file(path):
    file = open(path, "r")
    lines = file.readlines()
    id_dict = {}
    for line in lines:
        line = line.rstrip()
        ncbi, name = line.split("\t")
        id_dict[ncbi] = name

    return id_dict


def main():

    #################### handle user input #####################################
    version = '0.0.1'
    ################### initialize parser ######################################
    parser = argparse.ArgumentParser(description='You are running fdog.addAssembly version ' + str(version) + '.')
    ################## required arguments ######################################
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--fasta', help='Path to fasta file or folder', action='store', default='', required=True)
    required.add_argument('--out', help='Path to output folder.', action='store', default='', required=True)
    required.add_argument('--ncbi', help='NCBI number of species or mapping file', action='store', default='', required=True)
    required.add_argument('--ver', help='Version', action='store', default='', required=True)
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--link', help='link files not copy', action='store_true', default = True)

    args = parser.parse_args()

    fasta = args.fasta
    if check_path(fasta) == False:
        print("%s does not exists. Exiting ..."%(fasta))
        sys.exit()
    else:
        format = check_path(fasta)
    out_folder = args.out
    out_folder = os.path.abspath(out_folder) + '/'
    os.system('mkdir %s >/dev/null 2>&1' % (out_folder))
    ncbi = args.ncbi
    ver = args.ver
    ln = args.link
    id_dict = {}

    if check_path(ncbi) == False:
        if isdigit(ncbi) and format == "File":
            id_dict[ncbi] = fasta_file
        else:
            print("%s is no file or digit. Exiting ..."%(ncbi))
            sys.exit()
    elif check_path(ncbi) == "File":
        id_dict = parse_file(ncbi)
    else:
        print("%s is no file or digit. Exiting ..."%(ncbi))
        sys.exit()

    if format == "File":
        fa = id_dict[id]
        if check_fasta(fa):
            name = addTaxon_fn.generate_spec_name(id, "", ver)
            if ln == False:
                assembly_folder = out_folder + name
                os.system('mkdir %s >/dev/null 2>&1' % (assembly_folder))
                os.system("cp %fa %s/%s.fa" %(fa, assembly_folder, name))
            else:
                assembly_folder = out_folder + name
                os.system('mkdir %s >/dev/null 2>&1' % (assembly_folder))
                os.system("ln %fa %s/%s.fa" %(fa, assembly_folder, name))
        else:
            print("%s Fasta format not valid or header includes |"%(fa))

    for id in id_dict:
        fa = id_dict[id]
        fasta = os.path.abspath(fasta) + '/'
        if check_fasta(fasta + fa):
            name = addTaxon_fn.generate_spec_name(id, "", ver)
            if ln == False:
                assembly_folder = out_folder + name
                os.system('mkdir %s >/dev/null 2>&1' % (assembly_folder))
                os.system("cp %s/%fa %s/%s.fa" %(fasta, fa, assembly_folder, name))
            else:
                assembly_folder = out_folder + name
                os.system('mkdir %s >/dev/null 2>&1' % (assembly_folder))
                os.system("ln -s %s/%s %s/%s.fa" %(fasta, fa, assembly_folder, name))
        else:
            print("%s Fasta format not valid or header includes |"%(fa))

    print("DONE, files can be found: %s"%(out_folder))

main()
