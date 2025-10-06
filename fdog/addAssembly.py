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
    file.close()
    return id_dict


def main():
    print("#################################")
    #################### handle user input #####################################
    version = '0.0.3'
    ################### initialize parser ######################################
    parser = argparse.ArgumentParser(description='You are running fdog.addAssembly version ' + str(version) + '.')
    ################## required arguments ######################################
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--fasta', help='Path to fasta file or folder', action='store', default='', required=True)
    required.add_argument('--out', help='Path to output folder.', action='store', default='', required=True)
    required.add_argument('--ncbi', help='NCBI ID of species or a mapping file (tab separated) containing the NCBI ID and the corresponding file name placed in the folder given by --fasta. ', action='store', default='', required=True)
    required.add_argument('--ver', help='Version', action='store', default='', required=True)
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--link', help='links fasta files instead of copying them', action='store_true', default = False)

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
        if ncbi.isdigit() and format == "File":
            id_dict[ncbi] = fasta
        else:
            print("%s is no file or digit. Exiting ..."%(ncbi))
            sys.exit()
    elif check_path(ncbi) == "File":
        print("Parsing mapping file ...")
        id_dict = parse_file(ncbi)
        print("... done")
    else:
        print("%s is no file or digit. Exiting ..."%(ncbi))
        sys.exit()
    #print(format)
    #print(fasta)
    if format == "File":
        fa = id_dict[ncbi]
        if check_fasta(fa):
            name = addTaxon_fn.generate_spec_name(ncbi, "", ver)
            if ln == False:
                assembly_folder = out_folder + name
                os.system('mkdir %s >/dev/null 2>&1' % (assembly_folder))
                os.system("cp %s %s/%s.fa" %(fa, assembly_folder, name))
            else:
                assembly_folder = out_folder + name
                os.system('mkdir %s >/dev/null 2>&1' % (assembly_folder))
                os.system("ln %s %s/%s.fa" %(fa, assembly_folder, name))
        else:
            print("%s Fasta format not valid or header includes |"%(fa))
    
    else:
        for sp in id_dict:
            print("Adding species %s"%(sp))
            #print(id_dict)
            fa = id_dict[sp]
            fasta = os.path.abspath(fasta) + '/'
            #print(fa)
            #print(fasta)
            fasta_path = fasta + fa
            if check_fasta(fasta_path):
                name = addTaxon_fn.generate_spec_name(sp, "", ver)
                if ln == False:
                    assembly_folder = out_folder + name
                    os.system('mkdir %s >/dev/null 2>&1' % (assembly_folder))
                    os.system("cp %s %s/%s.fa" %(fasta_path, assembly_folder, name))
                else:
                    assembly_folder = out_folder + name
                    os.system('mkdir %s >/dev/null 2>&1' % (assembly_folder))
                    os.system("ln -s %s %s/%s.fa" %(fasta_path, assembly_folder, name))
            else:
                print("%s Fasta format not valid or header includes |"%(fasta_path))

    print("DONE, files can be found: %s"%(out_folder))

if __name__ == '__main__':
    main()
