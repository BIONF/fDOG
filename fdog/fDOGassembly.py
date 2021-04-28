############################ imports ###########################################
import os
import os.path
import sys
from Bio import SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
import argparse
import yaml
import subprocess
########################### functions ##########################################
def load_config(config_file):
    with open(config_file, 'r') as stream:
        try:
            return yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

def starting_subprocess(cmd, mode):
    if mode == 'debug':
        result = subprocess.run(cmd, shell=True)
    elif mode == 'silent':
        result = subprocess.run(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell=True)
    elif mode == 'normal':
        result = subprocess.run(cmd, stdout = subprocess.PIPE, shell=True)

def merge(blast_results, insert_length):
    #merging overlapping and contigous candidate regions
    number_regions = 0
    insert_length = int(insert_length)
    for key in blast_results:
        locations = blast_results[key]
        locations = sorted(locations, key = lambda x: int(x[3]))
        #print("test")
        #print(locations)
        size_list = len(locations)
        j = 0
        while j < size_list-1:
            i = j + 1
            while i < size_list:
                if ((locations[j][0] < locations[i][0]) and (locations[j][1] > locations[i][0]) and (locations[j][5] == locations[i][5]) and (locations[i][5] == '+')):
                    #merge overlapping regions plus strand
                    locations[j][1] = max(locations[j][1], locations[i][1])
                    locations[j][2] = min(locations[j][2], locations[i][2])
                    locations.pop(i)
                    size_list -= 1
                    i -= 1
                elif ((locations[j][1] > locations[i][1]) and (locations[j][0] < locations[i][1]) and (locations[j][5] == locations[i][5]) and (locations[i][5] == '-')):
                    #merge overlapping regions minus strand
                    locations[j][0] = min(locations[j][0], locations[i][0])
                    locations[j][2] = min(locations[j][2], locations[i][2])
                    locations.pop(i)
                    size_list -= 1
                    i -= 1
                elif ((locations[j][0] < locations[i][0]) and (locations[i][0] - locations[j][1] <= 2*insert_length) and (locations[j][5] == locations[i][5]) and (locations[i][5] == '+')):
                    #merging consecutive regions, the distance between booth is not longer than a cutoff, plus strand
                    locations[j][1] = max(locations[j][1], locations[i][1])
                    locations[j][2] = min(locations[j][2], locations[i][2])
                    locations.pop(i)
                    size_list -= 1
                    i -=1
                elif ((locations[j][1] > locations[i][1]) and (locations[j][0] - locations[i][1] <= 2* insert_length) and (locations[j][5] == locations[i][5]) and (locations[i][5] == '-')):
                    #merging consecutive regions, the distance between booth is not longer than a cutoff, minus strand
                    locations[j][0] = min(locations[j][0], locations[i][0])
                    locations[j][2] = min(locations[j][2], locations[i][2])
                    locations.pop(i)
                    size_list -= 1
                    i -=1
                i += 1
            j += 1

        number_regions += len(locations)
        blast_results[key] = locations

    return blast_results, number_regions

def parse_blast(line, blast_results, cutoff):
    # format blast line:  <contig> <sstart> <send> <evalue> <qstart> <qend>
    # format dictionary: {node_name: [(<start>,<send>,evalue, <qstart>,<qend>,<strand>)]}
    line = line.replace("\n", "")
    line_info = line.split("\t")
    evalue = float(line_info[3])
    #cut off
    if evalue > cutoff:
        return blast_results, evalue
    #add region to dictionary
    else:
        node_name, sstart, send, qstart, qend = line_info[0], int(line_info[1]), int(line_info[2]), int(line_info[4]), int(line_info[5])
        split = node_name.split("|")
        # finding out on which strand tBLASTn found a hit
        if sstart < send:
            strand = "+"
        else:
            sstart = int(line_info[2])
            send = int(line_info[1])
            strand = "-"
        #creating a dictionary that inlcudes every tBLASTn that is better as the evalue cut-off
        if len(split) > 1:
            node_name = split[1]
        if node_name in blast_results:
            list = blast_results[node_name]
            list.append([int(sstart),int(send), evalue, int(qstart), int(qend), strand])
            blast_results[node_name] = list
        else:
            blast_results[node_name] = [[int(sstart),int(send), evalue, int(qstart), int(qend), strand]]

    return blast_results, evalue

def candidate_regions(intron_length, cutoff_evalue, tmp_path):
    ###################### extracting candidate regions ########################
    # info about output blast http://www.metagenomics.wiki/tools/blast/blastn-output-format-6
    blast_file = open(tmp_path + "/blast_results.out", "r")
    evalue = 0
    blast_results = {}
    #parsing blast output
    while True:
        line = blast_file.readline()
        #end of file is reached
        if not line:
            break
        #parsing blast output
        blast_results, evalue = parse_blast(line, blast_results, cutoff_evalue)

    if blast_results == {}:
        return 0,0
    else:
        candidate_regions, number_regions = merge(blast_results, intron_length)

        return candidate_regions, number_regions

def extract_seq(region_dic, path, tmp_path, mode):

    for key in region_dic:
        #print("blastdbcmd -db " + path + " -dbtype 'nucl' -entry " + key + " -out tmp/" + key + ".fasta -outfmt %f")
        cmd = "blastdbcmd -db " + path + " -dbtype 'nucl' -entry " + key + " -out " + tmp_path + key + ".fasta -outfmt %f"
        starting_subprocess(cmd, mode)

def augustus_ppx(regions, candidatesOutFile, length_extension, profile_path, augustus_ref_species, ass_name, group, tmp_path, mode):
    output = open(candidatesOutFile, "w")

    for key in regions:
        locations = regions[key]
        counter = 0
        for i in locations:
            # some variables
            counter += 1
            start = str(i[0] - length_extension)
            end = str(i[1] + length_extension)
            name = key + "_" + str(counter)
            # augutus call
            cmd = "augustus --protein=1 --proteinprofile=" + profile_path + " --predictionStart=" + start + " --predictionEnd=" + end + " --species=" + augustus_ref_species + " " + tmp_path + key + ".fasta > " + tmp_path + name + ".gff"
            #result = subprocess.run(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell=True)
            starting_subprocess(cmd, mode)
            # transfer augustus output to as sequence
            cmd = "getAnnoFasta.pl --seqfile=" + tmp_path + key + ".fasta " + tmp_path + name + ".gff"
            starting_subprocess(cmd, mode)
            # parsing header and sequences
            try:
                sequence_file = open(tmp_path + name + ".aa", "r")
                lines = sequence_file.readlines()
                for line in lines:
                    if line[0] == ">":
                        id = line.replace(">", "")
                        header = ">" + group + "|" + ass_name + "|" + name + "_" + id
                        output.write(header)
                    else:
                        output.write(line)
                sequence_file.close()
            except FileNotFoundError:
                print("No gene found in region with ID:" + name + " , continuing with next region")
    output.close()

def searching_for_db(assembly_path):

    db_endings = ['.ndb', '.nhr', '.nin', '.nog', '.nos', '.not', '.nsq', '.ntf', '.nto']
    check = True
    for end in db_endings:
        check = check and os.path.exists(assembly_path + end)
    return check

def get_distance_biopython(file, matrix):
    aln = AlignIO.read(open(file), 'fasta')
    calculator = DistanceCalculator(matrix)
    dm = calculator.get_distance(aln)
    return dm

def readFasta(candidatesOutFile):
    seq_records = SeqIO.parse(candidatesOutFile, "fasta")
    return seq_records

def getSeedInfo(path):
    dic = {}
    seq_records = readFasta(path)
    for entry in seq_records:
        species = entry.id.split("|")[1]
        geneID = entry.id.split("|")[2]

        try:
            dic[species].append(geneID)
        except KeyError:
            dic[species] = [geneID]

    del seq_records
    return dic

def checkCoOrthologs(candidate_name, best_hit, ref, fdog_ref_species, candidatesOutFile, msaTool, matrix, dataPath, tmp_path):
    ###########getting sequences and write all in one file to make msa #########
    name_file = candidate_name + ".co"
    output_file = tmp_path + name_file + '.fasta'
    aln_file = tmp_path + name_file + '.aln'
    genome_dir_path = dataPath + '/genome_dir/%s/%s.fa'%(fdog_ref_species, fdog_ref_species)
    #print(searchTool)

    out = open(output_file, "w")
    inSeq = SeqIO.to_dict((SeqIO.parse(open(genome_dir_path), 'fasta')))
    out.write(">" + best_hit + "\n")
    out.write(str(inSeq[best_hit].seq) + "\n")
    out.write(">" + ref + "\n")
    out.write(str(inSeq[ref].seq )+ "\n")

    candidates = readFasta(candidatesOutFile)
    for record in candidates:
        if candidate_name in record.id:
            out.write(">" + candidate_name + "\n")
            out.write(str(record.seq) + "\n")
            break

    out.close()

    if msaTool == "muscle":
        os.system("muscle -quiet -in " + output_file + " -out " + aln_file)
        #print("muscle -quiet -in " + output_file + " -out " + aln_file)
    elif msaTool == "mafft-linsi":
        #print("mafft-linsi")
        os.system('mafft --maxiterate 1000 --localpair --anysymbol --quiet ' + output_file + ' > ' + aln_file)

    distances = get_distance_biopython(aln_file, matrix)

    distance_hit_query = distances[best_hit, candidate_name]
    distance_ref_hit = distances[best_hit, ref]

    if distance_ref_hit < distance_hit_query:
        #accepted
        return 1, distance_ref_hit, distance_hit_query

    else:
        #rejected
        return 0, distance_ref_hit, distance_hit_query

def backward_search(candidatesOutFile, fasta_path, strict, fdog_ref_species, evalue_cut_off, taxa, searchTool, checkCo, msaTool, matrix, dataPath, filter, tmp_path, mode):
    # the backward search uses the genes predicted from augustus and makes a blastp search
    #the blastp search is against all species that are part of the core_ortholog group if the option --strict was chosen or only against the ref taxa
    seedDic = getSeedInfo(fasta_path)
    #print(fasta_path)
    orthologs = []
    #print(seedDic)
    blast_dir_path = dataPath + "/blast_dir/"
    if strict != True:
        seed = [fdog_ref_species]
        try:
            id_ref = seedDic[fdog_ref_species]
        except KeyError:
            print("The fDOG reference species isn't part of the core ortholog group, ... exciting")
            return 0, seed
        if searchTool == "blast":
            cmd = "blastp -db " + blast_dir_path + fdog_ref_species + "/" + fdog_ref_species + " -outfmt '6 sseqid qseqid evalue' -max_target_seqs 10 -out " + tmp_path + "blast_" + fdog_ref_species + " -evalue " + str(evalue_cut_off) + " -query " + candidatesOutFile
            starting_subprocess(cmd, mode)
        else:
            print("diamonds are the girls best friends")
            ##### diamond call

        alg_file = open(tmp_path + "blast_" + fdog_ref_species, "r")
        lines = alg_file.readlines()
        alg_file.close()
        old_name = None
        min = 10
        for line in lines:
            id, gene, evalue = (line.replace("\n", "")).split("\t")
            gene_name = gene.split("|")[2]
            if gene_name != old_name:
                print("candidate:%s"%(gene_name))
                print("blast-hit:%s"%(id))
                min = float(evalue)
                if id in id_ref:
                    orthologs.append(gene)
                    print("\thitting\n")
                else:
                    if checkCo == True:
                        for i in id_ref:
                            print("Best hit %s differs from reference sequence %s! Doing further checks\n"%(id, i))
                            co_orthologs_result, distance_ref_hit, distance_hit_query = checkCoOrthologs(gene_name, id, i, fdog_ref_species, candidatesOutFile, msaTool, matrix, dataPath, tmp_path)
                            if co_orthologs_result == 1:
                                print("\t Distance query - blast hit: %6.4f, Distance blast hit - reference: %6.4f\tAccepting\n"%(distance_hit_query, distance_ref_hit))
                                orthologs.append(gene)
                            elif co_orthologs_result == 0:
                                print("\t Distance query - blast hit: %6.4f, Distance blast hit - reference: %6.4f\tRejecting\n"%(distance_hit_query, distance_ref_hit))
                    else:
                        print("\tnothitting\n")
            elif (gene_name == old_name) and float(evalue) == min and gene_name not in orthologs:
                if id in id_ref:
                    orthologs.append(gene)
                    print("\thitting\n")
                else:
                    if checkCo == True:
                        for i in id_ref:
                            print("Best hit %s differs from reference sequence %s! Doing further checks\n"%(id, i))
                            co_orthologs_result, distance_ref_hit, distance_hit_query = checkCoOrthologs(gene_name, id, i, fdog_ref_species, candidatesOutFile, msaTool, matrix, dataPath, tmp_path)
                            if co_orthologs_result == 1:
                                print("\t Distance query - blast hit: %6.4f, Distance blast hit - reference: %6.4f\tAccepting\n"%(distance_hit_query, distance_ref_hit))
                                orthologs.append(gene)
                            elif co_orthologs_result == 0:
                                print("\t Distance query - blast hit: %6.4f, Distance blast hit - reference: %6.4f\tRejecting\n"%(distance_hit_query, distance_ref_hit))
                    else:
                        print("\tnot hitting\n")
            old_name = gene_name


        if orthologs == []:
            print("No hit in the backward search, ...exciting")
            return 0, seed

    else:
        if taxa != []:
            seed = taxa
            try:
                i = seed.index(fdog_ref_species)
                seed.insert(0,seed.pop(i))
            except ValueError:
                seed.insert(0,fdog_ref_species)
            #print(seed)
            #print("with taxa list from user input")

        else:
            seed = []
            for key in seedDic:
                if key == fdog_ref_species:
                    seed.insert(0,key)
                else:
                    seed.append(key)

        orthologs = set({})

        for species in seed:
            print("backward search in species " + species + "\n")
            orthologs_new = set({})
            try:
                id_ref = seedDic[species]
            except KeyError:
                print("The species " + species + " isn't part of the core ortholog group, ... exciting")
                return 0, seed

            cmd = "blastp -db " + blast_dir_path + species + "/" + species + " -outfmt '6 sseqid qseqid evalue' -max_target_seqs 10 -seg " + filter + " -out " + tmp_path + "/blast_" + species + " -evalue " + str(evalue_cut_off) + " -query " + candidatesOutFile
            starting_subprocess(cmd, mode)
            alg_file = open(tmp_path + "/blast_" + species, "r")
            lines = alg_file.readlines()
            alg_file.close()
            old_name = None
            min = 10
            for line in lines:
                id, gene_name, evalue = (line.replace("\n", "")).split("\t")
                if gene_name != old_name:
                    min = float(evalue)
                    if id in id_ref:
                        orthologs_new.add(gene_name)

                elif (gene_name == old_name) and float(evalue) == min:
                    if id in id_ref:
                        orthologs_new.add(gene_name)

            #print(species)
            #print(orthologs_new)
            if species == fdog_ref_species:
                orthologs = orthologs_new
            else:
                orthologs = orthologs & orthologs_new
                if orthologs == {}:
                    print("No ortholog was found with option --strict")
                    return 0, seed



    #print(orthologs)
    orthologs = set(orthologs)
    return list(orthologs), seed

def addSequences(sequenceIds, candidate_fasta, core_fasta, output, name, species_list, refBool, tmp_path):

    output_file = open(output, "a+")
    if refBool == False:
        seq_records_core = readFasta(core_fasta)
        seq_records_core = list(seq_records_core)
        for species in species_list:
            for entry_core in seq_records_core:
                if species in entry_core.id:
                    output_file.write(">" + entry_core.id + "\n")
                    output_file.write(str(entry_core.seq) + "\n")

    if sequenceIds != 0:
        seq_records_candidate = readFasta(candidate_fasta)
        seq_records_candidate = list(seq_records_candidate)
        for entry_candidate in seq_records_candidate:
            if entry_candidate.id in sequenceIds:
                if entry_candidate.id == sequenceIds[0]:
                    print(entry_candidate.id)
                    output_file.write(">" + entry_candidate.id + "|1" + "\n")
                    output_file.write(str(entry_candidate.seq) + "\n")
                else:
                    output_file.write(">" + entry_candidate.id + "|0" + "\n")
                    output_file.write(str(entry_candidate.seq) + "\n")
    output_file.close()
    return 0

def createFasInput(orthologsOutFile, mappingFile):
    with open(orthologsOutFile, "r") as f:
        fas_seed_id = (f.readline())[1:-1]
        #fas_seed_id = fas_seed_id.split("|")[0]

    mappingFile = open(mappingFile, "a+")

    seq_records = readFasta(orthologsOutFile)
    for seq in seq_records:
        ncbi_id = (seq.id.split("@"))[1]
        mappingFile.write(seq.id + "\t" + "ncbi" + ncbi_id + "\n")


    return fas_seed_id

def cleanup(tmp, tmp_path):
    if tmp == False:
        os.system('rm -r ' + tmp_path)

def checkOptions():
    pass
    #muss ich unbedingt noch ergänzen wenn ich alle möglichen input Optionen implementiert habe!!!

def coorthologs(candidate_names, tmp_path, candidatesFile, fasta, fdog_ref_species, msaTool, matrix):
    if len(candidate_names) == 1:
        return candidate_names

    candidates = readFasta(candidatesFile)
    ref = readFasta(fasta)

    out = tmp_path + '/checkCoorthologs.fa'
    f = open(out,"w")

    aln_file = tmp_path + '/checkCoorthologs.aln'

    for record in ref:
        if fdog_ref_species in record.id:
            ref_id = record.id
            f.write(">" + record.id + "\n")
            f.write(str(record.seq) +  "\n")
            break

    for record in candidates:
        for name in candidate_names:
            if name in record.id:
                f.write(">" + name + "\n")
                f.write(str(record.seq) + "\n")
    f.close()

    if msaTool == "muscle":
        os.system("muscle -quiet -in " + out + " -out " + aln_file)
    elif msaTool == "mafft-linsi":
        os.system('mafft --maxiterate 1000 --localpair --anysymbol --quiet ' + out + ' > ' + aln_file)

    distances = get_distance_biopython(aln_file, matrix)

    min_dist = 10
    min_name = None

    for name in candidate_names:
        distance = distances[ref_id , name]
        if distance <= min_dist:
            min_dist = distance
            min_name = name

    checked = [min_name]

    for name in candidate_names:
        if name == min_name:
            pass
        elif distances[min_name , name] <= distances[min_name , ref_id]:
            checked.append(name)

    return checked

def changes_for_fas(file, header, mode):
    #def replace_first_line( src_filename, target_filename, replacement_line):
    f_in = open(file)
    first_line, remainder = f.readline(), f.read()
    line = first_line.split("|")[0]
    f_in.close()
    f_out = open(file + "s","w")
    f_out.write(line + "\n")
    f_out.write(remainder)
    f_out.close()

class Logger(object):
    def __init__(self, file):
        self.file = file
        self.terminal = sys.stdout
        self.log = self.file

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        pass


def main():

    #################### handle user input ########################################

    version = '0.0.1'

    parser = argparse.ArgumentParser(description='You are running fdog.assembly version ' + str(version) + '.')
    parser.add_argument('--version', action='version', version=str(version))

    required = parser.add_argument_group('Required arguments')
    required.add_argument('--gene', help='Core_ortholog group name. Folder inlcuding the fasta file, hmm file and aln file has to be located in core_orthologs/',
                            action='store', default='', required=True)
    required.add_argument('--augustusRefSpec', help='augustus reference species', action='store', default='', required=True)
    required.add_argument('--refSpec', help='Reference taxon for fDOG.', action='store', default='', required=True)

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--avIntron', help='average intron length of the assembly species in bp (default: 5000)',action='store', default=5000, type=int)
    optional.add_argument('--lengthExtension', help='length extension of the candidate regions in bp (default:5000)', action='store', default=5000, type=int)
    optional.add_argument('--assemblyPath', help='Path for the assembly directory', action='store', default='')
    optional.add_argument('--tmp', help='tmp files will not be deleted', action='store_true', default = False)
    optional.add_argument('--out', help='Output directory', action='store', default='')
    optional.add_argument('--dataPath', help='data directory', action='store', default='')
    optional.add_argument('--coregroupPath', help='core_ortholog directory', action='store', default='')
    optional.add_argument('--searchTool', help='Choose between blast and diamond as alignemnt search tool(default:blast)', action='store', choices=['blast', 'diamond'], default='blast')
    optional.add_argument('--evalBlast', help='E-value cut-off for the Blast search. (default: 0.00001)', action='store', default=0.00001, type=float)
    optional.add_argument('--strict', help='An ortholog is only then accepted when the reciprocity is fulfilled for each sequence in the core set', action='store_true', default=False)
    optional.add_argument('--msaTool', help='Choose between mafft-linsi or muscle for the multiple sequence alignment. DEFAULT: muscle', choices=['mafft-linsi', 'muscle'], action='store', default='muscle')
    optional.add_argument('--checkCoorthologsRef', help='During the final ortholog search, accept an ortholog also when its best hit in the reverse search is not the core ortholog itself, but a co-ortholog of it', action='store_true', default=False)
    optional.add_argument('--scoringmatrix', help='Choose a scoring matrix for the distance criteria used by the option --checkCoorthologsRef. DEFAULT: blosum62', choices=['identity', 'blastn', 'trans', 'benner6', 'benner22', 'benner74', 'blosum100', 'blosum30', 'blosum35', 'blosum40', 'blosum45', 'blosum50', 'blosum55', 'blosum60', 'blosum62', 'blosum65', 'blosum70', 'blosum75', 'blosum80', 'blosum85', 'blosum90', 'blosum95', 'feng', 'fitch', 'genetic', 'gonnet', 'grant', 'ident', 'johnson', 'levin', 'mclach', 'miyata', 'nwsgappep', 'pam120', 'pam180', 'pam250', 'pam30', 'pam300', 'pam60', 'pam90', 'rao', 'risler', 'structure'], action='store', default='blosum62')
    optional.add_argument('--coreTaxa', help='List of core taxa used during --strict', action='store', default='')
    optional.add_argument('--filter', help='Switch the low complexity filter for the blast search on.', action='store', default='no')
    optional.add_argument('--fasoff', help='Turn OFF FAS support', action='store_true', default=False)
    optional.add_argument('--pathFile', help='Config file contains paths to data folder (in yaml format)', action='store', default='')
    optional.add_argument('--searchTaxon', help='Search Taxon name', action='store', default='')
    optional.add_argument('--silent', help='Output will only be written into the log file', action='store_true', default=False)
    optional.add_argument('--debug', help='Stdout and Stderr from fdog.assembly and every used tool will be printed', action='store_true', default=False)


    args = parser.parse_args()

    # required
    group = args.gene
    augustus_ref_species = args.augustusRefSpec
    fdog_ref_species = args.refSpec
    #paths user input
    assemblyDir = args.assemblyPath
    dataPath = args.dataPath
    core_path = args.coregroupPath
    out = args.out
    pathFile = args.pathFile
    #I/O
    tmp = args.tmp
    strict = args.strict
    checkCoorthologs = args.checkCoorthologsRef
    filter = args.filter
    if filter == True or filter == 'yes':
        filter = 'yes'
    else:
        filter = 'no'
    #others
    average_intron_length = args.avIntron
    length_extension = args.lengthExtension
    searchTool = args.searchTool
    evalue = args.evalBlast
    msaTool = args.msaTool
    matrix = args.scoringmatrix
    taxa = args.coreTaxa
    if taxa == '':
        taxa =[]
    else:
        taxa = taxa.split(",")
    fasoff = args.fasoff
    searchTaxon = args.searchTaxon
    silent = args.silent
    debug = args.debug

    if debug == True and silent == True:
        print("It's not possible to use booth modes, please restart and use --debug or --silent")
        return 1
    else:
        if debug == True:
            mode = 'debug'
        elif silent == True:
            mode = 'silent'
        else:
            mode = 'normal'

    #checking paths
    if dataPath == '':
        fdogPath = os.path.realpath(__file__).replace('/fDOGassembly.py','')
        configFile = fdogPath + '/bin/pathconfig.txt'
        if not os.path.exists(configFile):
            sys.exit('No pathconfig.txt found. Please run fdog.setup (https://github.com/BIONF/fDOG/wiki/Installation#setup-fdog) or give a dataPath')
        if pathFile == '':
            with open(configFile) as f:
                dataPath = f.readline().strip()
        else:
            cfg = load_config(pathFile)
            try:
                dataPath = cfg['dataPath']
            except:
                dataPath = 'config'
    if core_path == '':
        core_path = out + '/core_orthologs/'
    else:
        if not core_path.endswith('/'):
            core_path = core_path + '/'

    if assemblyDir == '':
        assemblyDir = dataPath + '/assembly_dir/'
    if out == '':
        #print('test out \n')
        out = os.getcwd()
        os.system('mkdir ' + out + '/' + group + ' >/dev/null 2>&1')
        out = out + '/' + group + '/'

    try:
        f = open(out + "/fdog.log", "a+")
    except FileNotFoundError:
        f = open(out + "/fdog.log", "w")

    ################## How to handle std output and std error ##################

    if mode == 'silent':
        sys.stderr = f
        sys.stdout = f
    else:
        sys.stdout = Logger(f)

    # user input has to be checked here before fDOGassembly continues
    assembly_names = os.listdir(assemblyDir)

    ########################## some variables ##################################

    refBool = False # checks if sequences of reference species were already part of the extended.fa file

    ########### paths ###########

    msa_path = core_path + "/" + group +"/"+ group + ".aln"
    hmm_path = core_path + "/" + group +"/hmm_dir/"+ group + ".hmm"
    fasta_path = core_path + "/" + group +"/"+ group + ".fa"
    consensus_path = out + "/tmp/" + group + ".con"
    profile_path = out + "/tmp/" + group + ".prfl"

    ###################### create tmp folder ###################################

    cmd = 'mkdir ' + out + '/tmp'
    starting_subprocess(cmd, 'silent')

    ######################## consensus sequence ################################

    #make a majority-rule consensus sequence with the tool hmmemit from hmmer
    print("Building a consensus sequence for gene " + group + " \n")
    cmd = 'hmmemit -c -o' + consensus_path + ' ' + hmm_path
    starting_subprocess(cmd, mode)
    print("consensus sequence is finished\n")

    ######################## block profile #####################################

    print("Building a block profile for gene " + group + " \n")
    cmd = 'msa2prfl.pl ' + msa_path + ' --setname=' + group + ' >' + profile_path
    starting_subprocess(cmd, mode)

    if int(os.path.getsize(profile_path)) > 0:
        print("block profile is finished \n")
    else:
        print("Building block profiles failed. Using prepareAlign to convert alignment\n")
        new_path = core_path + group +"/"+ group + "_new.aln"
        #print(cmd)
        cmd = 'prepareAlign < ' + msa_path + ' > ' + new_path
        starting_subprocess(cmd, mode)
        cmd = 'msa2prfl.pl ' + new_path + ' --setname=' + group + ' >' + profile_path
        #print(cmd)
        starting_subprocess(cmd, mode)
        print("block profile is finished \n")

    searchBool = False

    #################### fDOG assembly computation for all species #############
    for asName in assembly_names:
        if searchBool == True:
            break
        if searchTaxon != '' and searchBool == False:
            asName = searchTaxon
            searchBool = True

        ################### path definitions ###################################

        cmd = 'mkdir ' + out + '/tmp/' + asName
        starting_subprocess(cmd, 'silent')
        tmp_path = out + "/tmp/" + asName + "/"
        candidatesOutFile = tmp_path + group + ".candidates.fa"
        if searchTaxon != '':
            orthologsOutFile = out + "/" + group + "_" + asName + ".extended.fa"
            fasOutFile = out + "/" + group + "_" + asName
            mappingFile = tmp_path + group + "_" + asName + ".mapping.txt"
        else:
            orthologsOutFile = out + "/" + group + ".extended.fa"
            fasOutFile = out + "/" + group
            mappingFile = out + "/tmp/" + group + ".mapping.txt"

        print("Searching in species " + asName + "\n")
        assembly_path = assemblyDir + "/" + asName + "/" + asName + ".fa"
        db_path = assemblyDir + "/" + asName + "/blast_dir/" + asName + ".fa"

    ######################## tBLASTn ###########################################
        #checks if data base exists already
        db_check = searching_for_db(db_path)
        if db_check == 0:
            print("creating a blast data base \n")
            cmd = 'makeblastdb -in ' + assembly_path + ' -dbtype nucl -parse_seqids -out ' + db_path
            starting_subprocess(cmd, mode)
            print("database is finished \n")
        else:
            print('blast data base exists already, continuing...')

        #makes a tBLASTn search against the new database
        #codon table argument [-db_gencode int_value], table available ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt
        print("tBLASTn search against data base")
        cmd = 'tblastn -db ' + db_path + ' -query ' + consensus_path + ' -outfmt "6 sseqid sstart send evalue qstart qend " -evalue ' + str(evalue) + ' -out ' + tmp_path + '/blast_results.out'
        starting_subprocess(cmd, mode)
        print("tBLASTn search is finished")

    ################### search for candidate regions and extract seq ###########
    # parse blast and filter for candiate regions
        regions, number_regions = candidate_regions(average_intron_length, evalue, tmp_path)

        if regions == 0:
            #no candidat region are available, no ortholog can be found
            if refBool == True:
                print("No candidate region found")
                continue
        else:
            print(str(number_regions) + " candiate regions were found. Extracting sequences...")
            extract_seq(regions, db_path, tmp_path, mode)

    ############### make Augustus PPX search ###################################

            print("starting augustus ppx \n")
            augustus_ppx(regions, candidatesOutFile, length_extension, profile_path, augustus_ref_species, asName, group, tmp_path, mode)
            print("augustus is finished \n")

    ################# backward search to filter for orthologs###################
            if int(os.path.getsize(candidatesOutFile)) <= 0:
                print("No genes found at candidate regions\n")
                if searchTaxon == '' and refBool == True:
                    continue
                else:
                    reciprocal_sequences = 0
            else:
                reciprocal_sequences, taxa = backward_search(candidatesOutFile, fasta_path, strict, fdog_ref_species, evalue, taxa, searchTool, checkCoorthologs, msaTool, matrix, dataPath, filter, tmp_path, mode)


    ################## checking accepted genes for co-orthologs ################
        if reciprocal_sequences == 0:
            print("No ortholog fulfilled the reciprocity criteria")
            if searchTaxon == '' and refBool == True:
                continue
            else:
                reciprocal_sequences = 0
        else:
            if regions != 0
                reciprocal_sequences = coorthologs(reciprocal_sequences, tmp_path, candidatesOutFile, fasta_path, fdog_ref_species, msaTool, matrix)
            else:
                reciprocal_sequences = 0

    ################ add sequences to extended.fa in the output folder##########

        addSequences(reciprocal_sequences, candidatesOutFile, fasta_path, orthologsOutFile, group, taxa, refBool, tmp_path)
        refBool = True

    ############### make Annotation with FAS ###################################
        # if we want to search in only one Taxon
        if searchTaxon != '' and fasoff == False:
            print("Calculating FAS scores")
            fas_seed_id = createFasInput(orthologsOutFile, mappingFile)
            # bug in calcFAS when using --tsv, have to wait till it's fixed before I can use the option
            cmd = 'mkdir ' + tmp_path + 'anno_dir'
            starting_subprocess(cmd, 'silent')
            cmd = 'calcFAS --seed ' + fasta_path + ' --query ' + orthologsOutFile + ' --annotation_dir ' + tmp_path + 'anno_dir --bidirectional --phyloprofile ' + mappingFile + ' --seed_id "' + fas_seed_id + '" --out_dir ' + out + ' --out_name ' + group + '_' + asName
            starting_subprocess(cmd, mode)
    #if we searched in more than one Taxon and no ortholog was found

    if refBool == False and searchTaxon == '':
        print("No orthologs found. Exciting ...")
        cleanup(tmp, tmp_path)
        return 1
    #if we searched in more than one taxon
    if fasoff == False and searchTaxon == '':
        print("Calculating FAS scores")
        tmp_path = out + '/tmp/'
        fas_seed_id = createFasInput(orthologsOutFile, mappingFile)
        # bug in calcFAS when using --tsv, have to wait till it's fixed before I can use the option
        cmd = 'calcFAS --seed ' + fasta_path + ' --query ' + orthologsOutFile + ' --annotation_dir ' + tmp_path + 'anno_dir --bidirectional --phyloprofile ' + mappingFile + ' --seed_id "' + fas_seed_id + '" --out_dir ' + out + ' --out_name ' + group
        starting_subprocess(cmd, mode)
    ################# remove tmp folder ########################################
    if searchTaxon != '':
        cleanup(tmp, tmp_path)
    else:
        cleanup(tmp, out + "/tmp/")

    f.close()

if __name__ == '__main__':
    main()
