# -*- coding: utf-8 -*-

#######################################################################
# Copyright (C) 2021 Hannah Muelbaier
#
#  This script is used to run fDOG-Assembly which performs targeted ortholog
#  searches on genome assemblies
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
import time
import shutil
import multiprocessing as mp

########################### functions ##########################################
def check_path(path):
    if not os.path.exists(path):
        print(path + " does not exist. Exciting ...")
        sys.exit()

def check_ref_sepc(species_list, fasta_file):
    file = open(fasta_file, "r")
    lines = file.readlines()
    species_file = []

    for line in lines:
        if line[0] == ">":
            species = line.split("|")[1]
            species_file.append(species)
    for species in species_list:
        if species in species_file:
            return species
    print("Reference species is not part of the ortholog group. Exciting ...")
    sys.exit()

def load_config(config_file):
    with open(config_file, 'r') as stream:
        try:
            return yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

def starting_subprocess(cmd, mode, time_out = None):

    try:
        if mode == 'debug':
            result = subprocess.run(cmd, shell=True, timeout = time_out)
        elif mode == 'silent':
            result = subprocess.run(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell=True, timeout = time_out)
        elif mode == 'normal':
            result = subprocess.run(cmd, stdout = subprocess.PIPE, shell=True, timeout = time_out)
    except subprocess.TimeoutExpired:
        return 1

def merge(blast_results, insert_length):
    #merging overlapping and contigous candidate regions
    #format dictionary: {node_name: [(<start>,<send>,evalue, <qstart>,<qend>,<strand>, <score>)]}
    number_regions = 0
    insert_length = int(insert_length)
    score_list = []
    for key in blast_results:
        locations = blast_results[key]
        locations = sorted(locations, key = lambda x: int(x[3]))
        size_list = len(locations)
        j = 0
        while j < size_list-1:
            i = j + 1
            while i < size_list:
                if ((locations[j][0] < locations[i][0]) and (locations[j][1] > locations[i][0]) and (locations[j][5] == locations[i][5]) and (locations[i][5] == '+')):
                    #merge overlapping regions plus strand
                    locations[j][1] = max(locations[j][1], locations[i][1])
                    locations[j][2] = min(locations[j][2], locations[i][2])
                    locations[j][4] = max(locations[j][4], locations[i][4])
                    locations[j][6] = max(locations[j][6], locations[i][6])
                    locations.pop(i)
                    size_list -= 1
                    i -= 1
                elif ((locations[j][1] > locations[i][1]) and (locations[j][0] < locations[i][1]) and (locations[j][5] == locations[i][5]) and (locations[i][5] == '-')):
                    #merge overlapping regions minus strand
                    locations[j][0] = min(locations[j][0], locations[i][0])
                    locations[j][2] = min(locations[j][2], locations[i][2])
                    locations[j][4] = max(locations[j][4], locations[i][4])
                    locations[j][6] = max(locations[j][6], locations[i][6])
                    locations.pop(i)
                    size_list -= 1
                    i -= 1
                elif ((locations[j][0] < locations[i][0]) and (locations[i][0] - locations[j][1] <= 2*insert_length) and (locations[j][5] == locations[i][5]) and (locations[i][5] == '+')):
                    #merging consecutive regions, the distance between booth is not longer than a cutoff, plus strand
                    locations[j][1] = max(locations[j][1], locations[i][1])
                    locations[j][2] = min(locations[j][2], locations[i][2])
                    locations[j][4] = max(locations[j][4], locations[i][4])
                    locations[j][6] = max(locations[j][6], locations[i][6])
                    locations.pop(i)
                    size_list -= 1
                    i -=1
                elif ((locations[j][1] > locations[i][1]) and (locations[j][0] - locations[i][1] <= 2* insert_length) and (locations[j][5] == locations[i][5]) and (locations[i][5] == '-')):
                    #merging consecutive regions, the distance between booth is not longer than a cutoff, minus strand
                    locations[j][0] = min(locations[j][0], locations[i][0])
                    locations[j][2] = min(locations[j][2], locations[i][2])
                    locations[j][4] = max(locations[j][4], locations[i][4])
                    locations[j][6] = max(locations[j][6], locations[i][6])
                    locations.pop(i)
                    size_list -= 1
                    i -=1
                i += 1
            j += 1

        for entry in locations:
            score_list.append(entry[6])
        number_regions += len(locations)
        blast_results[key] = locations

    return blast_results, number_regions, score_list

def parse_blast(line, blast_results, cutoff):
    # format blast line:  <contig> <sstart> <send> <evalue> <qstart> <qend> <score>
    # format dictionary: {node_name: [(<start>,<send>,evalue, <qstart>,<qend>,<strand>, <score>)]}
    line = line.replace("\n", "")
    line_info = line.split("\t")
    evalue = float(line_info[3])
    #cut off
    if evalue > cutoff:
        return blast_results, evalue
    #add region to dictionary
    else:
        node_name, sstart, send, qstart, qend, score = line_info[0], int(line_info[1]), int(line_info[2]), int(line_info[4]), int(line_info[5]), int(line_info[6])
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
            list.append([int(sstart),int(send), evalue, int(qstart), int(qend), strand, score])
            blast_results[node_name] = list
        else:
            blast_results[node_name] = [[int(sstart),int(send), evalue, int(qstart), int(qend), strand, score]]

    return blast_results, evalue

def get_x_results(blast_dic, x, score_list):

    new_dic = {}
    score_list.sort(reverse=True)
    min = score_list[x - 1]
    number_regions = 0

    for key in blast_dic:
        key_list = []
        entries = blast_dic[key]
        for i in entries:
            if i[6] >= min:
                key_list.append(i)
        if key_list != []:
            new_dic[key] = key_list
            number_regions += len(key_list)
    return new_dic, number_regions

def candidate_regions(intron_length, cutoff_evalue, tmp_path, x = 10):
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
        blast_file.close()
        return 0,0
    else:
        candidate_regions, number_regions, score_list = merge(blast_results, intron_length)
        blast_file.close()
        if number_regions > x:
            candidate_regions, number_regions = get_x_results(candidate_regions, x, score_list)
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
            #print(cmd)
            starting_subprocess(cmd, 'silent')
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
        if not os.path.exists(aln_file):
            print("Muscle failed for " + candidate_name + ". Making MSA with Mafft-linsi.")
            os.system('mafft --maxiterate 1000 --localpair --anysymbol --quiet ' + output_file + ' > ' + aln_file)

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
            #print("The fDOG reference species isn't part of the core ortholog group, ... exciting")
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
                print("candidate:%s"%(gene_name)) if mode == "debug" else ""
                print("blast-hit:%s"%(id)) if mode == "debug" else ""
                min = float(evalue)
                if id in id_ref:
                    orthologs.append(gene)
                    print("\thitting\n") if mode == "debug" else ""
                else:
                    if checkCo == True:
                        for i in id_ref:
                            print("Best hit %s differs from reference sequence %s! Doing further checks\n"%(id, i)) if mode == "debug" else ""
                            co_orthologs_result, distance_ref_hit, distance_hit_query = checkCoOrthologs(gene_name, id, i, fdog_ref_species, candidatesOutFile, msaTool, matrix, dataPath, tmp_path)
                            if co_orthologs_result == 1:
                                print("\t Distance query - blast hit: %6.4f, Distance blast hit - reference: %6.4f\tAccepting\n"%(distance_hit_query, distance_ref_hit)) if mode == "debug" else ""
                                orthologs.append(gene)
                            elif co_orthologs_result == 0:
                                print("\t Distance query - blast hit: %6.4f, Distance blast hit - reference: %6.4f\tRejecting\n"%(distance_hit_query, distance_ref_hit)) if mode == "debug" else ""
                    else:
                        print("\tnothitting\n") if mode == "debug" else ""
            elif (gene_name == old_name) and float(evalue) == min and gene_name not in orthologs:
                if id in id_ref:
                    orthologs.append(gene)
                    print("\thitting\n") if mode == "debug" else ""
                else:
                    if checkCo == True:
                        for i in id_ref:
                            print("Best hit %s differs from reference sequence %s! Doing further checks\n"%(id, i)) if mode == "debug" else ""
                            co_orthologs_result, distance_ref_hit, distance_hit_query = checkCoOrthologs(gene_name, id, i, fdog_ref_species, candidatesOutFile, msaTool, matrix, dataPath, tmp_path)
                            if co_orthologs_result == 1:
                                print("\t Distance query - blast hit: %6.4f, Distance blast hit - reference: %6.4f\tAccepting\n"%(distance_hit_query, distance_ref_hit)) if mode == "debug" else ""
                                orthologs.append(gene)
                            elif co_orthologs_result == 0:
                                print("\t Distance query - blast hit: %6.4f, Distance blast hit - reference: %6.4f\tRejecting\n"%(distance_hit_query, distance_ref_hit)) if mode == "debug" else ""
                    else:
                        print("\tnot hitting\n") if mode == "debug" else ""
            old_name = gene_name


        if orthologs == []:
            #print("No hit in the backward search, ...exciting")
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
            print("backward search in species %s\n" %species)
            orthologs_new = set({})
            try:
                id_ref = seedDic[species]
            except KeyError:
                #print("The species " + species + " isn't part of the core ortholog group, ... exciting")
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
            #print(orthologs)
            if species == fdog_ref_species:
                orthologs = orthologs_new
            else:
                orthologs = orthologs & orthologs_new
                if len(orthologs) == 0:
                    #print("No ortholog was found with option --strict")
                    return 0, seed



    #print(orthologs)
    orthologs = set(orthologs)
    return list(orthologs), seed

def addRef(output, core_fasta, species_list):
    #print(species_list)
    output_file = open(output, "a+")
    seq_records_core = readFasta(core_fasta)
    seq_records_core = list(seq_records_core)
    for species in species_list:
        for entry_core in seq_records_core:
            if species in entry_core.id:
                output_file.write(">" + entry_core.id + "\n")
                output_file.write(str(entry_core.seq) + "\n")
    output_file.close()

def addSeq(output, seq_list):
    output_file = open(output, "a+")

    for item in seq_list:
        #print(item)
        candidate_fasta = item[1]
        sequenceIds = item[0]
        if sequenceIds == 0 or sequenceIds == []:
            pass
        seq_records_candidate = readFasta(candidate_fasta)
        seq_records_candidate = list(seq_records_candidate)
        for entry_candidate in seq_records_candidate:
            if entry_candidate.id in sequenceIds:
                if entry_candidate.id == sequenceIds[0]:
                    output_file.write(">" + entry_candidate.id + "|1" + "\n")
                    output_file.write(str(entry_candidate.seq) + "\n")
                else:
                    output_file.write(">" + entry_candidate.id + "|0" + "\n")
                    output_file.write(str(entry_candidate.seq) + "\n")
    output_file.close()

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

    mappingFile.close()
    return fas_seed_id

def cleanup(tmp, tmp_path):
    if tmp == False:
        timeout = time.time() + 60*1
        while os.path.exists(tmp_path):
            shutil.rmtree(tmp_path, ignore_errors=True)
            if time.time() > timeout:
                print("tmp folder could not be removed!")
                break

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

def clean_fas(path, file_type):
    file = open(path, "r")
    lines = file.readlines()
    file.close()
    file = open(path,"w")

    for line in lines:
        if file_type == 'domains':
            long_id, remain = line.split("#")
            id = long_id.split("|")[0]
            new_line = id + "#" + remain
        else:
            long_id, remain = line.split("\t", 1)
            id = long_id.split("|")[0]
            new_line = id + "\t" + remain

        file.write(new_line)
    file.close()

def ortholog_search(args):
    (asName, out, assemblyDir, consensus_path, augustus_ref_species, group, length_extension, average_intron_length, evalue, strict, fdog_ref_species, msaTool, matrix, dataPath, filter, mode, fasta_path, profile_path, taxa, searchTool, checkCoorthologs) = args
    cmd = 'mkdir ' + out + '/tmp/' + asName
    starting_subprocess(cmd, 'silent')
    tmp_path = out + "tmp/" + asName + "/"
    candidatesOutFile = tmp_path + group + ".candidates.fa"
    #orthologsOutFile = out + "/" + group + ".extended.fa"
    fasOutFile = out + "/" + group
    #mappingFile = out + "/tmp/" + group + ".mapping.txt"

    print("Searching in species " + asName + "\n")
    assembly_path = assemblyDir + "/" + asName + "/" + asName + ".fa"
    db_path = assemblyDir + "/" + asName + "/blast_dir/" + asName + ".fa"
    db_check = searching_for_db(db_path)

    if db_check == 0:
        #print("Creating a blast data base...")
        cmd = 'makeblastdb -in ' + assembly_path + ' -dbtype nucl -parse_seqids -out ' + db_path
        starting_subprocess(cmd, mode)
        #print("\t ...finished \n")

    #makes a tBLASTn search against database
    #codon table argument [-db_gencode int_value], table available ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt
    #print("Starting tBLASTn search...")
    cmd = 'tblastn -db ' + db_path + ' -query ' + consensus_path + ' -outfmt "6 sseqid sstart send evalue qstart qend score " -evalue ' + str(evalue) + ' -out ' + tmp_path + '/blast_results.out'
    time_tblastn_start = time.time()
    exit_code = starting_subprocess(cmd, mode, 3600)
    time_tblastn_end = time.time()
    time_tblastn = time_tblastn_end - time_tblastn_start
    if exit_code == 1:
        print("The tblastn search takes too long for species %s. Exciting ..." % asName)
        f.close()
        cleanup(tmp, tmp_folder)
        sys.exit()
    #else:
        #print("\t ...finished")
    print("Time tblastn %s in species %s" % (str(time_tblastn), asName))

    regions, number_regions = candidate_regions(average_intron_length, evalue, tmp_path)
    if regions == 0:
        #no candidat region are available, no ortholog can be found
        print("No candidate region found for species %s!\n" % asName)
        return [], candidatesOutFile

    else:
        print(str(number_regions) + " candiate region(s) were found for species %s.\n" % asName)
        extract_seq(regions, db_path, tmp_path, mode)

    ############### make Augustus PPX search ###################################
    #print("Starting augustus ppx ...")
    time_augustus_start = time.time()
    augustus_ppx(regions, candidatesOutFile, length_extension, profile_path, augustus_ref_species, asName, group, tmp_path, mode)
    #print("\t ...finished \n")
    time_augustus_end = time.time()
    time_augustus = time_augustus_end - time_augustus_start
    print("Time augustus: %s species %s \n" % (str(time_augustus), asName))

    ################# backward search to filter for orthologs###################
    if int(os.path.getsize(candidatesOutFile)) <= 0:
        #print("No genes found at candidate regions\n")
        return [], candidatesOutFile

    reciprocal_sequences, taxa = backward_search(candidatesOutFile, fasta_path, strict, fdog_ref_species, evalue, taxa, searchTool, checkCoorthologs, msaTool, matrix, dataPath, filter, tmp_path, mode)

    if reciprocal_sequences == 0:
        if regions != 0:
            print("No ortholog fulfilled the reciprocity criteria for species %s.\n" % asName)
        return [], candidatesOutFile
    else:
        reciprocal_sequences = coorthologs(reciprocal_sequences, tmp_path, candidatesOutFile, fasta_path, fdog_ref_species, msaTool, matrix)

    return reciprocal_sequences, candidatesOutFile

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

    #################### handle user input #####################################

    start = time.time()

    version = '0.1.2'
    ################### initialize parser ######################################
    parser = argparse.ArgumentParser(description='You are running fdog.assembly version ' + str(version) + '.')
    parser.add_argument('--version', action='version', version=str(version))
    ################## required arguments ######################################
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--gene', help='Core_ortholog group name. Folder inlcuding the fasta file, hmm file and aln file has to be located in core_orthologs/',
                            action='store', default='', required=True)
    required.add_argument('--augustusRefSpec', help='augustus reference species', action='store', default='', required=True)
    required.add_argument('--refSpec', help='Reference taxon for fDOG.', action='store', nargs="+", default='', required=True)
    ################## optional arguments ######################################
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--avIntron', help='average intron length of the assembly species in bp (default: 50000)',action='store', default=50000, type=int)
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
    optional.add_argument('--coreTaxa', help='List of core taxa used during --strict', action='store', nargs="+", default=[])
    optional.add_argument('--filter', help='Switch the low complexity filter for the blast search on.', action='store', default='no')
    optional.add_argument('--fasoff', help='Turn OFF FAS support', action='store_true', default=False)
    optional.add_argument('--pathFile', help='Config file contains paths to data folder (in yaml format)', action='store', default='')
    optional.add_argument('--searchTaxa', help='Search Taxon name', action='store', nargs="+", default=[])
    optional.add_argument('--silent', help='Output will only be written into the log file', action='store_true', default=False)
    optional.add_argument('--debug', help='Stdout and Stderr from fdog.assembly and every used tool will be printed', action='store_true', default=False)
    optional.add_argument('--force', help='Overwrite existing output files', action='store_true', default=False)
    optional.add_argument('--append', help='Append the output to existing output files', action='store_true', default=False)
    optional.add_argument('--parallel', help= 'The ortholog search of multiple species will be done in parallel', action='store_true', default=False)
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
    fasoff = args.fasoff
    searchTaxa = args.searchTaxa
    silent = args.silent
    debug = args.debug
    force = args.force
    append = args.append
    parallel = args.parallel

    # output modes
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

    if out == '':
        out = os.getcwd()
    else:
        if out[-1] != "/":
            out = out + "/"
        check_path(out)

    if os.path.exists(out + '/' + group):
        if append != True and force != True:
            print("Output folder for group " + group + " exists already. Please choose --force or --append.")
            sys.exit()
        elif force == True:
            shutil.rmtree(out + '/' + group, ignore_errors=True)
            refBool = False
            os.system('mkdir ' + out + '/' + group + ' >/dev/null 2>&1')
            out = out + '/' + group + '/'
        elif append == True:
            out = out + '/' + group + '/'
            refBool = True
        else:
            refBool = False # checks if sequences of reference species were already part of the extended.fa file
    else:
        os.system('mkdir ' + out + '/' + group + ' >/dev/null 2>&1')
        out = out + '/' + group + '/'
        refBool = False

    if core_path == '':
        core_path = out + '/core_orthologs/'
    else:
        if not core_path.endswith('/'):
            core_path = core_path + '/'
        check_path(core_path)

    if assemblyDir == '':
        assemblyDir = dataPath + '/assembly_dir/'
    check_path(assemblyDir)

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

    ########################### other variables ################################
    if searchTaxa == []:
        assembly_names = os.listdir(assemblyDir)
    else:
        assembly_names = os.listdir(assemblyDir)
        for Taxon in searchTaxa:
            if Taxon not in assembly_names:
                print("Taxon %s is not in the assembly_dir" % Taxon)
                sys.exit()
        assembly_names = searchTaxa

    ################################# paths ####################################

    msa_path = core_path + "/" + group +"/"+ group + ".aln"
    check_path(msa_path)
    hmm_path = core_path + "/" + group +"/hmm_dir/"+ group + ".hmm"
    check_path(hmm_path)
    fasta_path = core_path + "/" + group +"/"+ group + ".fa"
    check_path(fasta_path)
    consensus_path = out + "/tmp/" + group + ".con"
    profile_path = out + "/tmp/" + group + ".prfl"
    tmp_folder = out + "/tmp"

    ########### is/are fDOG reference species part of ortholog group? ##########

    fdog_ref_species = check_ref_sepc(fdog_ref_species, fasta_path)

    ###################### create tmp folder ###################################

    cmd = 'mkdir ' + out + '/tmp'
    starting_subprocess(cmd, 'silent')

    print("Gene: " + group)
    print("fDOG reference species: " + fdog_ref_species + " \n")

    ######################## consensus sequence ################################
    group_computation_time_start = time.time()
    #make a majority-rule consensus sequence with the tool hmmemit from hmmer
    print("Building a consensus sequence")
    cmd = 'hmmemit -c -o' + consensus_path + ' ' + hmm_path
    starting_subprocess(cmd, mode)
    print("\t ...finished\n")

    ######################## block profile #####################################

    print("Building a block profile ...")
    cmd = 'msa2prfl.pl ' + msa_path + ' --setname=' + group + ' >' + profile_path
    starting_subprocess(cmd, 'silent')

    if int(os.path.getsize(profile_path)) > 0:
        print("\t ...finished \n")
    else:
        print("Building block profiles failed. Using prepareAlign to convert alignment\n")
        new_path = core_path + group +"/"+ group + "_new.aln"
        cmd = 'prepareAlign < ' + msa_path + ' > ' + new_path
        starting_subprocess(cmd, mode)
        cmd = 'msa2prfl.pl ' + new_path + ' --setname=' + group + ' >' + profile_path
        starting_subprocess(cmd, 'silent')
        print(" \t ...finished \n")

    group_computation_time_end = time.time()
    time_group = group_computation_time_end - group_computation_time_start

    ###################### ortholog search #####################################

    ortholog_sequences = []
    time_ortholog_start = time.time()
    if parallel == True:
        ##################### parallel compuataion #############################
        calls = []
        cpus = mp.cpu_count()
        pool = mp.Pool(cpus)
        for asName in assembly_names:
            calls.append([asName, out, assemblyDir, consensus_path, augustus_ref_species, group, length_extension, average_intron_length, evalue, strict, fdog_ref_species, msaTool, matrix, dataPath, filter, mode, fasta_path, profile_path, taxa, searchTool, checkCoorthologs])

        results = (pool.imap_unordered(ortholog_search, calls))
        pool.close()
        pool.join()
        for i in results:
            ortholog_sequences.append(i)
    else:
        ###################### computation species per species ################
        for asName in assembly_names:
            args = [asName, out, assemblyDir, consensus_path, augustus_ref_species, group, length_extension, average_intron_length, evalue, strict, fdog_ref_species, msaTool, matrix, dataPath, filter, mode, fasta_path, profile_path, taxa, searchTool, checkCoorthologs]
            reciprocal_sequences, candidatesOutFile = ortholog_search(args)
            ortholog_sequences.append([reciprocal_sequences, candidatesOutFile])

    ################## preparing output ########################################
    orthologsOutFile = out + "/" + group + ".extended.fa"
    time_ortholog_end = time.time()
    time_ortholog = time_ortholog_end - time_ortholog_start
    if taxa == []:
        taxa = [fdog_ref_species]
    if append == True:
        addSeq(orthologsOutFile, ortholog_sequences)
    else:
        addRef(orthologsOutFile, fasta_path, taxa)
        addSeq(orthologsOutFile, ortholog_sequences)
    mappingFile = out + "/tmp/" + group + ".mapping.txt"

    if fasoff == False:
        fas = time.time()
        print("Calculating FAS scores ...")
        tmp_path = out + '/tmp/'
        fas_seed_id = createFasInput(orthologsOutFile, mappingFile)
        cmd = 'fas.run --seed ' + fasta_path + ' --query ' + orthologsOutFile + ' --annotation_dir ' + tmp_path + 'anno_dir --bidirectional --tsv --phyloprofile ' + mappingFile + ' --seed_id "' + fas_seed_id + '" --out_dir ' + out + ' --out_name ' + group
        starting_subprocess(cmd, 'silent')
        clean_fas(out + group + "_forward.domains", 'domains')
        clean_fas(out + group + "_reverse.domains", 'domains')
        clean_fas(out + group + ".phyloprofile", 'phyloprofile')
        print("\t ...finished \n")
    ################# remove tmp folder ########################################
    end = time.time()
    time_fas = end - fas
    print("fDOG-Assembly finished completely in " + str(end-start) + "seconds.")
    print("Group preparation: %s \t Ortholog search: %s \t Fas: %s \n" % (str(time_group), str(time_ortholog), str(time_fas)))
    sys.stdout = sys.__stdout__

    f.close()
    cleanup(tmp, tmp_folder)

if __name__ == '__main__':
    main()
