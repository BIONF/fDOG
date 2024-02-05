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
import fdog.libs.alignment as align_fn
from tqdm import tqdm

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

def extract_sequence_from_to(name, file, start, end):
    #print(name)
    out = name + ".fasta"
    if int(start) < 0:
        start = 0
    with open(out,"w") as f:
        for seq_record in SeqIO.parse(file, "fasta"):
                f.write(">" + str(seq_record.id) + "\n")
                sequence_length = len(seq_record.seq)
                if int(end) > sequence_length:
                    end = sequence_length
                #for testing only
                #start = 0
                #end = len(seq_record.seq)
                f.write(str(seq_record.seq[int(start):int(end)]) + "\n")

    return out, start, end

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
                pass
                #print("No gene found in region with ID" + name + " in species " + ass_name + " , continuing with next region")
    output.close()

def metaeuk_single(regions, candidatesOutFile, length_extension, ass_name, group, tmp_path, mode, db):
    output = open(candidatesOutFile, "w")
    region = open(candidatesOutFile.replace(".candidates.fa", ".regions.txt"), "w")
    region.write("Contig/scaffold" + "\t" + "start" + "\t" + "end" + "\n")

    for key in regions:
        locations = regions[key]
        counter = 0
        for i in locations:
            #some variables
            counter += 1
            start = str(i[0] - length_extension)
            end = str(i[1] + length_extension)
            name = key + "_" + str(counter)
            file, start, end = extract_sequence_from_to(tmp_path + name, tmp_path + key + ".fasta", start, end)
            region.write(file + "\t" + str(start) + "\t" + str(end) + "\n")
            #metaeuk call
            cmd = "metaeuk easy-predict " + file + " " + db + " " + tmp_path + name + " " + tmp_path + "/metaeuk --min-exon-aa 5 --max-overlap 5 --min-intron 1 --overlap 1"
            #print(cmd)
            # other parameteres used by BUSCO with metazoa set--max-intron 130000 --max-seq-len 160000 --min-exon-aa 5 --max-overlap 5 --min-intron 1 --overlap 1
            starting_subprocess(cmd, mode)
            # parsing header and sequences
            try:
                sequence_file = open(tmp_path + name + ".fas", "r")
                lines = sequence_file.readlines()
                #print(lines)
                id = 0
                for line in lines:
                    if line[0] == ">":
                        id += 1
                        header = ">" + group + "|" + ass_name + "|" + name + "_" + str(id) + "\n"
                        output.write(header)
                    else:
                        output.write(line)
                sequence_file.close()

                gff_file = open(tmp_path + name + ".gff", "r")
                lines = gff_file.readlines()
                new_lines = []
                for line in lines:
                    values = line.split("\t")
                    values[3] = str(int(values[3]) + int(start))
                    values[4] = str(int(values[4]) + int(start))
                    new_lines.append("\t".join(values))
                gff_file.close()
                gff_file = open(tmp_path + name + ".gff", "w")
                for line in new_lines:
                    gff_file.write(line)
                gff_file.close()
            except FileNotFoundError:
                pass

    output.close()

def searching_for_db(assembly_path):

    db_endings = ['.ndb', '.nhr', '.nin', '.nog', '.nos', '.not', '.nsq', '.ntf', '.nto']
    check = True
    for end in db_endings:
        if not any(File.endswith(end) for File in os.listdir(assembly_path)):
            check = False
    return check

def get_distance_biopython(file, matrix):
    #print(file)
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

def checkCoOrthologs(candidate_name, best_hit, ref, fdog_ref_species, candidatesOutFile, msaTool, matrix, dataPath, tmp_path, mode='silent'):
    ###########getting sequences and write all in one file to make msa #########
    name_file = candidate_name + ".co"
    output_file = tmp_path + name_file + '.fa'
    aln_file = tmp_path + name_file + '.aln'
    genome_dir_path = dataPath + '/searchTaxa_dir/%s/%s.fa'%(fdog_ref_species, fdog_ref_species)
    if not os.path.exists(genome_dir_path):
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
        if align_fn.get_muscle_version(msaTool) == 'v3':
            cmd = "muscle -quiet -in " + output_file + " -out " + aln_file
        else:
            cmd = "muscle -align " + output_file + " -output " + aln_file
        starting_subprocess(cmd, mode)
        if not os.path.exists(aln_file):
            print("Muscle failed for %s. Making MSA with Mafft-linsi." % (candidate_name))
            cmd = 'mafft --maxiterate 1000 --localpair --anysymbol --quiet ' + output_file + ' > ' + aln_file
            starting_subprocess(cmd, mode)

    elif msaTool == "mafft-linsi":
        #print("mafft-linsi")
        cmd = 'mafft --maxiterate 1000 --localpair --anysymbol --quiet ' + output_file + ' > ' + aln_file
        starting_subprocess(cmd, mode)

    try:
        distances = get_distance_biopython(aln_file, matrix)
        distance_hit_query = distances[best_hit, candidate_name]
        distance_ref_hit = distances[best_hit, ref]
        #print(distances)
    except ValueError:
        pass
        #print("Failure in distance computation, Candidate  %s will be rejected" % candidate_name)
        return 0, "NaN", "NaN"

    #distance_hit_query = distances[best_hit, candidate_name]
    #distance_ref_hit = distances[best_hit, ref]

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
    blast_dir_path = dataPath + "/coreTaxa_dir/"
    #print(blast_dir_path)
    if not os.path.exists(blast_dir_path):
        blast_dir_path = dataPath + "/blast_dir/"
    #print(blast_dir_path)
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
                                if distance_ref_hit != "NaN":
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
            continue
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
    else:
        # clean up whole contigs
        for root, dirs, files in os.walk(tmp_path):
            for file in files:
                if file.endswith(".fasta"):
                    os.remove(os.path.join(root, file))

def coorthologs(candidate_names, tmp_path, candidatesFile, fasta, fdog_ref_species, msaTool, matrix, mode='silent'):
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

    already_written = []
    for record in candidates:
        for name in candidate_names:
            if name == record.id:
                if name not in already_written:
                    f.write(">" + record.id + "\n")
                    f.write(str(record.seq) + "\n")
                    already_written.append(name)
    f.close()

    if msaTool == "muscle":
        if align_fn.get_muscle_version(msaTool) == 'v3':
            cmd = "muscle -quiet -in %s -out %s" % (out, aln_file)
            #print("muscle -quiet -in " + output_file + " -out " + aln_file)
        else:
            cmd = "muscle -align %s -output %s" % (out, aln_file)
        starting_subprocess(cmd, mode)
        if not os.path.exists(aln_file):
            print("Muscle failed for %s. Making MSA with Mafft-linsi." % (aln_file))
            cmd = 'mafft --maxiterate 1000 --localpair --anysymbol --quiet ' + out + ' > ' + aln_file
            starting_subprocess(cmd, mode)
    elif msaTool == "mafft-linsi":
        cmd = 'mafft --maxiterate 1000 --localpair --anysymbol --quiet %s > %s'% (out, aln_file)
        starting_subprocess(cmd, mode)

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

def ortholog_search_tblastn(args):
    (asName, out, assemblyDir, consensus_path, augustus_ref_species, group, length_extension, average_intron_length, evalue, strict, fdog_ref_species, msaTool, matrix, dataPath, filter, mode, fasta_path, profile_path, taxa, searchTool, checkCoorthologs, gene_prediction, metaeuk_db) = args
    output = []
    cmd = 'mkdir ' + out + '/tmp/' + asName
    starting_subprocess(cmd, 'silent')
    tmp_path = out + "tmp/" + asName + "/"
    candidatesOutFile = tmp_path + group + ".candidates.fa"
    #orthologsOutFile = out + "/" + group + ".extended.fa"
    fasOutFile = out + "/" + group
    #mappingFile = out + "/tmp/" + group + ".mapping.txt"

    output.append("Searching in species " + asName + "\n")
    assembly_path = assemblyDir + "/" + asName + "/" + asName + ".fa"
    db_path = assemblyDir + "/" + asName + "/blast_dir/" + asName + ".fa"
    blast_dir_path = assemblyDir + "/" + asName + "/blast_dir/"
    if not os.path.exists(blast_dir_path):
        cmd = 'mkdir ' + blast_dir_path
        starting_subprocess(cmd, 'silent')
    db_check = searching_for_db(blast_dir_path)

    if db_check == 0:
        cmd = 'makeblastdb -in ' + assembly_path + ' -dbtype nucl -parse_seqids -out ' + db_path
        starting_subprocess(cmd, mode)

    #makes a tBLASTn search against database
    #codon table argument [-db_gencode int_value], table available ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt
    cmd = 'tblastn -db ' + db_path + ' -query ' + consensus_path + ' -outfmt "6 sseqid sstart send evalue qstart qend score " -evalue ' + str(evalue) + ' -out ' + tmp_path + '/blast_results.out'
    time_tblastn_start = time.time()
    exit_code = starting_subprocess(cmd, mode, 3600)
    time_tblastn_end = time.time()
    time_tblastn = time_tblastn_end - time_tblastn_start
    if exit_code == 1:
        output.append("The tblastn search takes too long for species %s. Skipping species ..." % asName)
        return [], candidatesOutFile, output

    output.append("Time tblastn %s in species %s" % (str(time_tblastn), asName))

    regions, number_regions = candidate_regions(average_intron_length, evalue, tmp_path)
    if regions == 0:
        #no candidat region are available, no ortholog can be found
        output.append("No candidate region found for species %s!\n" % asName)
        return [], candidatesOutFile, output

    else:
        output.append(str(number_regions) + " candiate region(s) were found for species %s.\n" % asName)
        extract_seq(regions, db_path, tmp_path, mode)

    if gene_prediction == "augustus":
        ############### make Augustus PPX search ###################################
        time_augustus_start = time.time()
        augustus_ppx(regions, candidatesOutFile, length_extension, profile_path, augustus_ref_species, asName, group, tmp_path, mode)
        time_augustus_end = time.time()
        time_augustus = time_augustus_end - time_augustus_start
        output.append("Time augustus: %s species %s \n" % (str(time_augustus), asName))
    else:
        time_metaeuk_start = time.time()
        if metaeuk_db == '':
            db = fasta_path
        else:
            db = metaeuk_db
        metaeuk_single(regions, candidatesOutFile, length_extension, asName, group, tmp_path, mode, db)
        time_metaeuk_end = time.time()
        time_metaeuk = time_metaeuk_end - time_metaeuk_start
        output.append("Time metaeuk: %s species %s \n" % (str(time_metaeuk), asName))

    ################# backward search to filter for orthologs###################
    if int(os.path.getsize(candidatesOutFile)) <= 0:
        #print("No genes found at candidate regions\n")
        return [], candidatesOutFile, output

    reciprocal_sequences, taxa = backward_search(candidatesOutFile, fasta_path, strict, fdog_ref_species, evalue, taxa, searchTool, checkCoorthologs, msaTool, matrix, dataPath, filter, tmp_path, mode)

    if reciprocal_sequences == 0:
        if regions != 0:
            output.append("No ortholog fulfilled the reciprocity criteria for species %s.\n" % asName)
        return [], candidatesOutFile, output
    else:
        reciprocal_sequences = coorthologs(reciprocal_sequences, tmp_path, candidatesOutFile, fasta_path, fdog_ref_species, msaTool, matrix)

    return reciprocal_sequences, candidatesOutFile, output

def blockProfiles(core_path, group, mode, out, msaTool):

    ######################## paths ################################
    msa_path = core_path + "/" + group +"/"+ group + ".aln"
    if not os.path.exists(msa_path):
        fasta_path = core_path + "/" + group +"/"+ group + ".fa"
        check_path(fasta_path)
        if msaTool == "muscle":
            if align_fn.get_muscle_version(msaTool) == 'v3':
                cmd= "muscle -quiet -in " + fasta_path + " -out " + msa_path
                #print("muscle -quiet -in " + output_file + " -out " + aln_file)
            else:
                cmd = "muscle -quiet -align" + fasta_path + " -out " + msa_path
        elif msaTool == "mafft-linsi":
            cmd = 'mafft --maxiterate 1000 --localpair --anysymbol --quiet ' + fasta_path + ' > ' + msa_path
        starting_subprocess(cmd, mode)

    profile_path = out + "/tmp/" + group + ".prfl"

    ######################## block profile #####################################

    print("Building a block profile ...", flush=True)
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
        print(" \t ...finished \n", flush=True)

    return profile_path

def consensusSequence(core_path, group, mode, out):

    ######################## paths ################################
    hmm_path = core_path + "/" + group +"/hmm_dir/"+ group + ".hmm"
    check_path(hmm_path)
    consensus_path = out + "/tmp/" + group + ".con"

    ######################## consensus sequence ################################
    #make a majority-rule consensus sequence with the tool hmmemit from hmmer
    print("Building a consensus sequence")
    cmd = 'hmmemit -c -o' + consensus_path + ' ' + hmm_path
    starting_subprocess(cmd, mode)
    print("\t ...finished\n")

    return consensus_path

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
    version = '0.1.4'
    ################### initialize parser ######################################
    parser = argparse.ArgumentParser(description='You are running fdog.assembly version ' + str(version) + '.')
    parser.add_argument('--version', action='version', version=str(version))
    ################## required arguments ######################################
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--gene', help='Core_ortholog group name. Folder inlcuding the fasta file, hmm file and aln file has to be located in core_orthologs/',
                            action='store', default='', required=True)
    required.add_argument('--refSpec', help='Reference taxon/taxa for fDOG.', action='store', nargs="+", default='', required=True)
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
    #optional.add_argument('--filter', help='Switch the low complexity filter for the blast search on.', action='store', default='no')
    optional.add_argument('--fasoff', help='Turn OFF FAS support', action='store_true', default=False)
    optional.add_argument('--pathFile', help='Config file contains paths to data folder (in yaml format)', action='store', default='')
    optional.add_argument('--searchTaxa', help='List of Taxa to search in', action='store', nargs="+", default=[])
    optional.add_argument('--silent', help='Output will only be written into the log file', action='store_true', default=False)
    optional.add_argument('--debug', help='Stdout and Stderr from fdog.assembly and every used tool will be printed', action='store_true', default=False)
    optional.add_argument('--force', help='Overwrite existing output files', action='store_true', default=False)
    optional.add_argument('--append', help='Append the output to existing output files', action='store_true', default=False)
    optional.add_argument('--parallel', help= 'The ortholog search of multiple species will be done in parallel', action='store_true', default=False)
    optional.add_argument('--augustus', help= 'Gene prediction is done by using the tool Augustus PPX', action='store_true', default=False)
    optional.add_argument('--augustusRefSpec', help='augustus reference species', action='store', default='')
    optional.add_argument('--metaeukDb', help='path to metaeuk reference database', action='store', default='')
    args = parser.parse_args()

    # required
    group = args.gene
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
    augustus_ref_species = args.augustusRefSpec
    metaeuk_db = args.metaeukDb

    #gene prediction tool
    augustus = args.augustus
    if augustus == True:

        if augustus_ref_species == '':
            print("Augustus reference species is required when using Augustus as gene prediction tool")
            return 1
        gene_prediction = "augustus"
    else:
        gene_prediction = "metaeuk"

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
            os.system('mkdir ' + out + '/' + group + ' >/dev/null 2>&1')
            out = out + '/' + group + '/'
        elif append == True:
            out = out + '/' + group + '/'

    else:
        os.system('mkdir ' + out + '/' + group + ' >/dev/null 2>&1')
        out = out + '/' + group + '/'

    if core_path == '':
        core_path = out + '/core_orthologs/'
    else:
        if not core_path.endswith('/'):
            core_path = core_path + '/'
        check_path(core_path)

    if assemblyDir == '':
        assemblyDir = dataPath + '/assembly_dir/'
    check_path(assemblyDir)

    if metaeuk_db != '':
        check_path(metaeuk_db)


    try:
        f = open(out + "/fdog.log", "a+")
    except FileNotFoundError:
        f = open(out + "/fdog.log", "w")

    ################## How to handle std output and std error ##################

    if mode == 'silent':
        sys.stderr = f
        sys.stdout = f
    else:
        pass
        #sys.stdout = Logger(f)

    ########################### other variables ################################
    if searchTaxa == []:
        assembly_names = os.listdir(assemblyDir)
    else:
        if len(searchTaxa) > 1:
            assembly_names = os.listdir(assemblyDir)
            for Taxon in searchTaxa:
                if Taxon not in assembly_names:
                    print("Taxon %s is not in the assembly_dir" % Taxon)
                    sys.exit()
            assembly_names = searchTaxa
        else:
            if searchTaxa[0] in os.listdir(assemblyDir):
                assembly_names = searchTaxa
            elif os.path.isfile(searchTaxa[0]):
                with open(searchTaxa[0]) as file:
                    lines = file.readlines()
                    assembly_names = [line.rstrip() for line in lines]
            else:
                print("Input %s for search Taxa is not in the assembly_dir or an existing file" % searchTaxa[0])

    ################################# paths ####################################

    fasta_path = core_path + "/" + group +"/"+ group + ".fa"
    check_path(fasta_path)
    tmp_folder = out + "/tmp"

    ########### is/are fDOG reference species part of ortholog group? ##########

    fdog_ref_species = check_ref_sepc(fdog_ref_species, fasta_path)

    ###################### create tmp folder ###################################

    cmd = 'mkdir ' + out + '/tmp'
    starting_subprocess(cmd, 'silent')

    print("Gene: " + group, flush=True)
    print("fDOG reference species: " + fdog_ref_species + " \n",flush=True)

    ###################### preparations ########################################

    if augustus == True:
        group_computation_time_start = time.time()
        consensus_path = consensusSequence(core_path, group, mode, out)
        profile_path = blockProfiles(core_path, group, mode, out, msaTool)
        group_computation_time_end = time.time()
        time_group = group_computation_time_end - group_computation_time_start
    else:
        #print("test")
        profile_path = ""
        group_computation_time_start = time.time()
        consensus_path = consensusSequence(core_path, group, mode, out)
        #concatinade core_group sequences if metaeuk should be run without tblastn
        group_computation_time_end = time.time()
        time_group = group_computation_time_end - group_computation_time_start


    ###################### ortholog search #####################################

    ortholog_sequences = []
    time_ortholog_start = time.time()
    
    if parallel == True:
        ##################### parallel computation #############################
        calls = []
        cpus = mp.cpu_count()
        pool = mp.Pool(cpus)
        for asName in assembly_names:
            calls.append([asName, out, assemblyDir, consensus_path, augustus_ref_species, group, length_extension, average_intron_length, evalue, strict, fdog_ref_species, msaTool, matrix, dataPath, filter, mode, fasta_path, profile_path, taxa, searchTool, checkCoorthologs, gene_prediction, metaeuk_db])

        #results = (pool.imap_unordered(ortholog_search_tblastn, calls))
        #pool.close()
        #pool.join()
      
        print("Searching for orthologs ...", flush=True)
        for i in tqdm(pool.imap_unordered(ortholog_search_tblastn, calls),total=len(calls)):
            ortholog_sequences.append([i[0], i[1]])
            if mode == 'debug':
                for k in i[2]:
                    print(k)
        #for i in results:
            #ortholog_sequences.append([i[0], i[1]])
            #for k in i[2]:
                #print(k)
        print("\t ...finished \n", flush=True)
        
    else:
        ###################### computation species wise ################
        for asName in tqdm(assembly_names):
            args = [asName, out, assemblyDir, consensus_path, augustus_ref_species, group, length_extension, average_intron_length, evalue, strict, fdog_ref_species, msaTool, matrix, dataPath, filter, mode, fasta_path, profile_path, taxa, searchTool, checkCoorthologs, gene_prediction, metaeuk_db]
            reciprocal_sequences, candidatesOutFile, output_ortholog_search = ortholog_search_tblastn(args)
            ortholog_sequences.append([reciprocal_sequences, candidatesOutFile])
            if mode == 'debug':
                for k in output_ortholog_search:
                    print(k)

    time_ortholog_end = time.time()
    time_ortholog = time_ortholog_end - time_ortholog_start

    ################## preparing output ########################################
    orthologsOutFile = out + "/" + group + "_og.fa"

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
        #print(cmd)
        starting_subprocess(cmd, 'silent')
        clean_fas(out + group + "_forward.domains", 'domains')
        clean_fas(out + group + "_reverse.domains", 'domains')
        clean_fas(out + group + ".phyloprofile", 'phyloprofile')
        print("\t ...finished \n")
        end = time.time()
        time_fas = end - fas
    else:
        end = time.time()
        time_fas = 0

    ################# remove tmp folder ########################################

    print("fDOG-Assembly finished completely in " + str(end-start) + "seconds.")
    print("Group preparation: %s \t Ortholog search: %s \t FAS: %s \n" % (str(time_group), str(time_ortholog), str(time_fas)))
    sys.stdout = sys.__stdout__

    f.close()
    cleanup(tmp, tmp_folder)

if __name__ == '__main__':
    main()
