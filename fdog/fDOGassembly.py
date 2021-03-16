############################ imports ###########################################
import os
import os.path
import sys
from Bio import SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
########################### functions ##########################################

def merge(blast_results, insert_length):
    number_regions = 0
    for key in blast_results:
        locations = blast_results[key]
        locations = sorted(locations, key = lambda x: int(x[3]))
        #print("test")
        #print(locations)
        size_list = len(locations)

        j = 0

        while j < size_list-1:
            i = 1
            while i < size_list-1:

                if ((locations[j][0] < locations[i][0]) and (locations[j][1] > locations[i][0]) and (locations[j][5] == locations[i][5])):
                    #merge overlapping regions
                    locations[j][1] = max(locations[j][1], locations[i][1])
                    locations[j][2] = min(locations[j][2], locations[i][2])
                    locations.pop(i)
                    size_list -= 1
                    i -= 1
                elif ((locations[j][0] < locations[i][0]) and (locations[i][0] - locations[j][1] <= 2* insert_length) and (locations[j][5] == locations[i][5])):
                    #print(j)
                    locations[j][1] = max(locations[j][1], locations[i][1])
                    locations[j][2] = min(locations[j][2], locations[i][2])
                    locations.pop(i)
                    size_list -= 1
                    i -=1
                i += 1
            j += 1

        number_regions += len(locations)
        blast_results[key] = locations

    #print(blast_results)
    return blast_results, number_regions

def parse_blast(line, blast_results):
    # format blast line:  <contig> <sstart> <send> <evalue> <qstart> <qend> <strand>
    #fomrat dictionary: {node_name: [(<start>,<end>)]}
    #print(line)
    line = line.replace("\n", "")
    line_info = line.split("\t")
    #print(line_info)
    evalue = float(line_info[3])

    #cut off
    if evalue > 0.00001:
        return blast_results, evalue
    #add region to dictionary
    else:
        node_name, sstart, send, qstart, qend = line_info[0], line_info[1], line_info[2], line_info[4], line_info[5]
        split = node_name.split("|")

        # finding out on which strand tBLASTn founded a hit
        if sstart < send:
            strand = "+"
        else:
            sstart = line_info[2]
            send = line_info[1]
            strand = "-"

        #creating a dictionary that inlcudes every tBLASTn that is better as the evalue cut-off of 0.00001
        if len(split) > 1:
            node_name = split[1]
        if node_name in blast_results:
            list = blast_results[node_name]
            list.append([int(sstart),int(send), evalue, int(qstart), int(qend), strand])
            blast_results[node_name] = list
        else:
            blast_results[node_name] = [[int(sstart),int(send), evalue, int(qstart), int(qend), strand]]

    return blast_results, evalue

def candidate_regions(intron_length, evalue):
    ###################### extracting candidate regions ########################
    # info about output blast http://www.metagenomics.wiki/tools/blast/blastn-output-format-6
    blast_file = open("tmp/blast_results.out", "r")
    evalue = 0
    blast_results = {}
    #parsing blast output
    while True:
        line = blast_file.readline()
        #end of file is reached
        if not line:
            break
        #parsing blast output
        blast_results, evalue = parse_blast(line, blast_results)
        #evalue cut-off
        if not evalue <= evalue:
            break
    if blast_results == {}:
        return 0,0
    else:
        candidate_regions, number_regions = merge(blast_results, intron_length)
        #candidate_regions, number_regions = merge_regions(blast_results, cut_off)
        #print(candidate_regions, number_regions)
        return candidate_regions, number_regions

def extract_seq(region_dic, path):
    #print(region_dic)
    for key in region_dic:
        #print("blastdbcmd -db " + path + " -dbtype 'nucl' -entry " + key + " -out tmp/" + key + ".fasta -outfmt %f")
        os.system("blastdbcmd -db " + path + " -dbtype 'nucl' -entry " + key + " -out tmp/" + key + ".fasta -outfmt %f")

def augustus_ppx(regions, candidatesOutFile, length_extension, profile_path, augustus_ref_species, ass_name, group):
    output = open(candidatesOutFile, "w")

    for key in regions:
        locations = regions[key]
        counter = 0
        for i in locations:
            counter += 1
            start = str(i[0] - length_extension)
            end = str(i[1] + length_extension)
            name = key + "_" + str(counter)
            #print("augustus --proteinprofile=" + profile_path + " --predictionStart=" + start + " --predictionEnd=" + end + " --species=" + augustus_ref_species + " tmp/" + key + ".fasta > tmp/" + key + ".gff")
            os.system("augustus --protein=1 --proteinprofile=" + profile_path + " --predictionStart=" + start + " --predictionEnd=" + end + " --species=" + augustus_ref_species + " tmp/" + key + ".fasta > tmp/" + name + ".gff")
            os.system("getAnnoFasta.pl --seqfile=tmp/" + key + ".fasta" + " tmp/" + name + ".gff")

            sequence_file = open("tmp/" + name + ".aa", "r")
            lines = sequence_file.readlines()
            for line in lines:
                if line[0] == ">":
                    id = line.replace(">", "")
                    header = ">" + group + "|" + ass_name + "|" + name + "_" + id
                    output.write(header)
                else:
                    output.write(line)
            sequence_file.close()

    output.close()

def searching_for_db(assembly_path):
    #print("test: " + str(assembly_path) + "\n")
    db_endings = ['.ndb', '.nhr', '.nin', '.nog', '.nos', '.not', '.nsq', '.ntf', '.nto']
    check = True
    for end in db_endings:
        #print(assembly_path + end + "\n")
        check = check and os.path.exists(assembly_path + end)
        #print(check)
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

def checkCoOrthologs(candidate_name, best_hit, ref, fdog_ref_species, candidatesOutFile, aligner, matrix):
    ###########getting sequences and write all in one file to make msa #########
    name_file = candidate_name + ".co"
    output_file = 'tmp/' + name_file + '.fasta'
    aln_file = 'tmp/' + name_file + '.aln'
    genome_dir_path = '../data/genome_dir/%s/%s.fa'%(fdog_ref_species, fdog_ref_species)
    #print(aligner)

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

    if aligner == "muscle":
        os.system("muscle -quiet -in " + output_file + " -out " + aln_file)
        #print("muscle -quiet -in " + output_file + " -out " + aln_file)
    elif aligner == "mafft-linsi":
        #print("mafft-linsi")
        os.system('mafft --maxiterate 1000 --localpair --anysymbol --quiet ' + output_file + ' > ' + aln_file)

    #d_ref = get_distance(aln_file, best_hit, ref)
    #d = get_distance(aln_file, best_hit, candidate_name)
    distances = get_distance_biopython(aln_file, matrix)

    distance_hit_query = distances[best_hit, candidate_name]
    distance_ref_hit = distances[best_hit, ref]

    if distance_ref_hit < distance_hit_query:
        return 1, distance_ref_hit, distance_hit_query
        #accepted
    else:
        return 0, distance_ref_hit, distance_hit_query
        #rejected

def backward_search(candidatesOutFile, fasta_path, strict, fdog_ref_species, evalue_cut_off, taxa, aligner, checkCo, msaTool, matrix):
    # the backward search uses the genes predicted from augustus and makes a blastp search
    #the blastp search is against all species that are part of the core_ortholog group if the option --strict was chosen or only against the ref taxa

    seedDic = getSeedInfo(fasta_path)
    orthologs = []
    #print(seedDic)
    blast_dir_path = "../data/blast_dir/"
    if strict != True:
        seed = [fdog_ref_species]
        try:
            id_ref = seedDic[fdog_ref_species]
        except KeyError:
            print("The fDOG reference species isn't part of the core ortholog group, ... exciting")
            return 0, seed
        if aligner == "blast":
            os.system("blastp -db " + blast_dir_path + fdog_ref_species + "/" + fdog_ref_species + " -outfmt '6 sseqid qseqid evalue' -max_target_seqs 10 -out tmp/blast_" + fdog_ref_species + " -evalue " + str(evalue_cut_off) + " -query " + candidatesOutFile)
        else:
            print("diamonds are the girls best friends")
            ##### diamond call

        alg_file = open("tmp/blast_" + fdog_ref_species, "r")
        lines = alg_file.readlines()
        alg_file.close()
        old_name = None
        min = 10
        for line in lines:
            id, gene_name, evalue = (line.replace("\n", "")).split("\t")
            gene_name = gene_name.split("|")[2]
            if gene_name != old_name:
                print("candidate:%s"%(gene_name))
                print("blast-hit:%s"%(id))
                min = float(evalue)
                if id in id_ref:
                    orthologs.append(gene_name)
                    print("\thitting\n")
                else:
                    if checkCo == True:
                        for i in id_ref:
                            print("Best hit %s differs from reference sequence %s! Doing further checks\n"%(id, i))
                            co_orthologs_result, distance_ref_hit, distance_hit_query = checkCoOrthologs(gene_name, id, i, fdog_ref_species, candidatesOutFile, msaTool, matrix)
                            if co_orthologs_result == 1:
                                print("\t Distance query - blast hit: %6.4f, Distance blast hit - reference: %6.4f\tAccepting\n"%(distance_hit_query, distance_ref_hit))
                                orthologs.append(gene_name)
                            elif co_orthologs_result == 0:
                                print("\t Distance query - blast hit: %6.4f, Distance blast hit - reference: %6.4f\tRejecting\n"%(distance_hit_query, distance_ref_hit))
                    else:
                        print("\tnothitting\n")
            elif (gene_name == old_name) and float(evalue) == min and gene_name not in orthologs:
                if id in id_ref:
                    orthologs.append(gene_name)
                    print("\thitting\n")
                else:
                    if checkCo == True:
                        for i in id_ref:
                            print("Best hit %s differs from reference sequence %s! Doing further checks\n"%(id, i))
                            co_orthologs_result, distance_ref_hit, distance_hit_query = checkCoOrthologs(gene_name, id, i, fdog_ref_species, candidatesOutFile, msaTool)
                            if co_orthologs_result == 1:
                                print("\t Distance query - blast hit: %6.4f, Distance blast hit - reference: %6.4f\tAccepting\n"%(distance_hit_query, distance_ref_hit))
                                orthologs.append(gene_name)
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

            os.system("blastp -db " + blast_dir_path + species + "/" + species + " -outfmt '6 sseqid qseqid evalue' -max_target_seqs 10 -out tmp/blast_" + species + " -evalue " + str(evalue_cut_off) + " -query " + candidatesOutFile)
            alg_file = open("tmp/blast_" + species, "r")
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
    return orthologs, seed

def addSequences(sequenceIds, candidate_fasta, core_fasta, output, name, species_list):
    #print("addSequences")
    #print(sequenceIds)
    #print(species_list)
    seq_records_core = readFasta(core_fasta)
    output_file = open(output + "/" + name + ".extended.fa", "a+")

    seq_records_core = list(seq_records_core)
    for species in species_list:
        for entry_core in seq_records_core:
            if species in entry_core.id:
                output_file.write(">" + entry_core.id + "\n")
                output_file.write(str(entry_core.seq) + "\n")

    seq_records_candidate = readFasta(candidate_fasta)
    seq_records_candidate = list(seq_records_candidate)
    for entry_candidate in seq_records_candidate:
        #print(entry_candidate.id)
        #print(sequenceIds)
        if entry_candidate.id.split("|")[2] in sequenceIds:
            output_file.write(">" + entry_candidate.id + "\n")
            output_file.write(str(entry_candidate.seq) + "\n")
    output_file.close()
    return 0

def createFasInput(orthologsOutFile, mappingFile):
    with open(orthologsOutFile, "r") as f:
        fas_seed_id = (f.readline())[1:-1]

    mappingFile = open(mappingFile, "w")

    seq_records = readFasta(orthologsOutFile)
    for seq in seq_records:
        ncbi_id = (seq.id.split("@"))[1]
        mappingFile.write(seq.id + "\t" + "ncbi" + ncbi_id + "\n")


    return fas_seed_id

def cleanup(tmp):
    if tmp == False:
        os.system('rm -r tmp/')

def checkOptions():
    pass
    #muss ich unbedingt noch ergänzen wenn ich alle möglichen input Optionen implementiert habe!!!



def main():

    #################### some variables ########################################
    ### for testing only
    #core-ortholog group name
    group = "778452"
    #species name assemblie (folder name in assemby folder)
    species_name = "L.pustulata"
    #assembly species_name
    assembly_name = "contigs.fa"
    assembly_path = "../data/assembly_dir/"+ species_name + "/" + assembly_name
    augustus_ref_species = "saccharomyces_cerevisiae_S288C"
    #cut_off_merging_candidates = 500
    average_intron_length = 5000
    length_extension = 5000
    tmp = False
    strict = False
    evalue = 0.00001
    out = os.getcwd()
    taxa = []
    aligner = "blast"
    checkCoorthologs = False
    msaTool = "muscle"
    matrix = 'blosum62'
    hmmpath = "./data/core_orthologs/"

    ########################### handle user input ##############################
    #user input core_ortholog group
    # I nedd a function that checks if user input is allowed
    #print(sys.argv)
    input = sys.argv

    for i in range(1,len(input)):
        if input[i] == "--assembly":
            assembly_path = input[i+1]
        elif input[i] == "--gene":
            group = input[i+1]
        elif input[i] == "--augRefSpecies":
            augustus_ref_species = input[i+1]
        elif input[i] == "--avIntron":
            average_intron_length = int(input[i+1])
        elif input[i] == "--lengthExtension":
            length_extension = int(input[i+1])
        elif input[i] == "--tmp":
            tmp = True
        elif input[i] == "--name":
            fdog_name =  input[i+1]
        elif input[i] == "--out":
            out = input[i+1]
        elif input[i] == "--refSpecies":
            fdog_ref_species = input[i+1]
        elif input[i] == "--strict":
            strict = True
        elif input[i] == "--evalue":
            evalue = input[i+1]
        elif input[i] == "--searchTaxa":
            taxa = input[i+1].split(",")
            #print(taxa)
        elif input[i] == "--aligner":
            aligner = input[i+1]
        elif input[i] == "--checkCoOrthologsRef":
            checkCoorthologs = True
        elif input[i] == "--msaTool":
            msaTool = input[i+1]
        elif input[i] == "--scoringmatrix":
            matrix = input[i+1]
        elif input[i] == "--hmmpath":
            hmmpath = input[i+1]
        elif input[i] == "--help":
            print("Parameters: \n")
            print("--assembly: path to assembly input file in fasta format \n")
            print("--gene: core_ortholog group name. Has to be located in data/core_orthologs\n")
            print("--augRefSpecies: reference species for augustus\n")
            print("--avIntron: average intron length of the selected species in bp (default: 5000)\n")
            print("--lengthExtension: length extension of the candidate regions in bp (default:5000)\n")
            print("--tmp: tmp files will not be deleted")
            print("--name: Species name according to the fdog naming schema [Species acronym]@[NCBI ID]@[Proteome version]")
            print("--refSpecies: fDOG reference species")
            print("--out: path to the output folder")
            print("--aligner: can use blast or diamond as an aligner in the backward search (default: blast)")
            print("--msaTool: {mafft-linsi, muscle}. Choose between mafft-linsi and muscle for the multiple seqeunce alignemnts. (default:muscle)")
            print("--evalue: evalue cut off for every blast search, (default: 0.00001)")
            print("--scoringmatrix: ..., (default:blosum62)")
            return 0


    ########################## paths ###########################################

    #open core_ortholog group
    msa_path = hmmpath + group +"/"+ group + ".aln"
    hmm_path = hmmpath + group +"/hmm_dir/"+ group + ".hmm"
    fasta_path = hmmpath + group +"/"+ group + ".fa"
    consensus_path = "tmp/" + group + ".con"
    profile_path = "tmp/" + group + ".prfl"
    path_assembly = assembly_path
    candidatesOutFile = "tmp/" + group + ".candidates.fa"
    orthologsOutFile = out + "/" + group + ".extended.fa"
    fasOutFile = out + "/" + group
    mappingFile = out + "/tmp/" + group + ".mapping.txt"

    os.system('mkdir tmp')

    ######################## check user input ##################################

    #the user input has to be checked!!

    ######################## consensus sequence ################################

    #make a majority-rule consensus seqeunce with the tool hmmemit from hmmer
    print("Building a consensus sequence \n")
    os.system('hmmemit -c -o' + consensus_path + ' ' + hmm_path)
    print("consensus sequence is finished\n")

    ######################## block profile #####################################
    print("Building a block profile \n")

    os.system('msa2prfl.pl ' + msa_path + ' --setname=' + group + ' >' + profile_path)
    #print(os.path.getsize(profile_path))
    if int(os.path.getsize(profile_path)) > 0:
        print("block profile is finished \n")
    else:
        print("Building block profiles failed. Using prepareAlign to convert alignment\n")
        new_path = "../data/core_orthologs/" + group +"/"+ group + "_new.aln"
        os.system('prepareAlign < ' + msa_path + ' > ' + new_path)
        os.system('msa2prfl.pl ' + new_path + ' --setname=' + group + ' >' + profile_path)
        print("block profile is finished \n")

    ######################## tBLASTn ###########################################

    #database anlegen

    db_check = searching_for_db(path_assembly)
    if db_check == 0:
        print("creating a blast data base \n")
        os.system('makeblastdb -in ' + path_assembly + ' -dbtype nucl -parse_seqids -out ' + path_assembly)
        print("database is finished \n")
    else:
        print('blast data base exists already, continuing...')


    #make a tBLASTn search against the new database
    #codon table argument [-db_gencode int_value], table available ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt

    print("tBLASTn search against new created data base")
    os.system('tblastn -db ' + path_assembly + ' -query ' + consensus_path + ' -outfmt "6 sseqid sstart send evalue qstart qend " -out tmp/blast_results.out')
    print("tBLASTn search is finished")

    ################### search for candidate regions and extract seq ###########

    # parse blast and filter for candiate regions
    regions, number_regions = candidate_regions(average_intron_length, evalue)

    if regions == 0:
        #no candidat region are available, no ortholog can be found
        print("No candidate region found")
        cleanup(tmp)
        return 0

    else:
        print(str(number_regions) + " candiate regions were found. Extracting sequences.")
        extract_seq(regions, path_assembly)

    ############### make Augustus PPX search ###################################
    print("starting augustus ppx \n")
    augustus_ppx(regions, candidatesOutFile, length_extension, profile_path, augustus_ref_species, fdog_name, group)
    print("augustus is finished \n")

    ################# bachward search to filter for orthologs###################
    #verschiede Modi beachten!
    reciprocal_sequences, taxa = backward_search(candidatesOutFile, fasta_path, strict, fdog_ref_species, evalue, taxa, aligner, checkCoorthologs, msaTool, matrix)
    if reciprocal_sequences == 0:
        cleanup(tmp)
        return 0

    ################ add sequences to extended.fa in the output folder##########
    addSequences(reciprocal_sequences, candidatesOutFile, fasta_path, out, group, taxa)

    ############### make Annotation with FAS ###################################
    fas_seed_id = createFasInput(orthologsOutFile, mappingFile)

    os.system('mkdir tmp/anno_dir')
    os.system('calcFAS --seed ' + fasta_path + ' --query ' + orthologsOutFile + ' --annotation_dir tmp/anno_dir --bidirectional --phyloprofile ' + mappingFile + ' --seed_id "' + fas_seed_id + '" --out_dir ' + out + ' --out_name ' + group )


    ################# remove tmp folder ########################################

    cleanup(tmp)


if __name__ == '__main__':
    main()
