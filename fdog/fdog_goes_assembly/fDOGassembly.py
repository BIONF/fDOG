############################ imports ###########################################
import os
import os.path
import sys
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



    print(blast_results)
    return blast_results, number_regions


def merge_regions(blast_results, cut_off):
    number_regions = 0
    for key in blast_results:
        locations = blast_results[key]
        size_list = len(locations)
        i = 0
        j = 1
        old_size = 0
        while size_list != old_size and i < size_list:
            old_size = size_list
            #print(locations)
            while j < size_list:

                # breakup point? or we have to skip this j
                if (i == j) and (j + 1 < size_list):
                    j+=1
                elif (i == j):
                    break

                if (locations[i][0] < locations[j][0]) and (locations[i][1] > locations[j][0]):
                    # start is between start and end -> merge
                    locations[i][1] = max(locations[j][1], locations[i][1])
                    locations[i][2] = min(locations[j][2], locations[i][2])
                    locations.pop(j)
                    j -= 1
                elif (locations[i][0] < locations[j][1]) and (locations[i][1] > locations[j][1]):
                    #end is between start and end -> merge
                    locations[i][0] = min(locations[j][0], locations[i][0])
                    locations[i][2] = min(locations[j][2], locations[i][2])
                    locations.pop(j)
                    j -= 1
                elif (locations[i][0] > locations[j][1]) and (locations[i][0] - locations[j][1] <= cut_off):
                    # end is not more than cut-off distanced
                    locations[i][0] = locations[j][0]
                    locations[i][2] = min(locations[j][2], locations[i][2])
                    locations.pop(j)
                    j -= 1
                elif (locations[i][1] < locations[j][0] and locations[j][0] - locations[i][1] <= cut_off):
                    # start is not more than cut-off distanced
                    locations[i][0] = locations[j][0]
                    locations[i][2] = min(locations[j][2], locations[i][2])
                    locations.pop(j)
                    j -= 1
                j += 1
                size_list = len(locations)

            i += 1
            j = 0
        number_regions += size_list

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


def candidate_regions(cut_off, intron_length):
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
        if not evalue <= 0.00001:
            break
    if blast_results == {}:
        return 1,0
    else:
        candidate_regions, number_regions = merge(blast_results, intron_length)
        #candidate_regions, number_regions = merge_regions(blast_results, cut_off)
        #print(candidate_regions, number_regions)
        return candidate_regions, number_regions


def extract_seq(region_dic, path):
    #print(region_dic)
    for key in region_dic:
        print("blastdbcmd -db " + path + " -dbtype 'nucl' -entry " + key + " -out tmp/" + key + ".fasta -outfmt %f")
        os.system("blastdbcmd -db " + path + " -dbtype 'nucl' -entry " + key + " -out tmp/" + key + ".fasta -outfmt %f")


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
    cut_off_merging_candidates = 500
    average_intron_length = 5000
    length_extension = 5000
    tmp = False

    ########################### handle user input ##############################
    #user input core_ortholog group
    #have to add an input option
    #print(sys.argv)
    input = sys.argv

    for i in range(1,len(input)):
        if input[i] == "--assembly":
            assembly_path = input[i+1]
        elif input[i] == "--gene":
            group = input[i+1]
        elif input[i] == "--refSpecies":
            augustus_ref_species = input[i+1]
        elif input[i] == "--avIntron":
            average_intron_length = int(input[i+1])
        elif input[i] == "--lengthExtension":
            length_extension = int(input[i+1])
        elif input[i] == "--tmp":
            tmp = True
        elif input[i] == "--help":
            print("Parameters: \n")
            print("--assembly: path to assembly input file in fasta format \n")
            print("--gene: core_ortholog group name. Has to be located in data/core_orthologs\n")
            print("--refSpecies: reference species for augustus\n")
            print("--avIntron: average intron length of the selected species in bp (default: 5000)\n")
            print("--lengthExtension: length extension of the candidate regions in bp (default:5000)\n")
            print("--tmp: tmp files will not be deleted")
            return 0


    ########################## paths ###########################################

    #open core_ortholog group
    msa_path = "../data/core_orthologs/" + group +"/"+ group + ".aln"
    hmm_path = "../data/core_orthologs/" + group +"/hmm_dir/"+ group + ".hmm"
    consensus_path = "tmp/" + group + ".con"
    profile_path = "tmp/" + group + ".prfl"
    path_assembly = assembly_path
    outfile = group + ".candidates.fa"

    os.system('mkdir tmp')

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
        os.system('msa2prfl.pl ' + msa_path + ' --setname=' + group + ' >' + profile_path)
        print("block profile is finished \n")

    ######################## tBLASTn ###########################################

    #database anlegen
    print("creating a blast database \n")
    os.system('makeblastdb -in ' + path_assembly + ' -dbtype nucl -parse_seqids -out ' + path_assembly)
    print("database is finished \n")

    #make a tBLASTn search against the new database

    print("tBLASTn search against new created database")
    os.system('tblastn -db ' + path_assembly + ' -query ' + consensus_path + ' -outfmt "6 sseqid sstart send evalue qstart qend " -out tmp/blast_results.out')
    print("tBLASTn search is finished")

    ################### search for candidate regions and extract seq ###########

    # parse blast and filter for candiate regions
    regions, number_regions = candidate_regions(cut_off_merging_candidates, average_intron_length)

    if regions == 1:
        #no candidat region are available, no ortholog can be found
        print("No candidate region found")
        os.system('rm -r tmp/')
        return 1

    else:
        print(str(number_regions) + " candiate regions were found. Extracting sequences.")
        extract_seq(regions, path_assembly)

    ############### make Augustus PPX search ###################################

    output = open(outfile, "w")

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
                    header = name + " " + id
                    output.write(header)
                else:
                    output.write(line)
            sequence_file.close()

    output.close()




    ################# remove tmp folder ########################################

    if tmp == False:
        os.system('rm -r tmp/')


if __name__ == '__main__':
    main()
