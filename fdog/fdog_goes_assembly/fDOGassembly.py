############################ imports ###########################################
import os
import sys
########################### functions ##########################################


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
            start = locations[i][0]
            end = locations[i][1]

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
    # format blast line:  <contig> <start> <end> <evalue> <score>
    #fomrat dictionary: {node_name: [(<start>,<end>)]}
    #print(line)
    line = line.replace("\n", "")
    line_info = line.split("\t")
    #print(line_info)
    evalue = float(line_info[3])

    #cut off
    if evalue > 0.0001:
        return blast_results, evalue
    #add region to dictionary
    else:
        node_name, start, end = line_info[0], line_info[1], line_info[2]
        if node_name in blast_results:
            list = blast_results[node_name]
            list.append([int(start),int(end), evalue])
            blast_results[node_name] = list
        else:
            blast_results[node_name] = [[int(start),int(end), evalue]]

    return blast_results, evalue


def candidate_regions(cut_off):
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
        return 1
    else:
        candidate_regions, number_regions = merge_regions(blast_results, cut_off)
        #print(candidate_regions, number_regions)
        return candidate_regions, number_regions


def extract_seq(region_dic, path):
    #print(region_dic)
    for key in region_dic:
        os.system("blastdbcmd -db " + path + " -dbtype 'nucl' -entry " + key + " -out tmp/" + key + ".fasta -outfmt %f")


def main():

    ########################### handle user input ##############################
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


    #user input core_ortholog group
    #have to add an input option
    print(sys.argv)
    input = sys.argv

    for i in range(1,len(input)):
        if input[i] == "--assembly":
            assembly_path = input[i+1]
        elif input[i] == "--gene":
            group == input[i+1]
        elif input[i] == "--refSpecies":
            augustus_ref_species = input[i+1]
        elif input[i] == "--cut_off":
            cut_off_merging_candidates = input[i+1]





    ########################## paths ###########################################

    #open core_ortholog group
    msa_path = "../data/core_orthologs/" + group +"/"+ group + ".aln"
    hmm_path = "../data/core_orthologs/" + group +"/hmm_dir/"+ group + ".hmm"
    consensus_path = "tmp/" + group + ".con"
    profile_path = "tmp/" + group + ".prfl"
    path_assembly = assembly_path

    os.system('mkdir tmp')


    ######################## consensus sequence ################################

    #make a majority-rule consensus seqeunce with the tool hmmemit from hmmer
    print("Building a consensus sequence \n")
    os.system('hmmemit -c -o' + consensus_path + ' ' + hmm_path)
    print("consensus seqeunce is finished\n")

    ######################## block profile #####################################
    print("Building a block profile \n")

    os.system('msa2prfl.pl ' + msa_path + ' --setname=' + group + ' >' + profile_path)
    print("block profile is finished \n")
    ######################## tBLASTn ###########################################

    #database anlegen
    print("creating a blast database \n")
    os.system('makeblastdb -in ' + path_assembly + ' -dbtype nucl -parse_seqids -out ' + path_assembly)
    print("database is finished \n")

    #make a tBLASTn search against the new database

    os.system('tblastn -db ' + path_assembly + ' -query ' + consensus_path + ' -outfmt "6 sseqid sstart send evalue bitscore" -out tmp/blast_results.out')

    ################### search for candidate regions and extract seq ###########

    # parse blast and filter for candiate regions
    regions, number_regions = candidate_regions(cut_off_merging_candidates)

    if regions == 1:
        #no candidat region are available, no ortholog can be found
        print("No candidate region found")
        os.system('rm -r tmp/')
        return 1

    else:
        print(str(number_regions) + " candiate regions were found. Extracting sequences.")
        extract_seq(regions, path_assembly)

    ############### make Augustus PPX search ####################################
    for key in regions:
        locations = regions[key]
        counter = 0
        for i in locations:
            counter += 1
            start = str(i[0])
            end = str(i[1])
            if start < end:
            #print("augustus --proteinprofile=" + profile_path + " --predictionStart=" + start + " --predictionEnd=" + end + " --species=" + augustus_ref_species + " tmp/" + key + ".fasta > tmp/" + key + ".gff")
                os.system("augustus --proteinprofile=" + profile_path + " --predictionStart=" + start + " --predictionEnd=" + end + " --species=" + augustus_ref_species + " tmp/" + key + ".fasta > tmp/" + key + "_" + str(counter) + ".gff")
            else:
                os.system("augustus --proteinprofile=" + profile_path + " --predictionStart=" + end + " --predictionEnd=" + start + " --species=" + augustus_ref_species + " tmp/" + key + ".fasta > tmp/" + key + "_" + str(counter) + ".gff")

    ################# remove tmp folder ########################################

    #have to be added after program ist finished, maybe use parametere so that the user can turn it off

if __name__ == '__main__':
    main()
