############################ imports ###########################################
import os

########################### handle user input ##################################

#user input core_ortholog group
#have to add an input option

#core-ortholog group name
group = "778452"

#species name assemblie (folder name in assemby folder)
species_name = "L.pustulata"

#assembly species_name
assembly_name = "contigs.fa"


########################## paths ###############################################

#open core_ortholog group
msa_path = "../data/core_orthologs/" + group +"/"+ group + ".aln"
hmm_path = "../data/core_orthologs/" + group +"/hmm_dir/"+ group + ".hmm"
consensus_path = "tmp/" + group + ".con"
profile_path = "tmp/" + group + ".prfl"
path_assembly = "../data/assembly_dir/" + species_name + "/" + assembly_name
augustus_path = "msa2prfl.pl"

os.system('mkdir tmp')


#msa = open("../data/core_orthologs/" + group +"/"+ group + ".aln", "r")
#lines = msa.readlines()
#msa.close()

######################## consensus sequence ####################################

#make a majority-rule consensus seqeunce with the tool hmmemit from hmmer
print("Building a consensus sequence \n")
os.system('hmmemit -c -o' + consensus_path + ' ' + hmm_path)
print("consensus seqeunce is finished\n")

######################## block profile #########################################
print("Building a block profile \n")
os.system('msa2prfl.pl ' + msa_path + ' --setname=' + group + ' >' + profile_path)
print("block profile is finished \n")
######################## tBLASTn ###############################################

#database anlegen
print("creating a blast database \n")
os.system('makeblastdb -in ' + path_assembly + ' -dbtype nucl -out ' + path_assembly)
print("database is finished \n")

#make a tBLASTn search against the new database

os.system('tblastn -db ' + path_assembly + ' -query ' + consensus_path + ' -outfmt "6 sseqid sstart send evalue bitscore" -out tmp/blast_results.out')

###################### extracting candidate regions ############################
# info about output blast http://www.metagenomics.wiki/tools/blast/blastn-output-format-6

#parsing blast output
# format <contig> <start> <end> <evalue> <score>
