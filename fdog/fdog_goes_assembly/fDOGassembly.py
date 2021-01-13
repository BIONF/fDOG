############################ imports ###########################################
import os

########################### handle user input ##################################

#user input core_ortholog group
#have to add an input option

group = "778452"

########################## paths ###############################################

#open core_ortholog group
msa_path = "../data/core_orthologs/" + group +"/"+ group + ".aln"
hmm_path = "../data/core_orthologs/" + group +"/hmm_dir/"+ group + ".hmm"
consensus_path = "../data/core_orthologs/" + group + "/" + group + ".con"
profile_path = "../data/core_orthologs/" + group + "/" + group + ".prfl"
augustus_path = "msa2prfl.pl"


#msa = open("../data/core_orthologs/" + group +"/"+ group + ".aln", "r")
#lines = msa.readlines()
#msa.close()

######################## consensus sequence ####################################

#make a majority-rule consensus seqeunce with the tool hmmemit from hmmer

os.system('hmmemit -c -o' + consensus_path + ' ' + hmm_path)

######################## block profile #########################################

os.system('msa2prfl.pl ' + msa_path + ' --setname=' + group + ' >' + profile_path)

######################## tBLASTn ###############################################
