####
### Copy right by Zhaoming Yin
### This code is to do some preprocess for optimization for distance computation
####

import sys
import os
from graph import read_graph
from graph import vis_graph
from graph import delete_reg_CC
#construct the graph
#it's a bi-direction graph



####output number of reg CC to a speicific file
def write_interim_result(num_reg_CC, interim_file):
	writer = open(interim_file, "w")
	writer.write("num_reg_CC: "+str(num_reg_CC)+"\n")

####main entrance####
#####################
input_file = sys.argv[1]
out_dir = sys.argv[2]
interim_file = sys.argv[3]
vis_file = sys.argv[4] + ".dot"

adj = read_graph(input_file)
#vis_graph(adj, vis_file)
num_reg_CC,count = delete_reg_CC(adj, out_dir, "")
write_interim_result(num_reg_CC, interim_file)

##additional debug info
#adj = read_graph("/Users/zhaomingyin/Dropbox/research/code/optec-code/data/dist/irrCC/200_0.1_0.1_0.1/1.0_0.1_0.1_0/0.gr")
#vis_graph(adj, "/Users/zhaomingyin/Dropbox/research/code/optec-code/vis0.dot")
