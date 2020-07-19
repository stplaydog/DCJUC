import sys
import math
import copy
import random
import copy

from order_gen import gen_genome
from order_gen import mutate
from order_gen import delete
from order_gen import duplicate
from order_gen import get_dup_list
from order_gen import pre_dup_process
from order_gen import assign_gf_identity


###########################################
# this is the main process
###########################################
#this is the number of gene family
gene_len = int(sys.argv[1])
#this is the mutation rate
thetabase = float(sys.argv[2])
#this is the indel rate
gamma = float(sys.argv[3])
#this is the duplication rate
phi = float(sys.argv[4])
#this is where to store gene sequences
seq_folder = sys.argv[5]
#this is evolutionary mode 'dual' or 'sing'
evo_mode = sys.argv[6]
#this is to indicate whether two genomes are 'balanced' or 'imbalanced'
str_mode = sys.argv[7]
#this is to set up the bijection result folder
bijection_folder = sys.argv[8]

for i in range(10):
    theta = (i+1)*thetabase
    for j in range(10):
        # prepare output files
        file1 = seq_folder + str(theta) +"_"+str(gamma)+"_"+str(phi)+"_"+str(j) 
        writer1 = open(file1, "w")
        file2 = bijection_folder + str(theta) +"_"+str(gamma)+"_"+str(phi)+"_"+str(j) 
        print file2
        writer2 = open(file2, "w")
        # prepare some data structures
        genomes = [[]]*(2)
        dup_list,dup_dic = get_dup_list(gene_len, int(gene_len*phi), str_mode)
        dup_id = [[],[]]
        # perform evolutionary operations
        for k in range(2):
            genomes[k] = list(xrange(gene_len))
            for l in range(gene_len):
                genomes[k][l]+=1
            #perform duplication first
            pre_dup_process(genomes[k], dup_list[k])
            #perform assigning indentity to each gene in each gene family
            assign_gf_identity(genomes[k], gene_len, dup_id[k], dup_dic[k])
            #real evolution part
            if evo_mode == 'sing' and k == 1:
                gen_genome (genomes[k], 1, 1, 1, theta, gamma, 0, dup_dic[k], dup_id[k])
            elif evo_mode == 'dual':
                gen_genome (genomes[k], 1, 1, 1, theta/2, gamma/2, 0, dup_dic[k], dup_id[k])
            #write things back
            writer1.write(">"+str(k)+"\n")
            for m in range(len(genomes[k])):
                writer1.write(str(genomes[k][m])+" ")
            writer1.write("\n")
            writer2.write(">"+str(k)+"\n")
            for m in range(len(dup_id[k])):
                writer2.write(str(dup_id[k][m][0])+":"+str(dup_id[k][m][1])+" ")
            writer2.write("\n")
