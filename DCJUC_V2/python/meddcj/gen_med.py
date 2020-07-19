import sys
sys.path.append('/scratch/zyin/optkit/optkit/python/disdcj/')
from order_gen import gen_genome
from order_gen import mutate
from order_gen import delete
from order_gen import duplicate
from order_gen import get_dup_list
from order_gen import pre_dup_process
from order_gen import assign_gf_identity
from o_2_g import o_to_g

##
# This function defines the process to generate three gene orders based on a identity 'parent' order
##
def gen_three_orders_n_graph(gene_len, thetabase, gamma, phi, seq_folder, str_mode, graph_folder):
    for i in range(10):
        # 10 different mutation rates
        theta = (i+1)*thetabase
        # each mutation rate has 10 cases
        for j in range(10):
            # prepare output files
            file = seq_folder + str(theta) +"_"+str(gamma)+"_"+str(phi)+"_"+str(j) 
            print file
            writer = open(file, "w")
            # prepare some data structures
            genomes = [[]]*(3)
            dup_list,dup_dic = get_dup_list(gene_len, int(gene_len*phi), str_mode, 3)
            dup_id = [[],[],[]]
            # perform evolutionary operations
            for k in range(3):
                genomes[k] = list(xrange(gene_len))
                for l in range(gene_len):
                    genomes[k][l]+=1
                #perform duplication first
                pre_dup_process(genomes[k], dup_list[k])
                #perform assigning indentity to each gene in each gene family
                assign_gf_identity(genomes[k], gene_len, dup_id[k], dup_dic[k])
                #real evolution part
                gen_genome (genomes[k], gene_len, 1, 1, theta, gamma, 0, dup_dic[k], dup_id[k])
                #write things back
                writer.write(">"+str(k)+"\n")
                for m in range(len(genomes[k])):
                    writer.write(str(genomes[k][m])+" ")
                writer.write("\n")
            writer.close()
            # let's turn order into graph   
            order_reader = open(file)
            graph_file = graph_folder + str(theta) +"_"+str(gamma)+"_"+str(phi)+"_"+str(j) 
            graph_writer = open(graph_file, 'w')
            o_to_g(order_reader, graph_writer)

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
#this is to indicate whether two genomes are 'balanced' or 'imbalanced'
str_mode = sys.argv[6]
# 
graph_folder = sys.argv[7]
# call the function
gen_three_orders_n_graph(gene_len, thetabase, gamma, phi, seq_folder, str_mode, graph_folder)
