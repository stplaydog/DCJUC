import sys
import math
import copy
import random
import copy

#########################################
# perform mutation, indel and duplication 
# on a genome. 
#########################################
def gen_genome (genome, gene_len, len, max_len, theta, gamma, phi, dup_dic, dup_id):
	mut_rate = int(gene_len*(math.sqrt(float(len)/float(max_len))*theta))
	del_rate = int(gene_len*(math.sqrt(float(len)/float(max_len))*gamma))
	dup_rate = int(gene_len*(math.sqrt(float(len)/float(max_len))*phi))
	while (mut_rate+del_rate+dup_rate)>0: 
		which = random.randint(0, 3)
		if which == 0 and del_rate>0:
			delete (genome, 1, dup_dic, dup_id)
			del_rate -= 1
		if which == 1 and mut_rate>0:
			mutate (genome, 1, dup_id)
			mut_rate -= 1
		if which == 2 and dup_rate>0:
			duplicate (genome, 1, dup_id)
			dup_rate -= 1
	
############################################
# perform mutation, which are generally the 
# inversion operations.
############################################
def mutate (genome, num, dup_id):
	for i in range(num):
		start = random.randint(0, len(genome)-1)
		end = random.randint(start, len(genome)-1)
		sub_g = genome[start:end]
		sub_g.reverse()
        b = copy.deepcopy(dup_id[start:end])
        b = map(b.reverse(), b)
        dup_id[start:end] = b
        for i in range(end-start):
            genome[start+i] = sub_g[i]*-1

############################################
# perform indel operations
############################################
def delete (genome, num, dup_dic, dup_id):
	for i in range(num):
		pos = random.randint(0, len(genome)-1) 
		while abs(genome[pos]) in dup_dic: 
            # when 'balanced', can not remove duplicated genes
            # when 'unbalanced', there is no such restriction.
			pos = random.randint(0, len(genome)-1)
		genome.pop(pos)
        dup_id.pop(pos)

############################################
# perform duplication operations
############################################
def duplicate (genome, num, dup_id):
	for i in range(num):
		pos = random.randint(0, len(genome)-1)
		gene = genome[pos]
		genome.insert(pos+1, gene)	

#########################################
# pre decide which genes to be duplicated
#########################################
def get_dup_list(num_gene, num_dup, str_mode):
    dup_dic = []
    dup_list = []
    for i in range(num_dup):
        dup_dic.append({})
        dup_list.append([])
    for i in range(num_dup):
        gene = random.randint(1, num_gene)
        for j in range(num_dup):
            dup_list[j].append(gene)
            if str_mode == 'balanced':
                dup_dic[j][gene] = True
            elif str_mode == 'unbalanced':
                # There is no need to indicate that this gene is duplicated,
                # since it might be deleted in the future.
                dup_dic[j][gene] = False
    return dup_list,dup_dic

##########################################
# add duplication back
##########################################
def pre_dup_process(genome, dup_list):
	for item in dup_list:
		size = len(genome)
		for i in range(size):
			if genome[i] == item:
				genome.insert(i+1, item)
				break

#################################################################
# give identity to each genes in a gene family
# dup_identity[which_genome][which_gene][gene_name:gene_identity]
#################################################################
def assign_gf_identity(genome, num_gf, dup_id, dup_dic):
    gene_count = {}
    for i in range(num_gf):
        gene_count[i+1] = 1
    for i in range(len(genome)):
        gene = genome[i]
        gc = gene_count[gene]
        dup_id.append([gene, gc])
        gene_count[gene] += 1
