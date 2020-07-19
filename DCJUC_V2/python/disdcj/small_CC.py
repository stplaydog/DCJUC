####
#### copy right by zhaoming yin
#### this code can handle both regular and irregular components
####

from CC_type import determine_cc_type
import sys
import os
import itertools
import random

def check_key(key, item):
	result = True
	if len(item[0][key])==1 or (key in item[1] and len(item[1][key])==1):
		result = True
	else:
		result = False
	if len(item[1][key])==1 or (key in item[0] and len(item[0][key])==1):
		result = True
	else:
		result = False
	return result

def select_exem_edge(item, key, color):
	result = []
	for i in item[color][key]:
		# remove others except i
		tmp = item
		tmp[color][key] = i
		for j in tmp[color][key]: #change connection except for i
			if j != i:
				tmp_j = []
				for k in tmp[color][j]:
					if k != key:
						tmp_j.append(k)
				tmp[color][j] = tmp_j
		result.append(tmp)

def disambi_exem(key, item):
	result = []
	if check_key(key, item)==True:
		return False, result
	else:
		##edge one
		disambi_one = select_exem_edge(item, key, 0)
		for interim in disambi_one:
			disambi_two = select_exem_edge(interim, key, 1)
			for finished in disambi_two:
				result.append(finished)
		return result

def get_mapping(key, item):
	list_one = []
	list_two = []
	#get list_one first
	list_one = list(itertools.permutations(item[0][key]))
	#then get list_two
	list_two = list(itertools.permutations(item[1][key]))
	return [list_one, list_two]

def perform_matching(key, item, mappings):
	result = []
	for i in range(len(mappings[0])):
		tmp = item
		for j in range(1, len(mappings[0][i])):
			new_vet = str(random.randint(100000,200000))
			while new_vet in item[0] or new_vet in item[1]:
				new_vet = str(random.randint(100000,200000))
			#color 0
			for k in range(len(tmp[0][mappings[0][i][j]])):
				if tmp[0][mappings[0][i][j]][k] == key:
					tmp[0][mappings[0][i][j]][k] = new_vet
			#color 1
			for k in range(len(tmp[1][mappings[1][i][j]])):
				if tmp[1][mappings[1][i][j]][k] == key:
					tmp[1][mappings[1][i][j]][k] = new_vet
		#keep only the first edge
		tmp[0][key] = [mappings[0][i][0]]
		tmp[1][key] = [mappings[1][i][0]]
		result.append(tmp)
	return result	
	
def disambi_matc(key, item):
	if check_key(key, item)==True:
		return False, []
	else:
		##get all combinations
		##we need some randomization methods now :)
		#get a mapping first
		mappings = get_mapping(key, item)
		#perform the matching
		result = perform_matching(key, item, mappings)
	return result

def check_final(item):
	result = True
	for key in item[0]:
		if len(item[0][key])==1 or (key in item[1] and len(item[1][key])==1):
			result = True
		else:
			result = False
	for key in item[1]:
		if len(item[1][key])==1 or (key in item[0] and len(item[0][key])==1):
			result = True
		else:
			result = False
	return result

#disambiguate the graph and compute the number of component in the graph
#every combination is a tmp graph
def compute_small_cc(adj, algo):
	has_multi_comp = True
	lis = []
	lis_final = []
	lis.append(adj)
	while has_multi_comp:
		item = lis.pop()
		is_ambi = False
		new_item = []
		for key in item[0]:
			if algo=='exem':
				is_ambi, new_item = disambi_exem(key, item)
			elif algo == "matc":
				is_ambi, new_item = disambi_matc(key, item)
			if is_ambi == True:
				break
		for tmp in new_item:
			if check_final(tmp) == True:
				lis_final.append(tmp)
			else:
				lis.append(tmp)
	max_num_cc = 0
	for item in lis_final:
		num_cc = 0.0
		vet_visited = {}
		for key in item[0]:
			if key not in vet_visited:
				visited, num = determine_cc_type(item, key)	
				num_cc += num
				for v in visited:
					vet_visited[v] = True
		if num_cc > max_num_cc:
			max_num_cc = num_cc
	return max_num_cc
