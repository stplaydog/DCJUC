####
### Copy right by Zhaoming Yin
### This is the code to enumerate number of connected component in BPG (breakpoint graph)
### to prove that we can process the connected component separately
### we need to prove that in practical, there are actually enough number of CC to be processed.
####

import sys

reader = open(sys.argv[1])

lines = reader.readlines()

adj1 = {}
adj2 = {}
visited1 = {}
visited2 = {}

#construct the graph
#it's a bi-direction graph

def remove_duplicate(adj, vis):
	tmp = []
	tmp_vis = []
	for i in range(len(adj)):
		dup = False
		for j in range(i+1, len(adj)):
			if adj[i] == adj[j]:
				dup = True
		if dup == False:
			tmp.append(adj[i])
			tmp_vis.append(False)
	return tmp,tmp_vis

for i in range(1, len(lines)):
	items = lines[i].strip().split(" ")
	if items[2] == '1':
		if items[0] not in adj1:
			adj1[items[0]] = [items[1]]
			visited1[items[0]] = [False]
		else:
			adj1[items[0]].append(items[1])
			visited1[items[0]].append(False)
		if items[1] not in adj1:
			adj1[items[1]] = [items[0]]
			visited1[items[1]] = [False]
		else:
			adj1[items[1]].append(items[0])
			visited1[items[1]].append(False)
	if items[2] == '2':
		if items[0] not in adj2:
			adj2[items[0]] = [items[1]]
			visited2[items[0]] = [False]
		else:
			adj2[items[0]].append(items[1])
			visited2[items[0]].append(False)
		if items[1] not in adj2:
			adj2[items[1]] = [items[0]]
			visited2[items[1]] = [False]
		else:
			adj2[items[1]].append(items[0])
			visited2[items[1]].append(False)

for key in adj1:
	adj1[key], visited1[key] = remove_duplicate(adj1[key], visited1[key])
for key in adj2:
	adj2[key], visited2[key] = remove_duplicate(adj2[key], visited2[key])

#enumerate number of connected components
#bfs should be used
num_component = 0
vet_map = {}
count = 0
while True:
	# find the seed first
	seed = -1
	for key in adj1:
		if key not in vet_map:
			seed = key
			break
	if seed == -1:
		for key in adj2:
			if key not in vet_map:
				seed = key
				break
	# then using the queue to perform bfs
	if seed != -1:
		is_multi_component = False
		queue = [seed]
		vet_map[seed] = True
		count += 1
		print '->', seed,":",count
		while len(queue) != 0:
			idx = len(queue)-1	
			elem = queue[idx]
			queue.pop(idx)
			#print 'pop',elem
			if (elem in adj1 and len(adj1[elem]) > 1) or (elem in adj2 and len(adj2[elem]) > 1):
				is_multi_component = True
			if elem in adj1:
				for i in range(len(adj1[elem])):
					if adj1[elem][i] not in vet_map:
						queue.append(adj1[elem][i])
						vet_map[adj1[elem][i]] = True
						visited1[elem][i] = True
						print adj1[elem][i] 
			if elem in adj2:
				for i in range(len(adj2[elem])):
					if adj2[elem][i] not in vet_map:
						queue.append(adj2[elem][i])
						vet_map[adj2[elem][i]] = True
						visited2[elem][i] = True
						print adj2[elem][i]
		if is_multi_component == True:
			num_component += 1
	else:
		break

print "number of component is: ", num_component
