
def determin_vet_type(adj):
	vet_type = {}
	PI_OPEN = 0
	GA_OPEN = 1
	PA_OPEN = 2
	CYCLE = 4
	CAP = 65533
	#for color one
	for key in adj[0]:
		if len(adj[0][key])>0 and adj[0][key][0] != CAP and key in adj[1]:
			vet_type[key] = CYCLE
		elif len(adj[0][key])>0 and adj[0][key][0] ==CAP and key in adj[1]:
			vet_type[key] = PA_OPEN
		elif len(adj[0][key])>0 and adj[0][key][0] != CAP and key not in adj[1]:
			vet_type[key] = PI_OPEN
	#for color two
	for key in adj[1]:
		if len(adj[1][key])>0 and adj[1][key][0] != CAP and key in adj[0]:
			vet_type[key] = CYCLE
		elif len(adj[1][key])>0 and adj[1][key][0] ==CAP and key in adj[0]:
			vet_type[key] = PA_OPEN
		elif len(adj[1][key])>0 and adj[1][key][0] != CAP and key not in adj[0]:
			vet_type[key] = GA_OPEN
	return vet_type
	

#implementation of phillip's method 
def determine_cc_type(adj, vet_type, source):
	if (source not in adj[0] or len(adj[0][source])==0) and (source not in adj[1] or len(adj[1][source])==0):
		return 1
	#define some variables
	PI_OPEN = 0
	GA_OPEN = 1
	PA_OPEN = 2
	CYCLE = 4
	CAP = 65533
	#start computing
	visited = {}
	cc_type = 0.0
	start_type = vet_type[source]
	end_type = -1
	visited[source] = True
	n = source
	#because we only examine cycles or paths
	#just visit next by the first entry
	while True:
		if n in adj[0] and n not in visited:
			n = adj[0][n][0]
			visited[n] = True
		elif n in adj[1] and n not in visited:
			n = adj[1][n][0]
			visited[n] = True
		else:
			end_type = vet_type[n]
			break
	#determine the type now
	if start_type == CYCLE:
		cc_type = 1.0
	elif start_type == end_type:
		cc_type = 1.0
	else:
		cc_type = 0.5
		
	return cc_type
