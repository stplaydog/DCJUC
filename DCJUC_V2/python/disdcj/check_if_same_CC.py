import sys
from graph import read_graph

adj = read_graph(sys.argv[1])

for key in adj[0]:
	idx = int(key)-1 if int(key)%2==1 else int(key)+1
	if len(adj[0][key])>1 and str(idx) not in adj[0] and str(idx) not in adj[1]:
		print "yes not in the same component", key
print "no, they are all in the same component"
