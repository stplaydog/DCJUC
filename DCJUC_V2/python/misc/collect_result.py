import sys
from os import listdir
from os.path import isfile, join

dis_stats = [[],[],[],[],[],[],[],[],[],[],[]]

files = [ f for f in listdir(sys.argv[1]) if isfile(join(sys.argv[1],f)) ]

for f in files:
	reader = open(sys.argv[1]+f)
	pos = int(float(f.split("_")[0])*10)-1
	lines = reader.readlines()
	num_comp = 0
	for line in lines:
		items = line.split(" ")
		num_comp += int(float(items[1].strip()))
	dis_stats[pos].append(num_comp)

for i in range(len(dis_stats)):
	for j in range(len(dis_stats[i])):
		dis_stats[i][j] = int(sys.argv[2])-dis_stats[i][j]
print dis_stats

