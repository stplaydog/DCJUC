import sys
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

reader = open(sys.argv[1])

lines = reader.readlines()

genes = []

e_score = []
for i in range(10):
        e_score.append([])
o_score = []
for i in range(10):
        o_score.append([])
for i in range(len(lines)):
	if lines[i].find("_") != -1 and lines[i].find("buck") == -1 and lines[i].find("upper") == -1:
		i_now = i+1
		which = -1
		for j in range(i_now, len(lines)):
			if lines[j].find("distance") != -1: 
				which = j
				break
		if which != -1:
			pos = int(float(lines[i].split("_")[0])*10-1)
			o_score[pos].append(int(lines[which].split(" ")[0]))
			e_score[pos].append(int(lines[which].split(" ")[3]))

	

#for line in lines:
#        if line.find("distance") != -1:
#                genes.append(line.strip())
#
#e_score = []
#for i in range(10):
#        e_score.append([])
#o_score = []
#for i in range(10):
#        o_score.append([])
#
#for i in range(10):
#	for j in range(10):
#        	g = genes[i*10+j].split(" ")
#        	s = int(g[3])
#        	e_score[i].append(s)
#        	s = int(g[0])
#        	o_score[i].append(s)

x = [40, 60, 80, 100, 120, 140, 160, 180, 200, 220] 
for i in range(len(x)):
	x[i] = x[i]+int(sys.argv[3]) 
print x
y = []
erry_up = []
erry_low = []    
for i in range(len(e_score)):
        y.append(np.mean(e_score[i]))
	standard_deviation = np.std(e_score[i])
plt.errorbar(x, y, yerr=standard_deviation, color='k', linewidth=2, linestyle="-", label="DCJ-Indel-Exem-Matc")

x = [40, 60, 80, 100, 120, 140, 160, 180, 200, 220]
for i in range(len(x)):
	x[i] = x[i]+int(sys.argv[3])
y = []
erry_up = []
erry_low = []
for i in range(len(o_score)):
        y.append(np.mean(o_score[i]))
	standard_deviation = np.std(o_score[i])
plt.errorbar(x, y, yerr=standard_deviation, color='k', linewidth=2, linestyle=":", label="DCJ-Indel-Matc")


plt.plot([0, 220], [0, 220], color='k', linestyle='-', linewidth=1, label="perfect estimator")
plt.legend( loc='upper left', numpoints = 1, prop={'size':18} )
plt.tick_params(labelsize=18)
plt.ylabel("estimated distance", fontsize=25)
plt.xlabel("actual # evolutionary events", fontsize=25)
plt.savefig(sys.argv[2])
