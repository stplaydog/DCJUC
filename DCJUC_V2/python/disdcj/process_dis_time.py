import sys
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

reader = open(sys.argv[1])
lines = reader.readlines()


time_score = []
for i in range(10):
        time_score.append([])
for i in range(len(lines)):
	if lines[i].find("_") != -1 and lines[i].find("buck") == -1 and lines[i].find("upper") == -1:
		i_now = i+1
		which = -1
		is_good = False
		for j in range(i_now, len(lines)):
			if lines[j].find("distance") != -1:
				is_good = True
			if lines[j].find("time") != -1 and is_good == True: 
				which = j
				break
		if which != -1:
			pos = int(float(lines[i].split("_")[0])*10-1)
			time_score[pos].append(float(lines[which].split(" ")[1]))
t=[]
for m in time_score:
	t.append(np.mean(m))

x=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
plt.plot(x, t, color='k', linestyle='-', linewidth=1, label="BnB time")
#plt.legend( loc='upper left', numpoints = 1, prop={'size':18} )
plt.tick_params(labelsize=18)
plt.ylabel("process time (DCJ-Indel-Matc)", fontsize=25)
plt.xlabel("evolutionary rate", fontsize=25)
plt.savefig(sys.argv[2])
