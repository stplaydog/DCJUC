import sys

writer = open(sys.argv[3], "w")

deg_reader = open(sys.argv[1])

lines = deg_reader.readlines()
num_v = len(lines)
num_e = 0
num_g = 1
for line in lines:
	e_num = int(line.strip())
	num_e += e_num

writer.write(str(num_v)+" "+str(num_g)+" "+str(num_e*2)+"\n")

pair_reader = open(sys.argv[2])
pairs = pair_reader.readline().strip().split("\t")
for p in pairs:
	p = p.strip().replace("(", "").replace(")", "").split(",")
	writer.write(p[0]+" "+p[1]+" "+"1\n")
	writer.write(p[1]+" "+p[0]+" "+"1\n")
