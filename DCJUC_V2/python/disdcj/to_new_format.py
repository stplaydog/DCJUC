import sys

reader = open(sys.argv[1])
writer = open(sys.argv[2], "w") 

lines = reader.readlines()

genes1 = lines[1].strip().split(" ")
genes2 = lines[3].strip().split(" ")

num_gf = 0
max = 0
for gene in genes1:
	if int(gene) > max:
		max = int(gene)
for gene in genes2:
	if int(gene) > max:
		max = int(gene)
num_gf = max

num_g1 = len(genes1)
num_g2 = len(genes2)

writer.write(str(num_gf)+" "+str(num_g1)+" "+str(num_g2)+"\n")
for line in lines:
	writer.write(line)
	
		
	
