import sys

reader = open(sys.argv[1])
writer = open(sys.argv[2], "w")

names = []
genes = []

lines = reader.readlines()

for line in lines:
	if line.find(">") == -1:
		genes.append(line)
	else:
		names.append(line)

#find the gene_num
gene_num=0
chr_num = len(genes)
edge_num = 0
for i in range(len(genes)):
	g = genes[i].split(" ")
	for j in range(len(g)):
		if abs(int(g[j])) > gene_num:
			gene_num = abs(int(g[j]))
	edge_num = edge_num + len(g)*2
#print head here
writer.write(str(gene_num*2)+" "+str(chr_num)+" "+str(edge_num)+"\n")

#print edges
for i in range(len(genes)):
	g = genes[i].split(" ")
	for j in range(len(g)):
		if j == 0:
			head = (abs(int(g[j]))-1)*2 if int(g[j])>0 else (abs(int(g[j]))-1)*2+1
			writer.write(str(head) + " " + "65533 " + str(i+1)+"\n")
		elif j == len(g)-1:
			head = (abs(int(g[j]))-1)*2 if int(g[j])>0 else (abs(int(g[j]))-1)*2+1
			tail = (abs(int(g[j]))-1)*2 if int(g[j])<0 else (abs(int(g[j]))-1)*2+1
			pre_tail = (abs(int(g[j-1]))-1)*2 if int(g[j-1])<0 else (abs(int(g[j-1]))-1)*2+1
			writer.write(str(pre_tail) + " " + str(head) + " "+ str(i+1)+"\n")
			writer.write(str(head) + " " + str(pre_tail) + " "+ str(i+1)+"\n")
			writer.write(str(tail) + " " + "65533 " + str(i+1)+"\n")
		else:
			head = (abs(int(g[j]))-1)*2 if int(g[j])>0 else (abs(int(g[j]))-1)*2+1
			tail = (abs(int(g[j-1]))-1)*2 if int(g[j-1])<0 else (abs(int(g[j-1]))-1)*2+1
			writer.write(str(tail) + " " + str(head) + " "+ str(i+1)+"\n")
			writer.write(str(head) + " " + str(tail) + " "+ str(i+1)+"\n")
