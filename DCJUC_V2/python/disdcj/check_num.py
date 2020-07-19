import sys

reader = open(sys.argv[1])
lines = reader.readlines()

genes = lines[1].split(" ")

num_dup=0
for i in range(1, len(genes)):
	if genes[i] == genes[i-1]:
		num_dup += 1

print 'dup:', num_dup

genes1 = lines[3].strip().split(" ")
gene_exist = {}
for gene in genes1:
	gene_exist[abs(int(gene))] = True

num_del = 0
for i in range(1, 1000):
	if i not in gene_exist:
		num_del += 1
		

print 'del:', len(gene_exist), num_del 
