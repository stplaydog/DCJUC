import sys

reader = open(sys.argv[1])
writer1 = open(sys.argv[2]+"_1", "w")
writer2 = open(sys.argv[2]+"_2", "w")

lines = reader.readlines()

line_count=0
for line in lines:
	if line.find(">") == -1:
		unique=0
		items = line.strip().split(" ")
		line_count +=1
		for gene in items:
			#only circular chromosome are considered
			if(line_count==1):
				writer1.write("GENE"+str(unique)+" "+gene+" "+"CHR1 "+"2\n")
			else:
				writer2.write("GENE"+str(unique)+" "+gene+" "+"CHR1 "+"2\n")
			unique+=1
