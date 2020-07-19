import sys

reader = open(sys.argv[1])
writer = open(sys.argv[2], "w")

lines = reader.readlines()


num_e = 0
num_v = 0
for line in lines:
	if line[0] == 'p':
		items = line.strip().split(" ")
		num_v = items[2]
		num_e = int(items[3])*2
		writer.write(str(num_v)+" 1 "+str(num_e)+"\n")
	if line[0] == 'e':
		items = line.strip().split(" ")
		writer.write(items[1]+" "+items[2]+" 1\n")
		writer.write(items[2]+" "+items[1]+" 1\n")

		
