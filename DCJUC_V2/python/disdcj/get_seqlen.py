import sys

reader = open(sys.argv[1])
line =reader.readline() 
print line.split(" ")[0]
