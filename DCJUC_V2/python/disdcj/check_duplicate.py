import sys

reader = open(sys.argv[1])

lines = reader.readlines()

dic = {}
for line in lines:
	if line.strip() not in dic:
		dic[line.strip()] = True
	else:
		print 'there is duplication detected! which is:', line.strip()
		break
print 'there is no duplication!'
