import sys


reader = open("out1")
lines = reader.readlines()

time1=[]
time2=[]
time3=[]
time4=[]
time5=[]
time6=[]

for line in lines:
	if line.find("time") != -1:
		content = line.split(" ")
		time1.append(float(content[2]))
		time2.append(float(content[4]))
		time3.append(float(content[5]))
		time4.append(float(content[6]))
		time5.append(float(content[7]))
		time6.append(float(content[8]))

avg1=[]
avg2=[]
avg3=[]
avg4=[]
avg5=[]
avg6=[]

for i in range(4):
	total1=0;
	total2=0;
	total3=0;
	total4=0;
	total5=0;
	total6=0;
	for j in range(100):
		total1 += time1[i*100+j]
		total2 += time2[i*100+j]
		total3 += time3[i*100+j]
		total4 += time4[i*100+j]
		total5 += time5[i*100+j]
		total6 += time6[i*100+j]
	avg1.append(total1/100)
	avg2.append(total2/100)
	avg3.append(total3/100)
	avg4.append(total4/100)
	avg5.append(total5/100)
	avg6.append(total6/100)

print "get       ","to       ","add       ","from     "
print avg1[0],avg2[0],avg3[0],avg4[0],avg5[0],avg6[0]	
print avg1[1],avg2[1],avg3[1],avg4[1],avg5[1],avg6[1]	
print avg1[2],avg2[2],avg3[2],avg4[2],avg5[2],avg6[2]	
print avg1[3],avg2[3],avg3[3],avg4[3],avg5[3],avg6[3]	
