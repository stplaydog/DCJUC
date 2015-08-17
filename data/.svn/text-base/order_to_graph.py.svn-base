import sys

reader = open(sys.argv[1])
writer = open(sys.argv[2], "w")

name1 = reader.readline();
seq1 = reader.readline().split(" ");
name2 = reader.readline();
seq2 = reader.readline().split(" ");
name3 = reader.readline();
seq3 = reader.readline().split(" ");

count = 0

ind = [0]*65534
total_gene=0
max =0

for i in range(len(seq1)):
	num = abs(int(seq1[i]))
	if ind[num]==0:
		total_gene+=1
		ind[num] = 1
	if num>max:
		max = num
for i in range(len(seq2)):
	num = abs(int(seq2[i]))
	if ind[num]==0:
		total_gene+=1
		ind[num] = 1
	if num>max:
		max = num
for i in range(len(seq3)):
	num = abs(int(seq3[i]))
	if ind[num]==0:
		total_gene+=1
		ind[num] = 1
	if num>max:
		max = num
if max > total_gene:
	total_gene = max
total_edge = len(seq1)*2+len(seq2)*2+len(seq3)*2
writer.write(str(total_gene*2)+" "+str(total_edge)+"\n")

for i in range(len(seq1)):
	onum = int(seq1[i])
	num = abs(int(seq1[i]))
	prevonum = int(seq1[i-1]) if (i-1)>=0 else 32767
	prevnum = abs(int(seq1[i-1])) if (i-1)>=0 else 32767
	head = (num-1)*2 if onum>0 else (num-1)*2+1
	tail = (num-1)*2+1 if onum>0 else (num-1)*2
	prev_head = (prevnum-1)*2 if prevonum>0 else (prevnum-1)*2+1
	prev_tail = (prevnum-1)*2+1 if prevonum>0 else (prevnum-1)*2
	if head>prev_tail:
		writer.write(str(prev_tail)+" "+str(head)+" 1\n") 
		if head!=65533:
			writer.write(str(head)+" "+str(prev_tail)+" 1\n") 
	else:
		writer.write(str(head)+" "+str(prev_tail)+" 1\n")
		if prev_tail != 65533:
			writer.write(str(prev_tail)+" "+str(head)+" 1\n")
	if i == len(seq1)-1:
		writer.write(str(tail)+" 65533"+" 1\n")

for i in range(len(seq2)):
	onum = int(seq2[i])
	num = abs(int(seq2[i]))
	prevonum = int(seq2[i-1]) if (i-1)>=0 else 32767
	prevnum = abs(int(seq2[i-1])) if (i-1)>=0 else 32767
	head = (num-1)*2 if onum>0 else (num-1)*2+1
	tail = (num-1)*2+1 if onum>0 else (num-1)*2
	prev_head = (prevnum-1)*2 if prevonum>0 else (prevnum-1)*2+1
	prev_tail = (prevnum-1)*2+1 if prevonum>0 else (prevnum-1)*2
	if head>prev_tail:
		writer.write(str(prev_tail)+" "+str(head)+" 2\n") 
		if head!=65533:
			writer.write(str(head)+" "+str(prev_tail)+" 2\n") 
	else:
		writer.write(str(head)+" "+str(prev_tail)+" 2\n")
		if prev_tail != 65533:
			writer.write(str(prev_tail)+" "+str(head)+" 2\n")
	if i == len(seq2)-1:
		writer.write(str(tail)+" 65533"+" 2\n")

for i in range(len(seq3)):
	onum = int(seq3[i])
	num = abs(int(seq3[i]))
	prevonum = int(seq3[i-1]) if (i-1)>=0 else 32767
	prevnum = abs(int(seq3[i-1])) if (i-1)>=0 else 32767
	head = (num-1)*2 if onum>0 else (num-1)*2+1
	tail = (num-1)*2+1 if onum>0 else (num-1)*2
	prev_head = (prevnum-1)*2 if prevonum>0 else (prevnum-1)*2+1
	prev_tail = (prevnum-1)*2+1 if prevonum>0 else (prevnum-1)*2
	if head>prev_tail:
		writer.write(str(prev_tail)+" "+str(head)+" 3\n") 
		if head != 65533:
			writer.write(str(head)+" "+str(prev_tail)+" 3\n")
	else:
		writer.write(str(head)+" "+str(prev_tail)+" 3\n")
		if prev_tail != 65533:
			writer.write(str(prev_tail)+" "+str(head)+" 3\n")
	if i == len(seq3)-1:
		writer.write(str(tail)+" 65533"+" 3\n")

