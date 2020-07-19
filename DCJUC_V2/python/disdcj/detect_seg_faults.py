import sys

reader = open(sys.argv[1])

lines = reader.readlines()

turn_on=0
curline=''
for line in lines:
    if line.find("processing:") != -1:
        turn_on += 1
        if turn_on > 1: #error here
            print curline.strip()
            turn_on = 1
        curline=line
    if line.find("num_extern_comp") != -1:
        turn_on -= 1
