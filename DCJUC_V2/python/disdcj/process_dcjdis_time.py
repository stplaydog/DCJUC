import sys
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

#read file into lines
reader_optkit_exem = open(sys.argv[1])
reader_optkit_matc = open(sys.argv[2])
reader_gredo = open(sys.argv[3])
writer = open(sys.argv[4], "w")
type = sys.argv[5]

lines_optkit_exem = reader_optkit_exem.readlines()
lines_optkit_matc = reader_optkit_matc.readlines()
lines_gredo = reader_gredo.readlines()


#initialize containers for data
time_optkit_exem = {}
time_optkit_matc = {}
time_gredo       = {}
if type == "sim":
    time_optkit_exem = {1:[],2:[],3:[],4:[],5:[],6:[],7:[],8:[],9:[],10:[]}
    time_optkit_matc = {1:[],2:[],3:[],4:[],5:[],6:[],7:[],8:[],9:[],10:[]}
    time_gredo       = {1:[],2:[],3:[],4:[],5:[],6:[],7:[],8:[],9:[],10:[]}
elif type == "real":
    time_optkit_exem = {}
    time_optkit_matc = {}
    time_gredo       = {}

#process name inforamtion
def process_name(line):
    name = 0
    if type == "sim":
        name = int(float(line.replace("processing ", "").split("_")[0])*10)
    elif type == "real":
        name = line.replace("processing ", "").replace(" ", "_").strip()
    return name

#add into data structure
def add_to_dic(name, val, time_dic):
    if type == "sim":
        time_dic[name].append(val)
    elif type == "real":
        time_dic[name] = val
    
#function for processing data
def process_line(lines, time_dic):
    name = 0
    val = 0.0
    for line in lines:
        if line.find("processing ") != -1 and line.find("finish") == -1:
            name = process_name(line)
        if line.find("time ") != -1 and line.find("Total") == -1:
            val = float(line.replace("time ", ""))
            add_to_dic(name, val, time_dic)


#extract optkit_exem information
process_line(lines_optkit_exem, time_optkit_exem)
#extract optkit_matc information
process_line(lines_optkit_matc, time_optkit_matc)
#extract gredo information
process_line(lines_gredo, time_gredo)


#process data using python 
if type == "sim":
    writer.write("mrate, time, algo\n")
    for key in time_optkit_exem:
        for i in range(10):
            writer.write(str(key)+", "+str(time_optkit_exem[key][i])+", BnB_Exem\n")
            writer.write(str(key)+", "+str(time_optkit_matc[key][i])+", BnB_Matc\n")
            writer.write(str(key)+", "+str(time_gredo[key][i])+", GREDO\n")
elif type == "real":
    writer.write("mrate, time, algo\n")
    for key in time_optkit_exem:
        writer.write(str(key)+", "+str(time_optkit_exem[key])+", BnB_Exem\n")
        writer.write(str(key)+", "+str(time_optkit_matc[key])+", BnB_Matc\n")
        writer.write(str(key)+", "+str(time_gredo[key])+", GREDO\n")
