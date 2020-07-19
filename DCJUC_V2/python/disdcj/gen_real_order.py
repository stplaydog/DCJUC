##########################################################
# This file is used to process ensenble data, and produce
# two genomes with only duplications of the same size.
##########################################################

import sys
import collections


##########################################################
# sort code based on their start position
##########################################################
def sort_code(code):
    for i in range(2):
        sorted(code[i], key=code[i].get)
        for k in code[i]:
            sorted(code[i][k], key=code[i][k].get)
            code[i][k] = collections.OrderedDict(sorted(code[i][k].items()))



##########################################################
# process genomes
##########################################################
def process_genome (file1, file2):
    # read lines
    reader = []
    lines = []
    reader1 = open(file1)
    reader2 = open(file2)
    lines1 = reader1.readlines()
    lines2 = reader2.readlines()
    reader.append(reader1)
    reader.append(reader2)
    lines.append(lines1)
    lines.append(lines2)
    # prepare some basic data structures
    family = {}
    f_occur = [{}, {}]
    code = [{}, {}]
    o_code = [{}, {}]
    family_count = 1
    # start processing    
    for i in range(2):
        last = ""
        last_on = False
        for line in lines[i]:
            if  line.find("ENS") != -1:
                # process the line and get different item information
                line = line.replace("MT", "100").replace("X", "101").replace("Y", "102").replace("W", "103").replace("Z", "104")
                items = line.split(",")
                gene_id = items[0]
                transcript_id = items[1]
                chrom_id = int(items[2])
                strand = int(items[3])
                start = int(items[4])
                family_id = items[6]
                # find a new gene, process it
                if gene_id != last:
                    last = gene_id
                    f_id =0
                    last_on = False
                if last_on == False and len(family_id)>10:
                    last_on = True
                    if family_id in family:
                        f_id = family[family_id]
                        if f_id in f_occur[i]:
                            f_occur[i][f_id] += 1
                        else:
                            f_occur[i][f_id] = 1
                    else:
                        family[family_id] = family_count
                        f_id = family_count
                        f_occur[i][f_id] = 1
                        family_count += 1
                    ##add a new chromosome if there is none
                    if chrom_id not in code[i]:
                        code[i][chrom_id] = {}
                        o_code[i][chrom_id] = {}
                    code[i][chrom_id][start] = f_id*strand
                    o_code[i][chrom_id][start] = gene_id
    sort_code(code)
    sort_code(o_code)
    return code, o_code, f_occur


##########################################################
# write out results
# mode could be:
# 'balanced' : only gene families with the same number of occurances
# 'unbalanced' : gene families no need to have the same number of occurances
##########################################################
def write_out(f, mode, f1, f2, code, o_code, f_occur):
    writer = open(f, "w") 
    wri = []
    writer1 = open(f1, "w") 
    writer2 = open(f2, "w") 
    wri.append(writer1);
    wri.append(writer2);
    for i in range(2):
        writer.write(">genome"+str(i)+"\n")
        another = 1 if i==0 else 0
        for k in code[i]:
            chrom_written = 0
            for j in code[i][k]:
                g_id = abs(code[i][k][j])
                if mode == 'balanced':
                    if g_id in f_occur[i] and g_id in f_occur[another] and f_occur[i][g_id]==f_occur[another][g_id]:
                        chrom_written += 1
                        writer.write(str(code[i][k][j])+" ")
                        wri[i].write(str(o_code[i][k][j])+" "+str(code[i][k][j])+" "+str(k)+" 1 \n")
                elif mode == 'unbalanced':
                    chrom_written += 1
                    writer.write(str(code[i][k][j])+" ")
                    wri[i].write(str(o_code[i][k][j])+" "+k+" "+code[i][k][j]+"\n")
            if chrom_written > 0:
                writer.write("\n")

##########################################################
# this is the main entrance
##########################################################
code, o_code, f_occur = process_genome (sys.argv[1], sys.argv[2])
write_out(sys.argv[4], sys.argv[3], sys.argv[5], sys.argv[6], code, o_code, f_occur)
