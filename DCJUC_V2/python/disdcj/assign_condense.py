import sys
import os
from graph import read_graph
from os import listdir
from os.path import isfile, join
from small_CC import compute_small_cc
from graph import vis_graph
from graph import delete_reg_CC
from graph import rename_adj

def check_cycle_two(adj, key, c1, c2):
    to = adj[c1][key][0]
    is_cycle_two = False
    if to in adj[c2]:
        for i in range(len(adj[c2][to])):
            if adj[c2][to][i] == key: 
                # which is a cycle of two
                is_cycle_two = True
                break
    return is_cycle_two, to

def opt_assin(adj, c1, c2, w_map):
    for key in adj[c1]:
        if len(adj[c1][key]) == 1:
            # check if this vertex satisfy the 
            # optimization assignment conditions
            is_cycle_two, to = check_cycle_two(adj, key, c1, c2)
            #change other vertices
            if is_cycle_two == True:
                # here print the map to the file
                w_map.write(str(key)+" "+str(to)+":"+str(key)+" "+str(to)+"\n")
                # change the graph to reflect the changes
                for c in range(2):
                    for i in range(len(adj[c][to])):
                        toto = adj[c][to][i]
                        if toto != key:
                            tmp = []
                            for j in range(len(adj[c][toto])):
                                if adj[c][toto][j] != to:
                                    tmp.append(adj[c][toto][j])
                            adj[c][toto]  = tmp
                    adj[c][to] = [key]
    return adj

def condense(adj, c1, c2):
    while True:
        find = False
        for key in adj[c1]:
            if len(adj[c1][key])==1 and \
                    key in adj[c2] and \
                    len(adj[c2][key])==1:
                to_c1 = adj[c1][key][0]
                to_c2 = adj[c2][key][0]
                if to_c1 != to_c2 and \
                        len(adj[c1][to_c1])==1 and \
                        to_c2 in adj[c2] and \
                        len(adj[c2][to_c2])==1:
                    if to_c2 in adj[c1] and len(adj[c1][to_c2])==1:
                        to_c2_to = adj[c1][to_c2][0]
                        if to_c2_to in adj[c1] and \
                                len(adj[c1][to_c2_to])==1:
                            find = True
                            # reconnection
                            adj[c1][to_c1][0] = to_c2_to
                            adj[c1][to_c2_to][0] = to_c1
                            adj[c1].pop(key)
                            adj[c2].pop(key)
                            adj[c1].pop(to_c2)
                            adj[c2].pop(to_c2)
                            break
        if find== False:
            break
    return adj


def optimal_assign_n_condense(input_dir, input_file, out_dir, bij_file):
    # prepare the output directory
    small_cc_num = 0
    if os.path.isdir(out_dir) == False:
        os.mkdir(out_dir)
    # open the map file
    w_bij = open(bij_file, 'w')
    # read the graph into adj
    adj = read_graph(input_dir+"/"+input_file)
    # perform the optiaml assignment
    adj = opt_assin(adj, 0, 1, w_bij)
    # condense the graph
    adj = condense(adj, 0, 1)
    adj = condense(adj, 1, 0)
    # perform the optiaml assignment again
    adj = opt_assin(adj, 0, 1, w_bij)
    # vis_graph(adj, "after_condense.dot")
    return adj

def compute_dis_n_output(adj, algo, input_file, out_dir):
    small_cc_num = 0
    input_base = input_file.replace(".gr", "")+"_"
    small_cc_num,count = delete_reg_CC(adj, out_dir, input_base)
    return small_cc_num
    

####################
# real execution part
####################

# create to folder
input_dir = sys.argv[1]
out_dir = sys.argv[2]
opt_tmp_file = sys.argv[3]
algo = sys.argv[4]
vis_folder = sys.argv[5]
bij_file = sys.argv[6]
dict_dir = sys.argv[7]

if os.path.isdir(out_dir) == False:
    os.mkdir(out_dir)
# list all files in the input folder
files = [ f for f in listdir(input_dir) if isfile(join(input_dir,f)) ]
small_cc_num = 0
for f in files:
    adj = optimal_assign_n_condense(input_dir, f, out_dir, bij_file)
    small_cc_num += compute_dis_n_output(adj, algo, f, out_dir)

files = [ f for f in listdir(out_dir) if isfile(join(out_dir,f)) ]
for f in files:
    adj = read_graph(out_dir+"/"+f)
    # when rename, there should be a map file stored
    new_adj = rename_adj(adj, out_dir+"/"+f, dict_dir+"/"+f)
    # vis_graph(new_adj, "after_rename.dot")    

# append the number of small irregular component to the end of the file
writer = open(opt_tmp_file, 'a')
writer.write("num_small_irr_CC: "+ str(small_cc_num)+"\n")
