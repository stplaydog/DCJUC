import sys
import os

from CC_type import determin_vet_type
from CC_type import determine_cc_type


def remove_duplicate(adj):
    tmp = []
    for i in range(len(adj)):
        dup = False
        for j in range(i+1, len(adj)):
            if adj[i] == adj[j]:
                dup = True
        if dup == False:
            tmp.append(adj[i])
    return tmp

def add_edge(adj, items):
    color = 0 if items[2]=='1' else 1
    for c in range(2):
        cc = 1 if c==0 else 0
        if items[c] not in adj[color]:
            adj[color][items[c]] = [items[cc]]
        else:
            adj[color][items[c]].append(items[cc])

def read_graph(input_file):
    reader = open(input_file)
    lines = reader.readlines()
    adj=[{},{}]
    for i in range(1, len(lines)):
        items = lines[i].strip().split(" ")
        add_edge(adj, items)
    
    for c in range(2):
        for key in adj[c]:
            adj[c][key] = remove_duplicate(adj[c][key])
    return adj

def vis_graph(adj, f):
    writer = open(f, "w")
    writer.write("graph{\n")
    for c in range(2):
        c_str = "red" if c==0 else "blue"
        for key in adj[c]:
            for i in range(len(adj[c][key])):
                if int(key) > int(adj[c][key][i]):
                    writer.write(key +"--"+ adj[c][key][i] +"[color="+c_str+"];\n")
    writer.write("}\n")

def detect_seed(adj, vet_map):
    seed = -1
    for key in adj[0]:
        if key not in vet_map:
            seed = key
            break
    if seed == -1:
        for key in adj[1]:
            if key not in vet_map:
                seed = key
                break
    return seed
    

##delete the regular CCs############### 
def delete_reg_CC(adj, out_dir, input_base):
    num_reg_CC = 0
    vet_map = {}
    count = 0
    if os.path.isdir(out_dir) == False:
        os.mkdir(out_dir)
    vet_type = determin_vet_type(adj) 
    while True:
        # find the seed first
        seed = detect_seed(adj, vet_map)
        # then using the queue to perform bfs
        # each detected CC is written into a file,
        # if this CC is regualr, delete it, else keep it.
        if seed != -1:
            writer = open(out_dir+"/"+input_base+str(count)+".gr", "w")
            is_multi_component = False
            queue = [seed]
            vet_map[seed] = True
            CC_vet = {}
            CC_vet[seed] = True
            while len(queue) != 0:
                idx = len(queue)-1    
                elem = queue[idx]
                queue.pop(idx)
                #print 'pop',elem
                if (elem in adj[0] and len(adj[0][elem]) > 1 and int(elem) != 65533) or (elem in adj[1] and len(adj[1][elem]) > 1 and int(elem) != 65533):
                    is_multi_component = True
                for color in range(2):
                    c_str = '1\n' if color==0 else "2\n"
                    if elem in adj[color]:
                        for i in range(len(adj[color][elem])):
                            if adj[color][elem][i] not in vet_map:
                                queue.append(adj[color][elem][i])
                                vet_map[adj[color][elem][i]] = True
                                #writer.write(elem+" "+adj[color][elem][i]+" "+c_str)
                                CC_vet[adj[color][elem][i]] = True
            if is_multi_component != True:
                os.remove(out_dir+"/"+input_base+str(count)+".gr")
                #there should be some code to determine which kind of CC it is, but let's just skip this step first
                num_reg_CC+=determine_cc_type(adj, vet_type, seed)
            else: 
                # write the irrCC into file
                # there is no renaming at this time
                writer.write("\n")
                for c in range(2):
                    c_str = "1" if c==0 else "2"
                    for key in adj[c]:
                        if key in CC_vet:
                            for i in range(len(adj[c][key])):
                                if int(key)>int(adj[c][key][i]):
                                    writer.write(key+" "+adj[c][key][i]+" "+c_str+"\n")
                writer.close()
                #vis_after_file = sys.argv[4] + "_"+str(count)+".dot"
                #tmp_adj = read_graph(out_dir+"/"+str(count)+".gr")
                #vis_graph(tmp_adj, vis_after_file)
            
        else:
            break
        count += 1
    return num_reg_CC,count

####this is to rename the graph for cpp file
def rename_adj(adj, f, f_map):
    # get the dic
    new_name = {}
    old_name = {}
    rank = 0
    comb_count = 0;
    for c in range(2):
        for key in sorted(adj[c]):
            if key not in new_name:
                new_name[key] = str(rank)
                old_name[str(rank)] = key
                rank += 1
                #
                if int(key)%2==0:
                    another = str(int(key)+1)
                    if another in adj[0] or another in adj[1]:
                        comb_count += 1
    num_v = rank
    # get a new adj
    new_adj=[{},{}]
    num_e = 0
    for c in range(2):
        for key in sorted(adj[c]):
            for i in range(len(adj[c][key])): 
                fr = new_name[key]
                to = new_name[adj[c][key][i]]
                if fr in new_adj[c]:
                    new_adj[c][fr].append(to) 
                else:
                    new_adj[c][fr] = [to]
                num_e += 1
    # output to the file
    writer = open(f, "w")
    writer.write(str(num_v)+" 2 "+str(num_e)+" "+str(comb_count/2)+"\n")
    for c in range(2):
        c_str = "1" if c==0 else "2"
        for key in new_adj[c]:
            for i in range(len(new_adj[c][key])):
                writer.write(key+" "+new_adj[c][key][i]+" "+c_str+"\n")
    writer.close()
    # write the map file
    writer1 = open(f_map, "w")
    for key in old_name:
        writer1.write(key +" "+old_name[key]+"\n")
    writer1.close()
    return new_adj
