import sys

def read_bijection(f):
    # dict: {gene_id : [[which_occur, pre_gene, after_gene], ...[]]}
    d = [{}, {}]
    d_idx = -1

    # read file now
    reader = open(f, "r")
    lines = reader.readlines()
    for line in lines:
        if line.find(">") == -1:
            d_idx += 1
            genes = line.split(" ")
            for i in range(len(genes)):
                gene = genes[i]
                gene_id = gene.split(":")[0]
                which_occur = gene.split(":")[1]

                # if has pre and post
                if i > 0 and i < len(genes) -1:
                    pre_gene_id = genes[i-1].split(":")[0]
                    post_gene_id = genes[i+1].split(":")[0]
                # if only has pre
                elif i == len(genes) - 1:
                    pre_gene_id = genes[i-1].split(":")[0]
                    post_gene_id = "-1"
                # if only has post
                elif i == 0: 
                    pre_gene_id = "-1"
                    post_gene_id = genes[i+1].split(":")[0]

                # assign
                if gene_id in d:
                    d[d_idx][gene_id].append([which_occur, pre_gene_id, post_gene_id])
                else:
                    d[d_idx][gene_id] = [which_occur, pre_gene_id, post_gene_id]
    return d

 def cal_bij_rate(f1, f2)
    # get the dict for benchmark file
    d = read_bijection(f1)

    # then read the calculated bijection
    reader = open(f2, "r")
    lines = reader.readlines()
    total = 0.0
    correct = 0.0
    for line in lines:
        total += 1
        # judge the correctness here
        gene_one_id = line.split(" ")[0].split(":")[0]
        gene_one_nei = line.split(" ")[0].split(":")[1]
        gene_two_id = line.split(" ")[1].split(":")[0]
        gene_two_nei = line.split(" ")[1].split(":")[1]
        
        # get the occur
        occur = 0
        for item in d[0][gene_one_id]:
            if item[1] == gene_one_nei or item[2] == gene_one_nei:
                occur = item[0]
        
        # find another genome in the benchmark
        for item in d[1][gene_two_id]:
            if item[0] == occur and (item[1] == gene_two_nei or item[2] == gene_two_nei):
                correct += 1
        correct_rate = correct / total;
        print correct_rate



 # main entrance
 cal_bij_rate(sys.argv[1], sys.argv[2])

