#############################
####set up configurations####
#############################
theta=0.1               # the mutation rate
gamma=0.0               # the indel rate
phi=0.1                 # the duplication rate
seqlen=1000             # the sequence length
evo_mode='dual'         # mutation on both of the genomes 'dual', 
                        # on a single genome 'sing'
str_mode='balanced'     # equal num of duplications 'balanced',
                        # unequal num of duplications 'unbalanced'
config=${seqlen}_${theta}_${gamma}_${phi}_${evo_mode}_${str_mode}

#############################
####set up folder information####
#############################
# base folder might be dynamically changed
base_folder=/Users/zhaomingyin/gitlocal/DCJUC/DCJUC_V2/
# folder to store sequence
seq_folder=${base_folder}data/dist/seq/${config}/
# folder to store graph
graph_folder=${base_folder}data/dist/graph/${config}/
# folder to store bijection result
bijection_folder=${base_folder}data/dist/bijection/${config}/

#############################
####set up executable information####
#############################
# folder to call shell scripts
shell_folder=${base_folder}shell/disdcj/
# folder to call python scripts
python_folder=${base_folder}python/disdcj/
# python code to generate gene orders
py_gen_o=${python_folder}gen_2_orders.py
#python code to transfer orders to graphs
py_o2g=${python_folder}o_2_g.py 

#############################
####make the directory for seq and graph
#############################
# check the existance of sequence folder, 
# if no create one
if [ ! -d "${seq_folder}" ]; then
	mkdir ${seq_folder}
fi
# check the existance of graph folder, 
# if no create one
if [ ! -d "${graph_folder}" ]; then
	mkdir ${graph_folder}
fi
# check the existance of bijection folder, 
# if no create one
if [ ! -d "${bijection_folder}" ]; then
	mkdir ${bijection_folder}
fi

#############################
#generate seqeuence first
#############################
python ${py_gen_o} ${seqlen} ${theta} ${gamma} ${phi} ${seq_folder} ${evo_mode} ${str_mode} ${bijection_folder}


#############################
#then transfer orders into graph
#############################
for file in $(ls ${seq_folder})
do
	echo compute $file
	python ${py_o2g} ${seq_folder}$file ${graph_folder}$file
done
