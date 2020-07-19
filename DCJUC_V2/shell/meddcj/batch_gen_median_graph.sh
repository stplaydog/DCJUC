#############################
####set up configurations####
#############################
theta=0.1               # the mutation rate
gamma=0.05              # the indel rate
phi=0.05                # the duplication rate
seqlen=100              # the sequence length
str_mode='unbalanced'   # equal num of duplications 'balanced',
                        # unequal num of duplications 'unbalanced'
config=${seqlen}_${theta}_${gamma}_${phi}_${str_mode}

#############################
####set up folder information####
#############################
# base folder might be dynamically changed
base_folder=/scratch/zyin/optkit/optkit/
# folder to store sequence
seq_folder=${base_folder}data/median/seq/${config}/
# folder to store graph
graph_folder=${base_folder}data/median/graph/${config}/

#############################
####set up executable information####
#############################
# folder to call shell scripts
shell_folder=${base_folder}shell/meddcj/
# folder to call python scripts
python_folder=${base_folder}python/meddcj/
# python code to generate gene orders and graphs
py_gen_m=${python_folder}gen_med.py

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

#############################
#generate seqeuence and graph 
#############################
python ${py_gen_m} ${seqlen} ${theta} ${gamma} ${phi} ${seq_folder} ${str_mode} ${graph_folder}

