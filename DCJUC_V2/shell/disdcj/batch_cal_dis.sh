########set basic information#############
gamma=0.1
phi=0.1
theta=0.1
seqlen=200
config=${seqlen}_${theta}_${gamma}_${phi}
folder_now=/Users/zhaomingyin/Dropbox/research/code/optec-code/
#folder_now=/home/users/leonxin/Dropbox/research/code/optec-code/

########some folder information and exe information#############
graph_folder=${folder_now}data/dist/graph/${config}/
shell_folder=${folder_now}shell/
irrCC=irrCC.py
irrCC_folder=${folder_now}data/dist/irrCC/${config}/
opt_tmp_result_folder=${folder_now}data/dist/opt_tmp_result/${config}/
assign_condense=assign_condense.py
assign_folder=${folder_now}data/dist/assign/${config}/
algo=exem
vis_folder=${folder_now}data/dist/vis/${config}/
vis_file=${vis_folder}${file}

#########make directories##################
if [ ! -d "${irrCC_folder}" ]; then
        mkdir ${irrCC_folder}
fi

if [ ! -d "${opt_tmp_result_folder}" ]; then
        mkdir ${opt_tmp_result_folder}
fi

if [ ! -d "${assign_folder}" ]; then
        mkdir ${assign_folder}
fi

if [ ! -d "${vis_folder}" ]; then
        mkdir ${vis_folder}
fi

###########remove regular CC and compute distancess#############
parallel_mode=1
dis_mode=2
exe=optec
time=$(TIMEFORMAT=%R; (time sleep 0.1; echo a) 2>&1 > /dev/null)
for file in $(ls $graph_folder)
do
	echo $file
	python $shell_folder$irrCC $graph_folder$file $irrCC_folder$file $opt_tmp_result_folder$file ${vis_file}
	#perform opt assignment and condense
	python  $shell_folder$assign_condense $irrCC_folder$file $assign_folder$file $opt_tmp_result_folder$file $algo ${vis_file}
	#compute distance
	for f in $(ls ${assign_folder}${file})
	do
		${folder_now}${exe} --dis --dis_mode ${dis_mode} --input_file ${assign_folder}${file}/$f --p_mode $parallel_mode --opt_file $opt_tmp_result_folder$file --seq_len $seqlen 
	done
done
time+=+$(TIMEFORMAT=%R; (time sleep 0.1) 2>&1)
bc <<< $time

###########visualization##################
#./shell/batch_vis_irrCC.sh ${config}
