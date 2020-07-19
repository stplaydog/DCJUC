mode_seq=1
mode_par=2
mode_tpar=3

true=0
false=1

#./optec data/dist/seqnew/0.1_0.05_0.05_0 1
#./optec --dcj --input_file data/dist/seqnew/0.1_0.05_0.05_0 1


#./optec --knap --input_file data/knap/p07 --p_mode $mode_seq --num_t 1
#./optec --knap --input_file data/knap/p07 --p_mode $mode_par --num_t 64
#./optec --knap --input_file data/knap/p07 --p_mode $mode_tpar --num_t 128


#./optec --dcj --input_file data/dist/seqnew/1.0_0.05_0.1_1 --p_mode $mode_seq --num_t 1
#./optec --dcj --input_file data/dist/seqnew/1.0_0.05_0.1_1 --p_mode $mode_par --num_t 16
#./optec --dcj --input_file data/dist/seqnew/1.0_0.05_0.1_1 --p_mode $mode_tpar --num_t 16

#./optec --median --input_file data/med/50_0.1_0.1_0.1_2 --p_mode $mode_par --num_t 2 --heu_level 2 --term_move 3 --is_opt $false 
#./optec --median --inputfile data/med/50_0.1_0.1_0.1_0 --p_mode $mode_par --num_t 2 --heu_level 2 --term_move 3 --is_opt $true

#4/4/2014
#./optec --median --input_file data/med/50_0.1_0.1_0.1_1 --p_mode $mode_par --num_t 2 --heu_level 2 --term_move 3 --is_opt $false 
#./optec --median --input_file /home/leonxin/optec-code/data/med/50_0.1_0.05_0.05_0 --p_mode 2 --num_t 2 --heu_level 2 --term_move 3 --is_opt 1

#5/8/2014
#python shell/irrCC.py human_mouse.gr human_mouse.regc
#python shell/assign_condense.py human_mouse 

#####5/16/2014######
gamma=0.1
phi=0.1
theta=0.1
seqlen=200
config=${seqlen}_${theta}_${gamma}_${phi}
folder_now=/Users/zhaomingyin/Dropbox/research/code/optec-code/
#folder_now=/home/users/leonxin/Dropbox/research/code/optec-code/
graph_folder=${folder_now}data/dist/graph/${config}/
shell_folder=${folder_now}shell/
irrCC=irrCC.py
irrCC_folder=${folder_now}data/dist/irrCC/${config}/
opt_tmp_result_folder=${folder_now}data/dist/opt_tmp_result/${config}/
assign_condense=assign_condense.py
assign_folder=${folder_now}data/dist/assign/${config}/
algo=exem
#file=1.0_0.1_0.1_0
file=0.7_0.1_0.1_8
vis_folder=${folder_now}data/dist/vis/${config}/
vis_file=${vis_folder}${file}

#make directories
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

#remove regular CCs	
python $shell_folder$irrCC $graph_folder$file $irrCC_folder$file $opt_tmp_result_folder$file ${vis_file}
#perform opt assignment and condense
python  $shell_folder$assign_condense $irrCC_folder$file $assign_folder$file $opt_tmp_result_folder$file $algo ${vis_file}
#visualization
#./shell/batch_vis_irrCC.sh


#####5/18/2014######
debug=0
gamma=0.1
phi=0.1
theta=0.1
seqlen=200
config=${seqlen}_${theta}_${gamma}_${phi}
folder_now=/Users/zhaomingyin/Dropbox/research/code/optec-code/
#folder_now=/home/users/leonxin/Dropbox/research/code/optec-code/
assign_folder=${folder_now}data/dist/assign/${config}/
dis_mode=2
exe=optec
parallel_mode=1

#rm correct (unneccessary) files for debugging
rm ${assign_folder}${file}/8_1.gr
rm ${assign_folder}${file}/2_8.gr
rm ${assign_folder}${file}/2_22.gr
rm ${assign_folder}${file}/2_15.gr
rm ${assign_folder}${file}/2_11.gr

#execute the algorithm
for f in $(ls ${assign_folder}${file})
do
	if [ "$debug" -eq "1" ]
	then #enter debug mode
		echo lldb -- ${folder_now}${exe} --dis --dis_mode ${dis_mode} --input_file ${assign_folder}${file}/$f --p_mode $parallel_mode --opt_file $opt_tmp_result_folder$file --seq_len $seqlen 
	else #enter execution mode
		echo processing: ${assign_folder}${file}/$f dismode: ${dis_mode}
		${folder_now}${exe} --dis --dis_mode ${dis_mode} --input_file ${assign_folder}${file}/$f --p_mode $parallel_mode --opt_file $opt_tmp_result_folder$file --seq_len $seqlen 
	fi
done
#call python function to collect all results

#visualization
${folder_now}vis.sh
