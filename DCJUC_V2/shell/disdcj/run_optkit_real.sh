#!/bin/bash
# remember to input file as the first parameter
if [ $# -ne 3 ]; then 
    echo "Illegal number of parameters, please follow:"
    echo "run_dist_real.sh <folder_now> <seqfile> <seqlen>"
    exit $?
fi
# if running real experiment, remember to remove 
file=$2
debug=0
folder_now=$1
seqlen=$3
config=$4
shell_folder=${folder_now}shell/disdcj/
python_folder=${folder_now}python/disdcj/
data_folder=${folder_now}data/real/optkit/
#folder_now=/home/users/leonxin/Dropbox/research/code/optec-code/
graph_folder=${data_folder}dup_graph/
irrCC=irrCC.py
irrCC_folder=${data_folder}irrCC/
opt_tmp_result_folder=${data_folder}opt_tmp_result/
assign_condense=assign_condense.py
assign_folder=${data_folder}assign/
algo=matc
#algo=exem
vis_folder=${data_folder}vis/
vis_file=${vis_folder}${file}
f_map_folder=${folder_now}data/real/optkit/f_map/
f_map_file=${f_map_folder}${file}
map_folder=${f_map_file}_folder/

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

if [ ! -d "${f_map_folder}" ]; then
        mkdir ${f_map_folder}
fi

if [ ! -d "${map_folder}" ]; then
        mkdir ${map_folder}
fi

#remove regular CCs	
python $python_folder$irrCC $graph_folder$file $irrCC_folder$file $opt_tmp_result_folder$file ${vis_file}
#perform opt assignment and condense
python  $python_folder$assign_condense $irrCC_folder$file $assign_folder$file $opt_tmp_result_folder$file $algo ${vis_file} ${f_map_file} ${map_folder} 
#visualization
#./shell/batch_vis_irrCC.sh


#rm correct (unneccessary) files for debugging
#for f in $(ls ${assign_folder}${file})
#do
#    if [ "${assign_folder}${file}/${f}" != "/Users/zhaomingyin/gitlocal/optkit/data/real/optkit/assign/chicken_mouse/1_9.gr" ]
#    then
#        rm ${assign_folder}${file}/${f}
#    else
#        echo $f
#    fi
#done

#execute the algorithm
dis_mode=1
if [ ${algo} == "matc" ]; then
	dis_mode=2
fi
exe=optkit
parallel_mode=1
checkiferr=check_iferr.py
#remove this file to clean up the previous result
rm $opt_tmp_result_folder$file 
for f in $(ls ${assign_folder}${file})
do
    #for debug purposes
    #iferr=`python ${python_folder}${checkiferr} ${assign_folder}${file}/$f`
    #if [ "${iferr}" -eq "2" ] 
    #then
    #    continue
    #fi
    
	if [ "$debug" -eq "1" ]
	then #enter debug mode
		echo lldb -- ${folder_now}src/${exe} --dis --dis_mode ${dis_mode} --input_file ${assign_folder}${file}/$f --p_mode $parallel_mode --opt_file $opt_tmp_result_folder$file --seq_len $seqlen 
	else #enter execution mode
		echo processing: ${assign_folder}${file}/$f dismode: ${dis_mode}
		${folder_now}src/${exe} --dis --dis_mode ${dis_mode} --input_file ${assign_folder}${file}/$f --p_mode $parallel_mode --opt_file $opt_tmp_result_folder$file --seq_len $seqlen 
	fi
done
