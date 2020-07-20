#!/bin/bash

####################
# general parameters
####################
if [ $# -ne 8 ]; then 
    echo ""
    echo "Illegal number of parameters, please follow:"
    echo -n "run_dist_optkit.sh " 
    echo -n "<seqlen> "
    echo -n "<gamma> "
    echo -n "<phi> "
    echo -n "<file> "
    echo -n "<single|dual> "
    echo -n "<balanced|unbalanced> "
    echo -n "<is_bij(1)|not_bij(0)> "
    echo "<exem|matc> "
    echo ""
    exit $?
fi

####################
# general parameters
####################
debug=0                 ##< is debug mode or not
exe=optkit              ##< indicate the executable name
parallel_mode=1         ##< which parallel mode 1:seq 2:buck 3:thread
seqlen=$1               ##< length of the sequence
gamma=$2                ##< indel rate
phi=$3                  ##< duplication rate
theta=0.1               ##< mutation rate
file=$4                 ##< file name
is_single=$5            ##< is evolved in single sequence or not
is_balanced=$6          ##< if two sequence have balanced content
is_bij=$7               ##< is calculate bijection rate or not
                        ##< the general config based on input parameters
algo=$8                 ##< indicate algorithm exem or matc
dis_mode=1              
if [ ${algo} == "matc" ]; then
	dis_mode=2
fi                      ##< which distance mode 1:exem 2:matc
config=${seqlen}_${theta}_${gamma}_${phi}_${is_single}_${is_balanced}
gdb_str="gdb --args"

####################
# folder information
####################
#folder_now=/Users/zhaomingyin/gitlocal/optkit/
folder_now=/Users/zhaomingyin/gitlocal/DCJUC/DCJUC_V2/
                        ##< current folder
graph_folder=${folder_now}data/dist/graph/${config}/
                        ##< folder that stores graph
shell_folder=${folder_now}shell/disdcj/
                        ##< folder that stores shell scripts
python_folder=${folder_now}python/disdcj/
                        ##< folder that stores python scripts
assign_folder=${folder_now}data/dist/assign/${config}/
                        ##< folder to store optimal assignment files
vis_folder=${folder_now}data/dist/vis/${config}/
                        ##< folder to store visualization files
bij_folder=${folder_now}data/dist/bij/
                        ##< bijection folder
dict_folder=${folder_now}data/dist/dict/${config}/
                        ##< folder to keep dictionary files
irrCC_folder=${folder_now}data/dist/irrCC/${config}/
                        ##< to store irregular connected component graphs
opt_tmp_result_folder=${folder_now}data/dist/opt_tmp_result/${config}/
                        ##< keep temporary results

####################
# intermediate file info
####################
vis_file=${vis_folder}${file}
                        ##< visualization file
bij_file=${bij_folder}${config}
                        ##< bijection file

####################
# executable/script information
####################
irrCC=irrCC.py          ##< to calculate irregular connected components
assign_condense=assign_condense.py
                        ##< optimal assignment and condense the graph
collect_dis=collect_dis.py

####################
# make directories
####################
if [ ! -d "${irrCC_folder}" ]; then
        mkdir -p ${irrCC_folder}
fi

if [ ! -d "${opt_tmp_result_folder}" ]; then
        mkdir -p ${opt_tmp_result_folder}
fi

if [ ! -d "${assign_folder}" ]; then
        mkdir -p ${assign_folder}
fi

if [ ! -d "${vis_folder}" ]; then
        mkdir -p ${vis_folder}
fi

if [ ! -d "${bij_folder}" ]; then
        mkdir -p ${bij_folder}
fi

if [ ! -d "${dict_folder}" ]; then
        mkdir -p ${dict_folder}
fi

echo come her

####################
# clean up previous result	
####################
rm $opt_tmp_result_folder$file 

####################
# remove regular CCs	
####################
python  $python_folder$irrCC \
        $graph_folder$file \
        $irrCC_folder$file \
        $opt_tmp_result_folder$file \
        ${vis_file}

echo here

####################
# perform opt assignment and condense
####################
python  $python_folder$assign_condense \
        $irrCC_folder$file \
        $assign_folder$file \
        $opt_tmp_result_folder$file \
        $algo \
        ${vis_file} \
        ${bij_file} \
        ${dict_folder}


echo here

echo ${assign_folder}${file}

####################
# visualization
####################
#./shell/batch_vis_irrCC.sh


####################
# rm correct (unneccessary) files for debugging
####################
#keep_file=/scratch/zyin/optkit/optkit/data/dist/assign/25804_0.1_0_0_dual_balanced/brownrat_zebrafish/0_4.gr 
#for d_f in $(ls ${assign_folder}${file})
#do
#    if [ "${assign_folder}${file}/${d_f}" != "${keep_file}" ] 
#    then
#        rm ${assign_folder}${file}/${d_f}
#    fi
#done

####################
# call C++ executable to really execute the file
####################
for f in $(ls ${assign_folder}${file})
do
    echo "come here"
    # prepare executable strings
    bij_exe_str="${folder_now}src/${exe} --dis --dis_mode ${dis_mode} "
    bij_exe_str+="--input_file ${assign_folder}${file}/$f "
    bij_exe_str+="--p_mode $parallel_mode "
    bij_exe_str+="--opt_file $opt_tmp_result_folder$file "
    bij_exe_str+="--seq_len $seqlen --cal_bij "
    bij_exe_str+="--dict_file ${dict_folder}$f "
    bij_exe_str+="--bij_file ${bij_file}"
    
    nobij_exe_str="${folder_now}src/${exe} --dis --dis_mode ${dis_mode} "
    nobij_exe_str+="--input_file ${assign_folder}${file}/$f "
    nobij_exe_str+="--p_mode $parallel_mode "
    nobij_exe_str+="--opt_file $opt_tmp_result_folder$file "
    nobij_exe_str+="--seq_len $seqlen"
    # execution 
	if [ "$debug" -eq "1" ]
	then 
        #enter debug mode
        if [ "$is_bij" -eq "1" ]
        then
            echo ""
		    echo ${gdb_str} ${bij_exe_str}
            echo ""
        else
            echo ""
		    echo ${gdb_str} ${nobij_exe_str}
            echo ""
        fi
	else 
        #enter execution mode
		echo processing: ${assign_folder}${file}/$f dismode: ${dis_mode}
        if [ "$is_bij" -eq "1" ]
        then
            ${bij_exe_str}
        else
            ${nobij_exe_str}
        fi
	fi
    # processing distance results here
    echo -n "distance: "
    python ${python_folder}${collect_dis} $opt_tmp_result_folder$file ${seqlen}
done
