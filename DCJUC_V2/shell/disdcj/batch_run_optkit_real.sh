#home folder
#HOME=/Users/zhaomingyin/gitlocal/optkit/
HOME=/scratch/zyin/optkit/optkit/
#executable folder
SHELL_HOME=${HOME}shell/disdcj/
EXE_HOME=${HOME}shell/exe/
PYTHON_HOME=${HOME}python/disdcj/
#data folder
DATA_HOME=${HOME}data/real/optkit/
SOURCE_DATA_HOME=${DATA_HOME}/source/
GRAPH_HOME=${DATA_HOME}dup_graph/
#executables
O2G=o_2_g.py
DISTEXE=run_optkit_dist.sh
LENEXE=get_seqlen.py
# set if debug mode or not
DEBUG=0
DIST_GRAPH_HOME=${HOME}data/dist/graph/${SEQLEN}_0.1_0_0_dual_balanced

for algo in 'exem' 'matc' 
do
    RES_FILE=${HOME}data/results/distaccy/result_optkit_real_${algo}
    echo ${RES_FILE}
    rm ${RES_FILE}
    for file1 in $(ls ${SOURCE_DATA_HOME})
    do
    	for file2 in $(ls ${SOURCE_DATA_HOME})
    	do
    		if [ $(expr "$file1" \< "$file2") -eq 1 ];
    		then
    			echo processing $file1 $file2 >> ${RES_FILE}
    			echo processing $file1 $file2 
    
    			#get seq len
    			SEQLEN=`python ${PYTHON_HOME}${LENEXE} ${GRAPH_HOME}${file1}_${file2}`
    
                # copy graph files
                if [ ! -d "${DIST_GRAPH_HOME}" ]; then
                    mkdir ${DIST_GRAPH_HOME} 
                fi
                cp ${GRAPH_HOME}${file1}_${file2} ${DIST_GRAPH_HOME} 
    
    			start=`${EXE_HOME}get_time`
    			#real computation
                if [ ${DEBUG} == 1 ]
                then
                    echo ${SHELL_HOME}/${DISTEXE} \
                    ${SEQLEN} \
                    0 \
                    0 \
                    ${file1}_${file2} \
                    "dual" \
                    "balanced" \
                    "0" \
                    "exem" 
                else
                    ${SHELL_HOME}/${DISTEXE} \
                    ${SEQLEN} \
                    0 \
                    0 \
                    ${file1}_${file2} \
                    "dual" \
                    "balanced" \
                    "0" \
                    "${algo}" >> ${RES_FILE} 
    
                fi
    			end=`${EXE_HOME}get_time`
    			echo -n "time " >> ${RES_FILE}
    			bc -l<<<"(1.0*(${end}-${start}))/(2.4*10^9)" >> ${RES_FILE}
    		fi
    	done
    done
done
