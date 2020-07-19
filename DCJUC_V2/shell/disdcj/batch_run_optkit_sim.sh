HOME=/scratch/zyin/optkit/optkit/
SHELL_HOME=${HOME}shell/disdcj/
EXE_HOME=${HOME}shell/exe/
DATA_HOME=${HOME}data/dist/
IS_DUAL=dual
IS_BALANCED=unbalanced
GRAPH_HOME=${DATA_HOME}graph/1000_0.1_0.1_0.1_${IS_DUAL}_${IS_BALANCED}/
SEQLEN=1000
GAMMA=0.1
PHI=0.1

PYTHON_HOME=${HOME}python/disdcj/
COLLECT_DIS=collect_dis.py
for algo in "exem" "matc"
do
    RES_FILE=${HOME}data/results/distaccy/result_optkit_sim_${algo}
    rm ${RES_FILE}
    for file in $(ls ${GRAPH_HOME})
    do
    	echo processing $file >> ${RES_FILE}
    	echo processing $file 
    	start=`${EXE_HOME}get_time`
    	${SHELL_HOME}run_optkit_dist.sh  \
                $SEQLEN \
                $GAMMA \
                $PHI \
                $file \
                ${IS_DUAL} \
                ${IS_BALANCED} \
                "0" \
                ${algo} >> ${RES_FILE} 
    
    	end=`${EXE_HOME}get_time`
    	echo -n "time " >> ${RES_FILE}
    	bc -l<<<"(1.0*(${end}-${start}))/(2.4*10^9)" >> ${RES_FILE}
    	echo finish processing $file >> ${RES_FILE}
    done
done
