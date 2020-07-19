
HOME=/scratch/zyin/optkit/optkit/
GRAPH_HOME=${HOME}data/median/graph/100_0.1_0.0_0.1_balanced/
TMP_HOME=${HOME}data/median/tmp/

for file in $(ls ${GRAPH_HOME})
do
	echo $file
    timeout 1200 ./src/optkit --median --input_file ${GRAPH_HOME}${file} --p_mode 1 --heu_level 2 --term_move 3 --is_opt 1 --dis_mode 1 --tmp_folder ${TMP_HOME} 
done
