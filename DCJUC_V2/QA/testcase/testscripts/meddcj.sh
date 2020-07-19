#
# created 03/14/2015 by Zhaoming Yin
# this test case is for the test of median computation using LK method
# distance metric are DCJ-Indel-Exemplar and DCJ-Indel-Matching
#

HOME=/Users/zhaomingyin/gitlocal/optkit_github/
GRAPH_HOME=${HOME}data/meddcj/graph/
TMP_HOME=${HOME}data/median/tmp/
EXE=${HOME}src/optkit

for file in $(ls ${GRAPH_HOME})
do
	echo $file
    ${EXE} --median --input_file ${GRAPH_HOME}${file} --p_mode 1 --heu_level 2 --term_move 3 --is_opt 1 --dis_mode 1 --tmp_folder ${TMP_HOME} 
done

