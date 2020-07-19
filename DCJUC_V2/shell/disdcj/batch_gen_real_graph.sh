# assign folders
HOME=/Users/zhaomingyin/gitlocal/optkit/
OPTKIT_DATA_HOME=${HOME}data/real/optkit/
GREDO_DATA_HOME=${HOME}data/real/gredo/
SOURCE_DATA_HOME=${OPTKIT_DATA_HOME}source/
PY_HOME=${HOME}python/disdcj/

# specify executables
PY_EXE=gen_real_order.py
O2G=o_2_g.py
GENDATA=gen_real_dup_order.py

# specify destination folder 
MODE=balanced
if [ ${MODE} == 'balanced' ]
then
    GRAPH_HOME=${OPTKIT_DATA_HOME}dup_graph/
    OPTKIT_DEST_HOME=${OPTKIT_DATA_HOME}dup_pair/
    GREDO_DEST_HOME=${GREDO_DATA_HOME}dup_pair/
else
    GRAPH_HOME=${OPTKIT_DATA_HOME}indel_graph/
    OPTKIT_DEST_HOME=${OPTKIT_DATA_HOME}indel_pair/
    GREDO_DEST_HOME=${iGREDO_DATA_HOME}indel_pair/
fi

for file1 in $(ls $SOURCE_DATA_HOME)
do
    for file2 in $(ls $SOURCE_DATA_HOME)
    do
        if [ "$file1" \< "$file2" ]
        then
            echo $file1 $file2
            # generate gredo data pair and optkit sequence
            python ${PY_HOME}${PY_EXE} ${SOURCE_DATA_HOME}$file1 \
                   ${SOURCE_DATA_HOME}$file2 ${MODE} \
                   ${OPTKIT_DEST_HOME}${file1}_${file2} \
                   ${GREDO_DEST_HOME}${file1}${file2}_${file1} \
                   ${GREDO_DEST_HOME}${file1}${file2}_${file2}
            # generate optkit graph
            python ${PY_HOME}${O2G} ${OPTKIT_DEST_HOME}${file1}_${file2} \
                   ${GRAPH_HOME}${file1}_${file2} 
        fi
    done
done
