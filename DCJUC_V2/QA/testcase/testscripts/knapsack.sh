# this is the test case script for knapsack problem
# edited by Zhaoming Yin on 01-02-2015

# definition of some variables
export DATASET_FOLDER=/Users/zhaomingyin/gitlocal/optkit_github/data/knapsack/
export MODE_SEQ=1
export TGT_FOLDER=/Users/zhaomingyin/gitlocal/optkit_github/QA/testcase/tgt_logs/
export TGT_LOG=${TGT_FOLDER}knapsack.log
export RFC_FOLDER=/Users/zhaomingyin/gitlocal/optkit_github/QA/testcase/rfc_logs/
export RFC_LOG=${RFC_FOLDER}knapsack.log

# step one: run different data sets, and generate logs
# in this script every run is serial
echo -n "" > ${TGT_LOG}
for f in $(ls ${DATASET_FOLDER});
do
    ./optkit --knap --input_file ${DATASET_FOLDER}${f} --p_mode ${MODE_SEQ} --num_t 1 >> ${TGT_LOG} 
done

# step two: compare logs and generate diffs
diff -q ${RFC_LOG} ${TGT_LOG} 
