HOME=/Users/zhaomingyin/gitlocal/optkit/
SHELL_HOME=${HOME}shell/disdcj/
PYTHON_HOME=${HOME}python/disdcj/
DATA_HOME=${HOME}data/dist/
SEQ_HOME=${DATA_HOME}seq/1000_0.1_0.0_0.1_dual_balanced/
GREDO_HOME=${DATA_HOME}gredo_seq/1000_0.1_0.0_0.1_dual_balanced/

if [ ! -d "${GREDO_HOME}" ]; then
        mkdir ${GREDO_HOME}
fi

for file in $(ls ${SEQ_HOME})
do
	echo process $file
	python ${PYTHON_HOME}gr_to_gredo.py ${SEQ_HOME}$file ${GREDO_HOME}$file
done
