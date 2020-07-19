HOME=/Users/zhaomingyin/gitlocal/optkit/
EXE_HOME=${HOME}shell/exe/
DATA_HOME=${HOME}data/dist/
SEQ_HOME=${DATA_HOME}seq/1000_0.1_0.0_0.1_dual_balanced/
GREDO_HOME=${DATA_HOME}gredo_seq/1000_0.1_0.0_0.1_dual_balanced/

for file in $(ls ${SEQ_HOME})
do
	echo processing $file
	start=`${EXE_HOME}get_time`
	${EXE_HOME}gredo ${GREDO_HOME}${file}_1 ${GREDO_HOME}${file}_2 3600 
	end=`${EXE_HOME}get_time`
	echo -n "time "
	bc -l<<<"(1.0*(${end}-${start}))/(2.4*10^9)"
	echo finish processing $file
done
