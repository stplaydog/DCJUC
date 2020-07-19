# defines folders to store data and source code
HOME=/Users/zhaomingyin/gitlocal/optkit/
SHELL_HOME=${HOME}shell/disdcj/
DATA_HOME=${HOME}data/real/gredo/dup_pair/
SOURCE_DATA_HOME=${HOME}data/real/optkit/source/
EXE_HOME=${HOME}shell/exe/
BUILD_EXE=build.pl
TIME_LIMIT_SEC=3600

# defines executables

for file1 in $(ls ${SOURCE_DATA_HOME})
do
	for file2 in $(ls ${SOURCE_DATA_HOME})
	do
		if [ $(expr "$file1" \< "$file2") -eq 1 ];
		then
			echo processing $file1 $file2
			start=`${EXE_HOME}get_time`
			${EXE_HOME}gredo ${DATA_HOME}${file1}${file2}_${file1} ${DATA_HOME}${file1}${file2}_${file2} ${TIME_LIMIT_SEC} 
			end=`${EXE_HOME}get_time`
			echo -n "time "
			bc -l<<<"(1.0*(${end}-${start}))/(2.4*10^9)"
		fi
	done
done
