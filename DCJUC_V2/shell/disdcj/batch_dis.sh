datahome=/home/leonxin/optec-code/data/dist/seq/
shellhome=/home/leonxin/optec-code/shell/
newdatahome=/home/leonxin/optec-code/data/dist/seqnew/
exe=/home/leonxin/optec-code/optec
num_t=1

for file in $(ls $datahome)
do
	echo $file
	python ${shellhome}to_new_format.py ${datahome}$file ${newdatahome}$file
	$exe ${newdatahome}$file $num_t 
done
