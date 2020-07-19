mode_seq=1
mode_par=2
mode_tpar=3

file=data/dist/seqnew/1.0_0.05_0.1_1

for i in {1..10}
do
	echo test $i num_t 1
	timeout 20 ./optec --dcj --input_file $file --p_mode $mode_seq --num_t 1
done


echo ----------------------------------------------------------------
for i in 2 4 8 16
do 
	for j in {1..10}
	do
		echo test $j num_t $i
		timeout 20 ./optec --dcj --input_file $file --p_mode $mode_par --num_t $i
	done
done


echo ----------------------------------------------------------------
for i in 2 4 8 16
do 
	for j in {1..10}
	do
		echo test $j num_t $i
		timeout 20 ./optec --dcj --input_file $file --p_mode $mode_tpar --num_t $i
	done
done
