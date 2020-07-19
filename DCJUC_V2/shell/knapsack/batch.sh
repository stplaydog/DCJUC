for t in 2 4 8 16 
do
	for i in {0..99}
	do
		echo "num t" $t	 "file" $i
	    ./knap data/p07 $t
	done
done
