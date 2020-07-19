folder_now=/Users/zhaomingyin/Dropbox/research/code/optec-code/
pair_folder=${folder_now}data/real_pair/
graph_folder=${folder_now}data/pair_graph/
shell_folder=${folder_now}shell/
irrCC=irrCC.py
irrCC_folder=${folder_now}data/irrCC/
vis_folder=${folder_now}data/vis/
opt_tmp_result_folder=${folder_now}data/opt_tmp_result/
assign_condense=assign_condense.py
o_2_g=o_2_g.py
assign_folder=${folder_now}data/assign/
algo=exem

for file in $(ls $pair_folder)
do
	echo $file
	#turn order files into graphs first
	python $shell_folder$o_2_g $pair_folder$file $graph_folder$file 
	#remove regular CCs	
	python $shell_folder$irrCC $graph_folder$file $irrCC_folder$file $opt_tmp_result_folder$file ${vis_folder}$file
	#perform opt assignment and condense
	python  $shell_folder$assign_condense $irrCC_folder$file $assign_folder$file $opt_tmp_result_folder$file $algo
done
