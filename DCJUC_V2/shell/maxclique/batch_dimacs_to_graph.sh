home=/Users/zhaomingyin/gitlocal/optkit/
data_home=${home}data/MC/Additional_cliques/
out_home=${home}data/MC/DIMACS/
shell_home=${home}shell/

mkdir ${out_home}
for file in $(ls $data_home)
do
	python ${shell_home}/dimacs_to_graph.py ${data_home}$file ${out_home}${file}
done
