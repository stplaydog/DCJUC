base_folder=/Users/zhaomingyin/Dropbox/research/code/optec-code/data/dist/
folder=${base_folder}/vis/$1/
png_folder=${base_folder}/$1/

if [ ! -d "${png_folder}" ]; then
	mkdir $png_folder
fi

for file in $(ls $folder)
do
	dot -Tpng $folder$file > $png_folder$file.png 
done
