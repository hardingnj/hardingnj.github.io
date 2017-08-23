#! /bin/bash


NB=$1
filename=$(basename "$NB")
extension="${filename##*.}"
filename="${filename%.*}"
echo $filename

jupyter nbconvert --to markdown $NB

rm -r ../images/${filename}_files
mv ${filename}_files ../images 
sed -i 's/\[png\](/[png](\/images\//' ${filename}.md
