#!/bin/sh
set -x

names="index about notes plan lectures exercise_delivery default_project software_packages"
dest=../web
template="--html_template=uio.html"

for name in $names; do
doconce format html $name $template
cp $name.html $dest/
done

name=exam_candidates
doconce format html $name $template --encoding=utf-8
cp $name.html $dest/

name=improvements13
doconce format html $name
cp $name.html $dest/

# These require full math
names="oblig2 oblig3"
names="oblig3"
names=""
for name in $names; do
preprocess -DFORMAT=html ../slides/src/newcommands_keep.p.tex > newcommands_keep.tex
doconce format html $name

preprocess -DFORMAT=pdflatex ../slides/src/newcommands_keep.p.tex > newcommands_keep.tex
doconce format pdflatex $name --device=paper
doconce ptex2tex $name
pdflatex $name
pdflatex $name
cp $name.html $name.pdf $dest/
cp $name.html $dest/
done

#cp background_orig.do.txt background.do.txt
#name=background
#doconce format pdflatex $name
#doconce ptex2tex $name
#pdflatex $name
#doconce format html $name --html_style=bloodish
#cp $name.pdf $name.html $dest/

doconce clean  # must do before next line in order to clean background.*
rm -f background.do.txt




