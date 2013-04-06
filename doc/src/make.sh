#!/bin/sh
#names="about index notes plan handwritings background"
names="about index notes plan handwritings exercise_delivery oblig1 oblig2 default_project"
dest=../web
template="--html-template=uio.html"

#cp background_orig.do.txt background.do.txt
#name=background
#doconce format pdflatex $name
#doconce ptex2tex $name
#pdflatex $name
#cp $name.pdf $dest/

for name in $names; do
doconce format html $name $template
cp $name.html $dest/
done

doconce clean  # must do before next line in order to clean background.*
rm -f background.do.txt




