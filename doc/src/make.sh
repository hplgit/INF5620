#!/bin/sh
names="about index index notes plan lectures exercise_delivery oblig1 oblig2 default_project"
dest=../web
template="--html_template=uio.html"

for name in $names; do
doconce format html $name $template
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




