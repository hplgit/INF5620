#!/bin/sh
# Command-line argument 1: lecture_name.do.txt
dofile=$1
if [ ! -f $dofile ]; then
  echo "No such file: $dofile"
  exit 1
fi

filename=`echo $dofile | sed 's/\.do\.txt//'`

# LaTeX PDF for printing
# (Smart to make this first to detect latex errors - HTML/MathJax
# gives far less errors and warnings)
rm -f *.aux
preprocess -DFORMAT=pdflatex ../newcommands_keep.p.tex > newcommands_keep.tex
doconce format pdflatex $filename --device=paper -DWITH_TOC
if [ $? -ne 0 ]; then echo "doconce could not compile document $filename.do.txt - abort"; exit; fi

ptex2tex $filename
pdflatex $filename
makeindex $filename
pdflatex $filename
cp ${filename}.pdf ${filename}-4print.pdf

# HTML with solarized style and one big file
preprocess -DFORMAT=html ../newcommands_keep.p.tex > newcommands_keep.tex
doconce format html $filename --html_style=solarized --html_output=${filename}-solarized --pygments_html_style=perldoc --pygments_html_linenos  -DWITH_TOC
if [ $? -ne 0 ]; then echo "doconce could not compile document $filename.do.txt - abort"; exit; fi
doconce replace "<li>" "<p><li>" ${filename}-solarized.html

# reveal.js HTML5 slides
doconce format html $filename --pygments_html_style=native --keep_pygments_html_bg --html_output=${filename}-reveal
doconce slides_html ${filename}-reveal reveal --html_slide_theme=darkgray

# deck.js HTML5 slides
doconce format html $filename --pygments_html_style=perldoc --keep_pygments_html_bg --html_output=${filename}-deck
doconce slides_html ${filename}-deck deck --html_slide_theme=sandstone.default

# Plain HTML with everything in one file
doconce format html $filename --html_style=bloodish --html_output=${filename}-1 -DWITH_TOC
doconce replace "<li>" "<p><li>" ${filename}-1.html

# Plain HTML with one page per slide
doconce format html $filename --html_style=bloodish -DWITH_TOC
doconce replace "<li>" "<p><li>" ${filename}.html
doconce split_html ${filename}.html

# LaTeX Beamer
rm -f *.aux
preprocess -DFORMAT=pdflatex ../newcommands_keep.p.tex > newcommands_keep.tex
if [ $? -ne 0 ]; then echo "doconce could not compile document $filename.do.txt - abort"; exit; fi
doconce format pdflatex $filename
doconce ptex2tex $filename -DLATEX_HEADING=beamer envir=minted
doconce slides_beamer $filename --beamer_slide_theme=red_shadow
pdflatex -shell-escape $filename
if [ $? -ne 0 ]; then echo "pdflatex could not compile document $filename.do.txt- abort"; exit; fi
pdflatex -shell-escape $filename
cp ${filename}.pdf ${filename}-beamer.pdf
