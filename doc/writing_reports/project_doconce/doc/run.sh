#!/bin/sh
# Run experiment documented in reports

python decay_exper1_do.py 1.25 0.75 0.5 0.1

# ----- Make reports -----

# Make publish database for bibliography (from BibTeX file refs.bib)
publish import refs

# HTML
file=tmp_report
doconce format html report_wordpress.html
mv -f report_wordpress.html.html report_do.html

# LaTeX
doconce format pdflatex report_wordpress.html
doconce ptex2tex report_wordpress.html envir=minted
pdflatex -shell-escape report_wordpress.html
pdflatex -shell-escape report_wordpress.html
mv -f report_wordpress.html.pdf report.pdf

# Sphinx
doconce sphinx_dir theme=pyramid report
cp *.png sphinx-rootdir
python automake_sphinx.py

