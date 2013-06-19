#!/bin/sh -x
# Make automatically generated Doconce report from a program in
# a wide variety of formats (for demonstrating the layout of
# the various formats)

report=tmp_report
dt="1.25 0.75 0.5 0.1"
dir=reports
failures=""

# Drop cleaning and running simulations if command line arguments
# and figure files exist
if [ $# -gt 0 -a -f BE.png -a -d $dir ]; then
  echo 'No new simulations - just report generation'
else
sh clean.sh
rm -rf $dir
mkdir $dir
python decay_exper1_do.py $dt
fi

if [ $? -ne 0 ]; then failures="$failures:decay_exer1_do.py"; exit 1; fi
cp $report.do.txt $dir/report.do.txt
cp .publish_references.pub publish_config.py BE.* FE.* CN.* error.* $dir

cd $dir

# Solarized doconce html
doconce format html report --html_style=solarized --html_admon=apricot
if [ $? -ne 0 ]; then failures="$failures:doconce-html-solarized"; exit 1; fi
mv -f report.html report_solarized.html

# blogger.com (Google) blog
cp report.do.txt report2.do.txt
# Remove title, author and date
doconce subst -m '^TITLE:.+$' '' report2.do.txt
doconce subst -m '^AUTHOR:.+$' '' report2.do.txt
doconce subst -m '^DATE:.+$' '' report2.do.txt
# Figures must have URLs to where they are stored on the web
for figname in BE FE CN error; do
  doconce replace "[$figname," "[https://raw.github.com/hplgit/hplgit.github.com/master/teamods/writing_reports/_static/$figname.png," report2.do.txt
done
doconce format html report2
mv -f report2.html report_blogger.html
# Paste report_blogger.html text into a new blog on your Google account

# Wordpress blog
doconce format html report2 --wordpress
mv -f report2.html report_wordpress.html

# Plain doconce html blueish
doconce format html report --html_style=blueish --html_admon=yellow
if [ $? -ne 0 ]; then failures="$failures:doconce-html"; fi
mv -f report.html report_blueish.html

# Plain doconce html bloodish
doconce format html report --html_style=bloodish
if [ $? -ne 0 ]; then failures="$failures:doconce-html"; fi
mv -f report.html report_bloodish.html

# Plain doconce html bloodish with handwriting
# Clicker+Script Stalemate Architects+Daughter
doconce format html report --html_style=bloodish --html_body_font=Architects+Daughter --html_heading_font=Architects+Daughter
if [ $? -ne 0 ]; then failures="$failures:doconce-html"; fi
mv -f report.html report_bloodish_Architects_Daughter.html

doconce format html report --html_style=bloodish --html_body_font=Clicker+Script --html_heading_font=Clicker+Script
if [ $? -ne 0 ]; then failures="$failures:doconce-html"; fi
mv -f report.html report_bloodish_Clicker_Script.html

doconce format html report --html_style=bloodish --html_body_font=Stalemate --html_heading_font=Stalemate
if [ $? -ne 0 ]; then failures="$failures:doconce-html"; fi
mv -f report.html report_bloodish_Stalemate.html

# Vagrant doconce html
cp -r ~/vc/doconce/bundled/html_styles/style_vagrant .
cp -r style_vagrant/css css
# Customize the template
doconce replace LogoWord 'DiffEq' style_vagrant/template_vagrant.html
doconce replace withSubWord 'Project' style_vagrant/template_vagrant.html
doconce replace '<a href="">GO TO 1</a>' '<a href="http://wikipedia.org">Wikipedia</a>' style_vagrant/template_vagrant.html
doconce replace '<a href="">GO TO 2</a>' '<a href="http://wolframalpha.com">WolframAlpha</a>' style_vagrant/template_vagrant.html
doconce replace 'Here goes a footer, if desired, maybe with a Copyright &copy;' 'H. P. Langtangen &copy; 2013' style_vagrant/template_vagrant.html

# TOC, TITLE and AUTHOR are not appropriate with HTML templates
cp report.do.txt report2.do.txt
doconce replace 'TITLE: '  '#TITLE: ' report2.do.txt
doconce replace 'AUTHOR: '  '#AUTHOR: ' report2.do.txt
doconce replace 'TOC: '  '#TOC: ' report2.do.txt
doconce format html report2 --html_style=vagrant --html_template=style_vagrant/template_vagrant.html
if [ $? -ne 0 ]; then failures="$failures:doconce-html-vagrant"; fi

mv -f report2.html report_vagrant.html

# GitHub minimal theme
cp -r ~/vc/doconce/bundled/html_styles/style_github_minimal .
mv -f style_github_minimal/css/* css/
mv -f style_github_minimal/js .
doconce replace 'Main Permanent Header' 'Project Report' style_github_minimal/template_github_minimal.html
doconce replace 'Permanent SubHeader' 'A Differential Equation Problem' style_github_minimal/template_github_minimal.html

doconce format html report2 --html_template=style_github_minimal/template_github_minimal.html
if [ $? -ne 0 ]; then failures="$failures:doconce-html-github-minimal"; fi

mv -f report2.html report_github_minimal.html
rm -f report2*

# IPython notebook
doconce format ipynb report
if [ $? -ne 0 ]; then failures="$failures:doconce-ipynb"; fi
cp report.ipynb ~/vc/hplgit.github.com/store/  # store so it's reachable on web

# MediaWiki
doconce format mwiki report
if [ $? -ne 0 ]; then failures="$failures:doconce-mwiki"; fi
cp report.mwiki ~/vc/hplgit.github.com.wiki/Experiments-with-Schemes-for-Exponential-Decay.mediawiki

# Sphinx
doconce sphinx_dir theme=pyramid report
python automake_sphinx.py
cd sphinx-rootdir
doconce replace solarized '' make_themes.sh # have it, but it doesn't work...
sh make_themes.sh
if [ $? -ne 0 ]; then failures="$failures:make_themes.sh"; fi
mv -f sphinx-* ../
cd ..
mv -f sphinx-rootdir rootdir

# PDF for printing
doconce format pdflatex report --device=paper
doconce ptex2tex report envir=minted -DADMON=box -DPALATINO
rm -f *.aux
pdflatex -shell-escape report
bibtex report
pdflatex -shell-escape report
cp report.pdf report_4printing.pdf

# PDF for phone
doconce format pdflatex report
doconce ptex2tex report envir=minted -DA6PAPER -DPALATINO
rm -f *.aux
pdflatex -shell-escape report
bibtex report
pdflatex -shell-escape report
cp report.pdf report_4phone.pdf

# PDF for screen viewing
doconce format pdflatex report
doconce ptex2tex report envir=minted  -DPALATINO -ADMON=colors
# Do some polishing of report.tex for display of the latex source
# to the world
doconce subst -m '^%--+ begin preamble -+$' '' report.tex
doconce subst -m '^%--+ end preamble -+$' '' report.tex
doconce subst -m '^% --+ main content -+$' '' report.tex
doconce subst -m '^% --+ end of main content -+$' '' report.tex
#doconce replace "ptex2tex (extended LaTeX)" "LaTeX" report.tex
doconce subst -s "demonstrated\..+\\end\{abstract\}" "demonstrated.\n\end{abstract}" report.tex
doconce replace '\noindent' '' report.tex
rm -f *.aux
pdflatex -shell-escape report
bibtex report
pdflatex -shell-escape report

# Markdown
doconce format pandoc report
# Remove title, author, etc. (does not work well)
doconce subst '% .*' '' report.md
doconce md2html report
cp report.html report_md.html

# Native HTML and LaTeX formats
cd ..
python decay_exper1_html.py $dt
if [ $? -ne 0 ]; then failures="$failures:decay_exper1_html.py"; exit 1; fi
cp $report.html $dir/report_html.html

python decay_exper1_mathjax.py $dt
if [ $? -ne 0 ]; then failures="$failures:decay_exper1_mathjax.py"; exit 1; fi
cp $report.html $dir/report_mathjax.html

python decay_exper1_latex.py $dt
if [ $? -ne 0 ]; then failures="$failures:decay_exper1_latex.py"; exit 1; fi
pdflatex $report
bibtex $report
pdflatex $report
cp $report.pdf $dir/report_plain_latex.pdf
cp $report.tex $dir/report_plain_latex.tex

# Make top index.html file
cp decay_report_demo.do.txt $dir/tmp.do.txt
cd $dir
themes=`/bin/ls -d sphinx-*`
for theme in $themes; do
    doconce replace XXXXX "\"$theme\": \"_static/$theme/index.html\", XXXXX" tmp.do.txt
done
doconce replace ", XXXXX" "" tmp.do.txt
doconce format html tmp
if [ $? -ne 0 ]; then failures="$failures:doconce-reports/tmp.do.txt"; fi
mv -f tmp.html index.html

# Compile with shell instructions
doconce replace 'TITLE: Examples of scientific reports in different formats' 'TITLE: Examples of scientific reports in different formats and how they are made' tmp.do.txt
doconce format html tmp -DCODE
if [ $? -ne 0 ]; then failures="$failures:doconce-reports/tmp.do.txt"; fi
mv -f tmp.html doconce_commands.html

pyg="pygmentize -f html -O full,style=emacs"
for file in *.html; do
  $pyg -o $file.html -l html $file
  doconce subst 'body\s+\.err\s+\{ .+' 'body  .err { border: 0; } /* drop error */' $file.html
done
rm -f index.html.html  report.md.html.html # not of interest
$pyg -o report_sphinx.rst.html -l rst rootdir/report.rst
$pyg -o report.p.tex.html -l latex report.p.tex
$pyg -o report.tex.html -l latex report.tex
$pyg -o report.md.html -l latex report.md
$pyg -o report.ipynb.html -l text report.ipynb
$pyg -o report.mwiki.html -l text report.mwiki
$pyg -o report_latex.html -l latex report_plain_latex.tex

doconce pygmentize report.do.txt perldoc

rm -f *.aux *.dvi *.log *.idx *.out *.toc tmp* *~ automake* *.tex *.rst *.md
mkdir _static
mv -f *.png *.html *.pdf sphinx-* css js *.ipynb *.mwiki   _static
mv -f _static/index.html _static/doconce_commands.html .  # don't copy the index file

# Make project tree
proj=project_mathjax
rm -rf $proj
mkdir $proj
mkdir $proj/src
cp ../decay_exper1_mathjax.py $proj/src
mkdir $proj/doc
cp _static/report_mathjax.html $proj/doc/report.html
cp _static/*.png $proj/doc/
cat > $proj/doc/run.sh <<EOF
#!/bin/sh
# Run experiment documented in report.html

python decay_exper1_mathjax.py 1.25 0.75 0.5 0.1
EOF

proj=project_doconce
rm -rf $proj
mkdir $proj
mkdir $proj/src
cp ../decay_exper1_*.py $proj/src
mkdir $proj/doc
cp -r _static/report_mathjax.html _static/report_blueish.html _static/report*.pdf _static/sphinx-* _static/*.png $proj/doc/
cat > $proj/doc/run.sh <<EOF
#!/bin/sh
# Run experiment documented in reports

python decay_exper1_do.py 1.25 0.75 0.5 0.1

# ----- Make reports -----

# Make publish database for bibliography (from BibTeX file refs.bib)
publish import refs

# HTML
file=tmp_report
doconce format html $file
mv -f $file.html report_do.html

# LaTeX
doconce format pdflatex $file
doconce ptex2tex $file envir=minted
pdflatex -shell-escape $file
pdflatex -shell-escape $file
mv -f $file.pdf report.pdf

# Sphinx
doconce sphinx_dir theme=pyramid report
cp *.png sphinx-rootdir
python automake_sphinx.py

EOF


cd ..

# Archive
rm -rf archived-reports
cp -r $dir archived-reports
rm -rf archived-reports/rootdir  # not to be archived
#cp -r archived-reports/* ~/vc/INF5620/doc/writing_reports/

#sh clean.sh
echo "failures: $failures"
exit

# The script below is not used - it is decay_report_demo.do.txt that
# becomes the index.html file.
cat > tmp.do.txt <<EOF
TITLE: The scientific report in different formats

 * "HTML as written by `decay_exper1__mathjax.py`": "_static/report_mathjax.html"
 * "HTML": "_static/report_do.html" as generated from Doconce by `decay_exper1_do.py`
 * "PDF": "_static/report_plain.pdf" as generated by `decay_exper1_latex.py` ("LaTeX file": "_static/report_latex.html")
 * "PDF": "_static/report.pdf" for *electronic view* (as generated via LaTeX from Doconce by `decay_exper1_do.py`)
 * "PDF": "_static/report_4printing.pdf" for *printing* (as generated via LaTeX from Doconce by `decay_exper1_do.py`)
 * "PDF": "_static/report_4phone.pdf" for *viewing on phones* (as generated via LaTeX from Doconce by `decay_exper1_do.py`)
 * "Sphinx": "_static/sphinx-default/index.html" (default layout)

Here are numerous other Sphinx themes:

 * "agni": "_static/sphinx-agni/report.html"
 * "agogo": "_static/sphinx-agogo/report.html"
 * "basic": "_static/sphinx-basic/report.html"
 * "basicstrap": "_static/sphinx-basicstrap/report.html"
 * "bootstrap": "_static/sphinx-bootstrap/report.html"
 * "cbc": "_static/sphinx-cbc/report.html"
 * "classy": "_static/sphinx-classy/report.html"
 * "cloud": "_static/sphinx-cloud/report.html"
 * "default": "_static/sphinx-default/report.html"
 * "epub": "_static/sphinx-epub/report.html"
 * "fenics": "_static/sphinx-fenics/report.html"
 * "fenics_minimal": "_static/sphinx-fenics_minimal/report.html"
 * "flask": "_static/sphinx-flask/report.html"
 * "haiku": "_static/sphinx-haiku/report.html"
 * "jal": "_static/sphinx-jal/report.html"
 * "nature": "_static/sphinx-nature/report.html"
 * "pylons": "_static/sphinx-pylons/report.html"
 * "pyramid": "_static/sphinx-pyramid/report.html"
 * "redcloud": "_static/sphinx-redcloud/report.html"
 * "scrolls": "_static/sphinx-scrolls/report.html"
 * "slim-agogo": "_static/sphinx-slim-agogo/report.html"
 * "solarized": "_static/sphinx-solarized/report.html"
 * "sphinxdoc": "_static/sphinx-sphinxdoc/report.html"
 * "traditional": "_static/sphinx-traditional/report.html"
 * "vlinux-theme": "_static/sphinx-vlinux-theme/report.html"
EOF
doconce format html tmp
mv -f tmp.html $proj/doc/index.html
#rm tmp*
