#!/bin/sh -x
dir=reports
sh clean.sh
rm -rf $dir
mkdir $dir

failures=""


report=tmp_report
dt="1.25 0.75 0.5 0.1"
python decay_exper1_do.py $dt
if [ $? -ne 0 ]; then failures="$failures:decay_exer1_do.py"; fi
cp $report.do.txt $dir/report.do.txt
cp *.png $dir

cd $dir
doconce format html report
if [ $? -ne 0 ]; then failures="$failures:doconce-html"; fi
mv report.html report_do.html

doconce sphinx_dir theme=pyramid report
cp *.png sphinx-rootdir
python automake_sphinx.py
cd sphinx-rootdir
doconce replace solarized '' make_themes.sh # have it, but it doesn't work...
sh make_themes.sh
if [ $? -ne 0 ]; then failures="$failures:make_themes.sh"; fi
mv -f sphinx-* ../
cd ..
mv sphinx-rootdir rootdir

doconce format pdflatex report
doconce ptex2tex report envir=minted
# Do some polishing of report.tex for display of latex source
 doconce subst -m  -s '.*%--+ begin preamble -+$' '' report.tex
 doconce subst '.*%--+ end preamble -+' '' report.tex
#doconce replace "ptex2tex (extended LaTeX)" "LaTeX" report.tex
doconce subst -s "demonstrated\..+\\end\{abstract\}" "demonstrated.\n\end{abstract}" report.tex
doconce subst -s "% Purpose: section with computer.*\\end{itemize}" "\end{itemize}" report.tex
doconce subst "% Purpose: .*" "" report.tex
doconce replace '\noindent' '' report.tex
pdflatex -shell-escape report
pdflatex -shell-escape report

doconce format pandoc report
# Remove title, author, etc. (does not work well)
doconce subst '% .*' '' report.mkd
pandoc -s --mathjax -f markdown -t report.mkd.html -o tmp.html report.mkd
doconce subst '</style>' '</style>
<script type="text/x-mathjax-config">
MathJax.Hub.Config({
  TeX: {
     equationNumbers: {  autoNumber: "AMS"  },
     extensions: ["AMSmath.js", "AMSsymbols.js", "autobold.js"]
  }
});
</script>' report.mkd.html

cd ..
python decay_exper1_html.py $dt
if [ $? -ne 0 ]; then failures="$failures:decay_exper1_html.py"; fi
cp $report.html $dir/report_html.html
python decay_exper1_html_mathjax.py $dt
if [ $? -ne 0 ]; then failures="$failures:decay_exper1_html_mathjax.py"; fi
cp $report.html $dir/report_html_mathjax.html

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

pyg="pygmentize -f html -O full,style=emacs"
for file in *.html; do
  $pyg -o $file.html -l html $file
  doconce subst 'body\s+\.err\s+\{ .+' 'body  .err { border: 0; } /* drop error */' $file.html
done
rm -f index.html.html  report.mkd.html.html # not of interest
$pyg -o report_sphinx.rst.html -l rst rootdir/report.rst
$pyg -o report.p.tex.html -l latex report.p.tex
$pyg -o report.tex.html -l latex report.tex
$pyg -o report.mkd.html -l latex report.mkd
$pyg -o report.do.txt.html -l text report.do.txt

rm -f *.aux *.dvi *.log *.idx *.out *.toc tmp* *~ automake* *.tex *.rst *.mkd
mkdir _static
mv -f *.png *.html *.pdf sphinx-* _static
mv -f _static/index.html .

# Make project tree
proj=project_mathjax
rm -rf $proj
mkdir $proj
mkdir $proj/src
cp ../decay_exper1_html_mathjax.py $proj/src
mkdir $proj/doc
cp _static/report_html_mathjax.html _static/*.png $proj/doc/report.html
cat > $proj/doc/run.sh <<EOF
#!/bin/sh
# Run experiment documented in report.html

python decay_exper1_html_mathjax.py 1.25 0.75 0.5 0.1
EOF

proj=project_doconce
rm -rf $proj
mkdir $proj
mkdir $proj/src
cp ../decay_exper1_*.py $proj/src
mkdir $proj/doc
cp -r _static/report_html_mathjax.html _static/report_do.html _static/report.pdf _static/sphinx-* _static/*.png $proj/doc/
cat > $proj/doc/run.sh <<EOF
#!/bin/sh
# Run experiment documented in reports

python decay_exper1_do.py 1.25 0.75 0.5 0.1

# Make reports

# HTML
file=tmp_report
doconce format html $file
mv $file.html report_do.html

# LaTeX
doconce format pdflatex $file
doconce ptex2tex $file envir=minted
pdflatex -shell-escape $file
pdflatex -shell-escape $file
mv $file.pdf report.pdf

# Sphinx
doconce sphinx_dir theme=pyramid report
cp *.png sphinx-rootdir
python automake_sphinx.py

EOF

cat > tmp.do.txt <<EOF
TITLE: The scientific report in different formats

 * "HTML as written by `decay_exper1_html_mathjax.py`": "report_html_mathjax.html"
 * "HTML": "report_do.html" as generated from Doconce by `decay_exper1_do.py`
 * "PDF": "reportpdf" as generated via LaTeX from Doconce by `decay_exper1_do.py`
 * "Sphinx": "sphinx-default/index.html" (default layout)

Here are numerous other Sphinx themes:

 * "agni": "sphinx-agni/report.html"
 * "agogo": "sphinx-agogo/report.html"
 * "basic": "sphinx-basic/report.html"
 * "cbc": "sphinx-cbc/report.html"
 * "classy": "sphinx-classy/report.html"
 * "cloud": "sphinx-cloud/report.html"
 * "default": "sphinx-default/report.html"
 * "epub": "sphinx-epub/report.html"
 * "fenics": "sphinx-fenics/report.html"
 * "fenics_minimal": "sphinx-fenics_minimal/report.html"
 * "flask": "sphinx-flask/report.html"
 * "haiku": "sphinx-haiku/report.html"
 * "jal": "sphinx-jal/report.html"
 * "nature": "sphinx-nature/report.html"
 * "pylons": "sphinx-pylons/report.html"
 * "pyramid": "sphinx-pyramid/report.html"
 * "redcloud": "sphinx-redcloud/report.html"
 * "scrolls": "sphinx-scrolls/report.html"
 * "slim-agogo": "sphinx-slim-agogo/report.html"
 * "solarized": "sphinx-solarized/report.html"
 * "sphinxdoc": "sphinx-sphinxdoc/report.html"
 * "traditional": "sphinx-traditional/report.html"
 * "vlinux-theme: "sphinx-vlinux-theme/report.html"
EOF
doconce format html tmp
mv tmp.html $proj/doc/index.html
#rm tmp*

cd ..
sh clean.sh

# Archive
rm -rf archived-reports
cp -r reports archived-reports
rm -rf archived-reports/rootdir  # not to be archived
#cp -r archived-reports/* ~/vc/INF5620/doc/writing_reports/
rm -rf archived-projects
cp -r reports/project_mathjax archived-projects
cp -r reports/project_doconce archived-projects

echo "failures: $failures"