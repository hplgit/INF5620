#!/bin/sh
# Automate making the plots FE1.png, BE1.png, CN1.png from decay4.py
prog=dc_plot_mpl.py

#python $prog

versions="FE BE CN"
dts="0.4 0.04"
results=""
for v in $versions; do
  files=""
  for dt in $dts; do
    files="$files ${v}_$dt.png"
  done
  result="${v}1.png"
  montage -background white -geometry 100% -tile 2x $files $result
  results="$results $result"

  files=""
  for dt in $dts; do
    files="$files ${v}_$dt.pdf"
  done
  result="${v}1.pdf"
  pdftk $files output tmp.pdf
  pdfnup --nup 2x1 tmp.pdf
  mv tmp-nup.pdf $result
  results="$results $result"
done
echo "mv $results ../fig-decay/"

