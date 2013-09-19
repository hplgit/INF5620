#!/bin/sh
# Automate the run of four dt cases for the decay_odespy1.py script.

time_steps="1.05 0.75 0.5 0.25"
backend=matplotlib
#backend=gnuplot
plotfiles=""

for dt in $time_steps; do
   python decay_odespy1.py $dt --SCITOOLS_easyviz_backend $backend
   plotfiles="$plotfiles odespy1_dt_$dt.png"
done

montage -background white -geometry 100% -tile 2x $plotfiles decay_odespy1_png.png

