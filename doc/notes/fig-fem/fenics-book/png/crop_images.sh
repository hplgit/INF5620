#!/bin/bash
#
# Script for cropping images. The original images were generated using
# the script plot_book_elements.sh in DOLFIN in combination with
# fullscreen screenshots on my widescreen MacBook, so the images come
# out as widescreen which means some cropping is needed.
#
# Anders Logg, 2010-12-09

for f in originals/*.png; do
    g=`echo $f | cut -d'/' -f2`
    echo "Converting $g..."
    convert -crop 910x931+340+0 $f $g
done
