#!/bin/sh
prog=decay_exper0.py

python $prog 1.25 0.75 0.5 0.1
if [ $? -eq 0 ]; then
image_types='png pdf'
for ext in $image_types; do
/bin/mv -i FE.$ext ../../fig-decay/FE4c.$ext
/bin/mv -i BE.$ext ../../fig-decay/BE4c.$ext
/bin/mv -i CN.$ext ../../fig-decay/CN4c.$ext
done
else
echo Could not run $prog successfully!
fi
