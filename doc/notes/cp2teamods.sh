#!/bin/sh
names="cyode"
for name in $names; do
dest=../../../hplgit.github.com/teamods/$name
echo cp -r fig-$name main_$name.html main_$name.pdf $name-sphinx $dest
cp -r fig-$name main_$name.html main_$name.pdf $name-sphinx $dest
done
