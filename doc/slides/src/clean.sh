#!/bin/sh
# Remove all files that make_lectures.sh can recreate.

if [ ! -d lec-* ]; then
  echo 'Run clean.sh from a subdirectory for a topic'
  exit 1
fi

doconce clean
rm -rf automake_sphinx.py