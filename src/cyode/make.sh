#!/bin/sh
if [ $# -eq 2 ]; then
  problem=$1
  version=$2
else
  echo "Usage:   ./make.sh problem version"
  echo "Example: ./make.sh 1 1"
  exit 1
fi
python setup.py $problem $version build_ext --inplace
# Test if module works:
python -c "import ode${problem}_cy${version}"
if [ $? -eq 0 ]; then
  echo "Module ode${problem}_cy${version} successfully built"
fi
# Compile and view C code:
# cython -a ode$problem_cy$version.pyx