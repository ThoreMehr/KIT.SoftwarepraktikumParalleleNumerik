#!/bin/bash
filename=$1
base=$filename
base="${base%.*}"
echo "Compile $filename to $base"
nvcc -arch sm_20 --compiler-bindir=/usr/bin/g++-4.4 $filename -o $base -lm
