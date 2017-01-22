#!/bin/bash
filename=$1
base=$filename
base="${base%.*}"
echo "Compile $filename to $base"
nvcc --compiler-bindir=/usr/bin/g++-4.4 -use_fast_math -prec-div=false $filename -o $base -lm
