#!/bin/bash
for i in $(ls *.cu); do
	filename=$(basename "$i")
	sh compile.sh $filename
done
