#!/bin/bash
time -p bash -c "for ((i=0;i<10;i++));do $1 > /dev/null;done;"
