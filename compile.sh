#!/bin/bash

rm read
g++ geometry.C read.C analysis.C main.C `root-config --libs --cflags` -o read
#new comment
#compile for memory test
#g++ read.C `root-config --libs --cflags` -o read  -ltcmalloc
