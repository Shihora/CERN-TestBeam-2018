#!/bin/bash

rm read
g++ geometry.C read.C analysis.C main.C -rpath ${ROOTSYS}/lib `root-config --libs --cflags` -o read 
#compile for memory test
#g++ read.C `root-config --libs --cflags` -o read  -ltcmalloc
