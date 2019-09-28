#!/bin/sh
set -v
#for i in clang/8.0.1 gcc/9.2.0 intel/19.0.5 pgi/19.7
for i in intel/19.0.5
do
    module purge
    module load cmake $i
    cmake -DCMAKE_VECTOR_VERBOSE=1 .
    make 
    ./globalsums
    make clean
    make distclean
    module purge
done
