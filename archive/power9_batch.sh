#!/bin/sh
set -v
for i in clang/8.0.0 gcc/9.1.0 ibm/xlc-16.1.1.3-xlf-16.1.1.3 pgi/19.3
do
    module purge
    module load cmake $i
    cmake .
    #cmake -DCMAKE_VECTOR_VERBOSE .
    make 
    ./globalsums
    make clean
    make distclean
    module purge
done
