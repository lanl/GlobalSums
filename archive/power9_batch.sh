#!/bin/sh
set -v
for i in clang/8.0.1 gcc/9.2.0 ibm/xlc-16.1.1.3-xlf-16.1.1.3 pgi/19.7
do
    module purge
    module load cmake $i
    cmake .
    #cmake -DCMAKE_VECTOR_VERBOSE=1 .
    make 
    ./globalsums
    make clean
    make distclean
    module purge
done
