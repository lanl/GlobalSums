#!/bin/sh
set -v
for i in clang/8.0.1 gcc/9.1.0 intel/19.0.4 pgi/18.10
do
    module purge
    module load cmake $i
    cmake .
    make 
    ./globalsums
    make clean
    make distclean
    module purge
done
