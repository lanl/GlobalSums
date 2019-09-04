#!/bin/sh
set -v
for i in clang/8.0.0 gcc/9.1.0 ThunderX2CN99/RHEL/7/gcc-8.2.0/armpl/19.2.0
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
