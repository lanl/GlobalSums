cmake_minimum_required (VERSION 3.0)
project (globalsums)

set (CMAKE_C_STANDARD 99)

# Set vectorization flags for a few compilers
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O3")
if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang") # using Clang
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fstrict-aliasing -fvectorize -fopenmp-simd\
       -march=native -mtune=native -Rpass-analysis=loop-vectorize")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU") # using GCC
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fstrict-aliasing -ftree-vectorize -fopenmp-simd\
       -march=native -mtune=native -mprefer-vector-width=512")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopt-info-vec-optimized -fopt-info-vec-missed\
      -fopt-info-loop-optimized -fopt-info-loop-missed")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel") # using Intel C
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ansi-alias -xHost -qopenmp-simd\
       -qopt-report=5 -qopt-report-phase=vec,loop")
   if ("${CMAKE_C_COMPILER_VERSION}" VERSION_GREATER "17.0.4")
     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qopt-zmm-usage=high")
   endif ("${CMAKE_C_COMPILER_VERSION}" VERSION_GREATER "17.0.4")
elseif (CMAKE_C_COMPILER_ID MATCHES "MSVC")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /Qvec-report:2")
elseif (CMAKE_C_COMPILER_ID MATCHES "XL")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qalias=restrict -qreport")
elseif (CMAKE_C_COMPILER_ID MATCHES "Cray")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -h restrict -h vector3 -h msgs -h negmsgs")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3")
if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang") # using Clang
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -fopenmp-simd")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fstrict-aliasing -fvectorize\
       -march=native -mtune=native -Rpass-analysis=loop-vectorize")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU") # using GCC
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -fopenmp-simd")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -DMAX_VECTOR_SIZE=512 -fstrict-aliasing -ftree-vectorize -fopenmp-simd\
       -march=native -mtune=native -mprefer-vector-width=512")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopt-info-vec-optimized -fopt-info-vec-missed\
      -fopt-info-loop-optimized -fopt-info-loop-missed")
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel") # using Intel C
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ansi-alias -xHost\
       -qopt-report=5 -qopt-report-phase=vec,loop -qopenmp -qopenmp-simd")
   if ("${CMAKE_C_COMPILER_VERSION}" VERSION_GREATER "17.0.4")
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopt-zmm-usage=high")
   endif ("${CMAKE_C_COMPILER_VERSION}" VERSION_GREATER "17.0.4")
elseif (CMAKE_C_COMPILER_ID MATCHES "MSVC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Qvec-report:2 -openmp:experimental")
elseif (CMAKE_C_COMPILER_ID MATCHES "XL")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qalias=restrict -qsimd=auto -qreport -qopenmp")
elseif (CMAKE_C_COMPILER_ID MATCHES "Cray")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -h restrict -h vector3 -h msgs -h negmsgs -h omp")
endif()

find_package(OpenMP)

try_compile(HAVE_AVX512 "${CMAKE_BINARY_DIR}" "${CMAKE_SOURCE_DIR}/test.c")

if(HAVE_AVX512)
    message("-- Detecting AVX512 hardware -- found")
    add_definitions(-DHAVE_AVX512)
else()
    message("-- Detecting AVX512 hardware -- not present")
endif()

# Adds build target of globalsums with source code files
add_executable(globalsums globalsums.c do_full_qdsum.c do_full_qdsum_wtrunc.c do_kahan_sum.c do_kahan_sum_gcc_v.c
               do_kahan_sum_omp.c do_kahan_sum_omp_wbittrunc.c do_kahan_sum_v.c do_knuth_sum.c do_knuth_sum_v.c
               do_ldsum.c do_ldsum_wbittrunc.c do_ldsum_wdigittrunc.c do_pair_sum.c do_qdsum.c do_qdsum_wtrunc.c
               do_sum.c do_sum_omp.c do_sum_omp_wbittrunc.c do_sum_wbittrunc.c do_sum_wdigittrunc.c
               kahan_fog_vector.cpp kahan_fog_vector8.cpp kahan_gcc_vector8.c kahan_gcc_vector8.c kahan_intel_vector8.c)

target_include_directories(globalsums PRIVATE vectorclass)
if(OPENMP_FOUND)
   set_target_properties(globalsums PROPERTIES COMPILE_FLAGS ${OpenMP_C_FLAGS})
   set_target_properties(globalsums PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
   set_target_properties(globalsums PROPERTIES LINK_FLAGS "${OpenMP_C_FLAGS}")
endif(OPENMP_FOUND)
target_link_libraries(globalsums quadmath m)

# Cleanup
add_custom_target(distclean COMMAND rm -rf CMakeCache.txt CMakeFiles
                  Makefile cmake_install.cmake globalsums.dSYM ipo_out.optrpt)
