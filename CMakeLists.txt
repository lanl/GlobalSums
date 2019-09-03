cmake_minimum_required (VERSION 3.0)
project (globalsums)

set (CMAKE_C_STANDARD 99)

# Set optional intrinsic support
if ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
   if (NOT CMAKE_C_COMPILER_ID MATCHES "PGI")
      set(HAVE_X86_64_INTRINSICS true)
   endif (NOT CMAKE_C_COMPILER_ID MATCHES "PGI")
endif ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang") # using Clang
   set(HAVE_GCC_VECTOR_EXTENSIONS true)

elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU") # using GCC
   set(HAVE_QUADMATH true)
   set(HAVE_GCC_VECTOR_EXTENSIONS true)
   if ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
      set(HAVE_FOG_VECTOR_CLASS true)
   endif ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")

elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel") # using Intel C
   set(HAVE_QUADMATH true)
   set(HAVE_FOG_VECTOR_CLASS true)
   set(HAVE_GCC_VECTOR_EXTENSIONS true)

elseif (CMAKE_C_COMPILER_ID MATCHES "PGI")
   #set(HAVE_FOG_VECTOR_CLASS true)

endif()

# Set vectorization flags for a few compilers
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O3")

if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang") # using Clang
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fstrict-aliasing -march=native -mtune=native")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Rpass=loop-vectorize -Rpass-missed=loop-vectorize -Rpass-analysis=loop-vectorize")
   set(CMAKE_NOVEC_FLAGS "${CMAKE_NOVEC_FLAGS} -fno-vectorize")
   set(CMAKE_VEC_FLAGS "${CMAKE_VEC_FLAGS} -fvectorize")

elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU") # using GCC
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fstrict-aliasing")
   if ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -mtune=native")
   elseif ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "ppc64le")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mcpu=powerpc64le")
   endif ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopt-info-vec-optimized -fopt-info-vec-missed\
      -fopt-info-loop-optimized -fopt-info-loop-missed")
   set(CMAKE_NOVEC_FLAGS "${CMAKE_NOVEC_FLAGS} -fno-tree-vectorize")
   set(CMAKE_VEC_FLAGS "${CMAKE_VEC_FLAGS} -ftree-vectorize")
   if ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
      if ("${CMAKE_C_COMPILER_VERSION}" VERSION_GREATER "7.4.0")
         set(CMAKE_VEC_FLAGS "${CMAKE_VEC_FLAGS} -mprefer-vector-width=512")
      endif ("${CMAKE_C_COMPILER_VERSION}" VERSION_GREATER "7.4.0")
   endif ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")

elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel") # using Intel C
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ansi-alias -fp-model:precise")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qopt-report=5 -qopt-report-phase=openmp,loop,vec")
   set(CMAKE_NOVEC_FLAGS "${CMAKE_NOVEC_FLAGS} -no-vec")
   set(CMAKE_VEC_FLAGS "${CMAKE_VEC_FLAGS} -xHOST")
   if ("${CMAKE_C_COMPILER_VERSION}" VERSION_GREATER "17.0.4")
     set(CMAKE_VEC_FLAGS "${CMAKE_VEC_FLAGS} -qopt-zmm-usage=high")
   endif ("${CMAKE_C_COMPILER_VERSION}" VERSION_GREATER "17.0.4")

elseif (CMAKE_C_COMPILER_ID MATCHES "PGI")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -alias=ansi")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Minfo=loop,inline,vect")
   set(CMAKE_NOVEC_FLAGS "${CMAKE_NOVEC_FLAGS} -Mnovect ")
   set(CMAKE_VEC_FLAGS "${CMAKE_VEC_FLAGS} -Mvect=simd")

elseif (CMAKE_C_COMPILER_ID MATCHES "MSVC")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Qvec-report:2")
   set(CMAKE_NOVEC_FLAGS "${CMAKE_NOVEC_FLAGS}" " ")
   set(CMAKE_VEC_FLAGS "${CMAKE_VEC_FLAGS}" " ")

elseif (CMAKE_C_COMPILER_ID MATCHES "XL")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qalias=restrict -qstrict -qhot -qarch=auto -qtune=auto")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -qreport")
   set(CMAKE_NOVEC_FLAGS "${CMAKE_NOVEC_FLAGS} -qsimd=noauto")
   # "long vector" optimizations
   #set(CMAKE_NOVEC_FLAGS "${CMAKE_NOVEC_FLAGS} -qhot=novector")
   set(CMAKE_VEC_FLAGS "${CMAKE_VEC_FLAGS} -qsimd=auto")

elseif (CMAKE_C_COMPILER_ID MATCHES "Cray")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -h restrict=a")
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -h msgs -h negmsgs -h list=a")
   set(CMAKE_NOVEC_FLAGS "${CMAKE_NOVEC_FLAGS} -h vector=0")
   set(CMAKE_VEC_FLAGS "${CMAKE_VEC_FLAGS} -h vector=3")

endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3")
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") # using Clang
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fstrict-aliasing -fvectorize -march=native -mtune=native")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Rpass=loop-vectorize -Rpass-missed=loop-vectorize -Rpass-analysis=loop-vectorize") 

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") # using GCC
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMAX_VECTOR_SIZE=512 -fstrict-aliasing -ftree-vectorize -fopenmp-simd ")
   if ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mtune=native")
   elseif ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "ppc64le")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mcpu=powerpc64le")
   endif ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopt-info-vec-optimized -fopt-info-vec-missed\
      -fopt-info-loop-optimized -fopt-info-loop-missed")
   if ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
      if ("${CMAKE_CXX_COMPILER_VERSION}" VERSION_GREATER "7.4.0")
         set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mprefer-vector-width=512")
      endif ("${CMAKE_CXX_COMPILER_VERSION}" VERSION_GREATER "7.4.0")
   endif ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel") # using Intel C
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMAX_VECTOR_SIZE=512 -fp-model:precise -ansi-alias -xHost\
       -qopt-report=5 -qopt-report-phase=vec,loop")
   if ("${CMAKE_CXX_COMPILER_VERSION}" VERSION_GREATER "17.0.4")
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopt-zmm-usage=high")
   endif ("${CMAKE_CXX_COMPILER_VERSION}" VERSION_GREATER "17.0.4")

elseif (CMAKE_CXX_COMPILER_ID MATCHES "PGI")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -alias=ansi")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Minfo=loop,inline,vect")

elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Qvec-report:2")

elseif (CMAKE_CXX_COMPILER_ID MATCHES "XL")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qalias=restrict -qstrict -qhot -qarch=auto -qtune=auto")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qreport")

elseif (CMAKE_CXX_COMPILER_ID MATCHES "Cray")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -h restrict -h vector3 -h msgs -h negmsgs")

endif()

find_package(OpenMP)

try_compile(HAVE_AVX512 "${CMAKE_BINARY_DIR}" "${CMAKE_SOURCE_DIR}/test.c")

if(HAVE_AVX512)
    message("-- Detecting AVX512 hardware -- found")
    add_definitions(-DHAVE_AVX512)
else()
    message("-- Deecting AVX512 hardware -- not present")
endif()

# Adds build target of globalsums with source code files
set (FILE_LIST globalsums.c
               do_kahan_sum.c
               serial_sum_novec.c
               do_kahan_sum_omp.c
               do_kahan_sum_omp_wbittrunc.c
               do_knuth_sum.c
               do_ldsum.c
               do_ldsum_wbittrunc.c
               do_ldsum_wdigittrunc.c
               do_pair_sum.c
               do_sum.c
               do_sum_omp.c
               do_sum_omp_wbittrunc.c
               do_sum_wbittrunc.c
               do_sum_wdigittrunc.c)

if (HAVE_FOG_VECTOR_CLASS)
   add_definitions(-DHAVE_FOG_VECTOR_CLASS)
   set (FILE_LIST ${FILE_LIST} 
               kahan_fog_vector.cpp
               kahan_fog_vector8.cpp
               knuth_fog_vector.cpp
               knuth_fog_vector8.cpp)
endif (HAVE_FOG_VECTOR_CLASS)

if (HAVE_X86_64_INTRINSICS)
   add_definitions(-DHAVE_X86_64_INTRINSICS)
   set (FILE_LIST ${FILE_LIST} 
               kahan_intel_vector.c
               kahan_intel_vector8.c
               knuth_intel_vector.c
               knuth_intel_vector8.c)
endif (HAVE_X86_64_INTRINSICS)

if (HAVE_GCC_VECTOR_EXTENSIONS)
   add_definitions(-DHAVE_GCC_VECTOR_EXTENSIONS)
   set (FILE_LIST ${FILE_LIST} 
               kahan_gcc_vector.c
               kahan_gcc_vector8.c
               knuth_gcc_vector.c
               knuth_gcc_vector8.c)
endif (HAVE_GCC_VECTOR_EXTENSIONS)

if (HAVE_QUADMATH)
   add_definitions(-DHAVE_QUADMATH)
   set (FILE_LIST ${FILE_LIST} 
               do_full_qdsum.c
               do_full_qdsum_wtrunc.c
               do_qdsum.c
               do_qdsum_wtrunc.c)
endif (HAVE_QUADMATH)

add_executable(globalsums ${FILE_LIST})

target_include_directories(globalsums PRIVATE vectorclass)
set_source_files_properties(globalsums.c PROPERTIES COMPILE_FLAGS ${CMAKE_VEC_FLAGS})
set_source_files_properties(do_full_qdsum.c PROPERTIES COMPILE_FLAGS ${CMAKE_VEC_FLAGS})
set_source_files_properties(do_full_qdsum_wtrunc.c PROPERTIES COMPILE_FLAGS ${CMAKE_VEC_FLAGS})
set_source_files_properties(do_kahan_sum.c PROPERTIES COMPILE_FLAGS ${CMAKE_VEC_FLAGS})
set_source_files_properties(serial_sum_novec.c PROPERTIES COMPILE_FLAGS ${CMAKE_NOVEC_FLAGS})
set_source_files_properties(do_kahan_sum_omp.c PROPERTIES COMPILE_FLAGS ${CMAKE_VEC_FLAGS})
set_source_files_properties(do_kahan_sum_omp_wbittrunc.c PROPERTIES COMPILE_FLAGS ${CMAKE_VEC_FLAGS})
set_source_files_properties(do_knuth_sum.c PROPERTIES COMPILE_FLAGS ${CMAKE_VEC_FLAGS})
set_source_files_properties(do_ldsum.c PROPERTIES COMPILE_FLAGS ${CMAKE_VEC_FLAGS})
set_source_files_properties(do_ldsum_wbittrunc.c PROPERTIES COMPILE_FLAGS ${CMAKE_VEC_FLAGS})
set_source_files_properties(do_ldsum_wdigittrunc.c PROPERTIES COMPILE_FLAGS ${CMAKE_VEC_FLAGS})
set_source_files_properties(do_pair_sum.c PROPERTIES COMPILE_FLAGS ${CMAKE_VEC_FLAGS})
set_source_files_properties(do_qdsum.c PROPERTIES COMPILE_FLAGS ${CMAKE_VEC_FLAGS})
set_source_files_properties(do_qdsum_wtrunc.c PROPERTIES COMPILE_FLAGS ${CMAKE_VEC_FLAGS})
set_source_files_properties(do_sum.c PROPERTIES COMPILE_FLAGS ${CMAKE_VEC_FLAGS})
set_source_files_properties(do_sum_omp.c PROPERTIES COMPILE_FLAGS ${CMAKE_VEC_FLAGS})
set_source_files_properties(do_sum_omp_wbittrunc.c PROPERTIES COMPILE_FLAGS ${CMAKE_VEC_FLAGS})
set_source_files_properties(do_sum_wbittrunc.c PROPERTIES COMPILE_FLAGS ${CMAKE_VEC_FLAGS})
set_source_files_properties(do_sum_wdigittrunc.c PROPERTIES COMPILE_FLAGS ${CMAKE_VEC_FLAGS})
set_source_files_properties(kahan_fog_vector.cpp PROPERTIES COMPILE_FLAGS ${CMAKE_VEC_FLAGS})
set_source_files_properties(kahan_fog_vector8.cpp PROPERTIES COMPILE_FLAGS ${CMAKE_VEC_FLAGS})
set_source_files_properties(kahan_gcc_vector.c PROPERTIES COMPILE_FLAGS ${CMAKE_VEC_FLAGS})
set_source_files_properties(kahan_gcc_vector8.c PROPERTIES COMPILE_FLAGS ${CMAKE_VEC_FLAGS})
set_source_files_properties(kahan_intel_vector.c PROPERTIES COMPILE_FLAGS ${CMAKE_VEC_FLAGS})
set_source_files_properties(kahan_intel_vector8.c PROPERTIES COMPILE_FLAGS ${CMAKE_VEC_FLAGS})
set_source_files_properties(knuth_fog_vector.cpp PROPERTIES COMPILE_FLAGS ${CMAKE_VEC_FLAGS})
set_source_files_properties(knuth_fog_vector8.cpp PROPERTIES COMPILE_FLAGS ${CMAKE_VEC_FLAGS})
set_source_files_properties(knuth_gcc_vector.c PROPERTIES COMPILE_FLAGS ${CMAKE_VEC_FLAGS})
set_source_files_properties(knuth_gcc_vector8.c PROPERTIES COMPILE_FLAGS ${CMAKE_VEC_FLAGS})
set_source_files_properties(knuth_intel_vector.c PROPERTIES COMPILE_FLAGS ${CMAKE_VEC_FLAGS})
set_source_files_properties(knuth_intel_vector8.c PROPERTIES COMPILE_FLAGS ${CMAKE_VEC_FLAGS})
if(OPENMP_FOUND)
   set_target_properties(globalsums PROPERTIES COMPILE_FLAGS ${OpenMP_C_FLAGS})
   set_target_properties(globalsums PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
   set_target_properties(globalsums PROPERTIES LINK_FLAGS "${OpenMP_C_FLAGS}")
endif(OPENMP_FOUND)
if (NOT CMAKE_C_COMPILER_ID MATCHES "XL" AND NOT CMAKE_C_COMPILER_ID MATCHES "PGI")
target_link_libraries(globalsums quadmath)
endif (NOT CMAKE_C_COMPILER_ID MATCHES "XL" AND NOT CMAKE_C_COMPILER_ID MATCHES "PGI")
target_link_libraries(globalsums m)

# Cleanup
add_custom_target(distclean COMMAND rm -rf CMakeCache.txt CMakeFiles
                  Makefile cmake_install.cmake globalsums.dSYM ipo_out.optrpt)
