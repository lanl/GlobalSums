cmake_minimum_required (VERSION 3.0)
project (globalsums)

set (CMAKE_C_STANDARD 99)

# Set default optimization level
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3")

set(CMAKE_MODULE_PATH cmake/Modules)

find_package(Vector)
#message("VECTOR_C_FLAGS is ${VECTOR_C_FLAGS}")
#message("VECTOR_C_VERBOSE is ${VECTOR_C_VERBOSE}")
#message("VECTOR_NOVEC_C_FLAGS is ${VECTOR_NOVEC_C_FLAGS}")
#message("VECTOR_ALIASING_C_FLAGS is ${VECTOR_ALIASING_FLAGS}")
#message("VECTOR_ARCH_C_FLAGS is ${VECTOR_ARCH_C_FLAGS}")
#message("VECTOR_FPMODEL_C_FLAGS is ${VECTOR_FPMODEL_C_FLAGS}")
#message("VECTOR_OPENMP_SIMD_C_FLAGS is ${VECTOR_OPENMP_SIMD_C_FLAGS}")
find_package(OpenMP)

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

#Pushing vector flags onto compiler flags for try compiles
set(CMAKE_C_FLAGS_SAVE ${CMAKE_C_FLAGS})
set(CMAKE_CXX_FLAGS_SAVE ${CMAKE_CXX_FLAGS})
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${VECTOR_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${VECTOR_C_FLAGS}")

try_compile(HAVE_FOG_VECTOR_CLASS "${CMAKE_BINARY_DIR}" "${CMAKE_SOURCE_DIR}/fogvectorclasstest.cpp")

if(HAVE_FOG_VECTOR_CLASS)
    message(STATUS "Trying Fog Vector Class -- works")
    add_definitions(-DHAVE_FOG_VECTOR_CLASS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMAX_VECTOR_SIZE=512")
    set (FILE_LIST ${FILE_LIST} 
               kahan_fog_vector.cpp
               kahan_fog_vector8.cpp
               knuth_fog_vector.cpp
               knuth_fog_vector8.cpp)
else()
    message(STATUS "Trying Fog Vector Class -- fails")
endif()


try_compile(HAVE_X86_64_INTRINSICS "${CMAKE_BINARY_DIR}" "${CMAKE_SOURCE_DIR}/x86_64_intrinsics_test.c")

if(HAVE_X86_64_INTRINSICS)
    message(STATUS "Trying X86_64 intrinsics -- works")
    add_definitions(-DHAVE_X86_64_INTRINSICS -DMAX_VECTOR_SIZE=512)
    set (FILE_LIST ${FILE_LIST} 
               kahan_intel_vector.c
               kahan_intel_vector8.c
               knuth_intel_vector.c
               knuth_intel_vector8.c)

    try_compile(HAVE_AVX512 "${CMAKE_BINARY_DIR}" "${CMAKE_SOURCE_DIR}/AVX512test.c")

    if(HAVE_AVX512)
        message(STATUS "Detecting AVX512 hardware -- found")
        add_definitions(-DHAVE_AVX512)
    else()
        message(STATUS "Detecting AVX512 hardware -- not present")
    endif()

else()
    message(STATUS "Trying X86_64 intrinsics -- fails")
endif()


try_compile(HAVE_GCC_VECTOR_EXTENSIONS "${CMAKE_BINARY_DIR}" "${CMAKE_SOURCE_DIR}/gcc_vector_extrinsics_test.c")

if(HAVE_GCC_VECTOR_EXTENSIONS)
    message(STATUS "Trying GCC Vector Extensions -- works")
    add_definitions(-DHAVE_GCC_VECTOR_EXTENSIONS)
    set (FILE_LIST ${FILE_LIST} 
               kahan_gcc_vector.c
               kahan_gcc_vector8.c
               knuth_gcc_vector.c
               knuth_gcc_vector8.c)
else()
    message(STATUS "Trying GCC Vector Extensions -- fails")
endif (HAVE_GCC_VECTOR_EXTENSIONS)


try_compile(HAVE_QUADMATH "${CMAKE_BINARY_DIR}" "${CMAKE_SOURCE_DIR}/quadmath_test.c")

if(HAVE_QUADMATH)
    message(STATUS "Trying Quad math library -- works")
    add_definitions(-DHAVE_QUADMATH)
    set (FILE_LIST ${FILE_LIST} 
               do_full_qdsum.c
               do_full_qdsum_wtrunc.c
               do_qdsum.c
               do_qdsum_wtrunc.c)
else()
    message(STATUS "Trying Quad math library -- fails")
endif (HAVE_QUADMATH)

#Popping compiler flags off and back to original state
set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS_SAVE})
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_SAVE})

add_executable(globalsums ${FILE_LIST})

if (CMAKE_VECTOR_VERBOSE)
    set(VECTOR_C_FLAGS "${VECTOR_C_FLAGS} ${VECTOR_C_VERBOSE}")
    set(VECTOR_CXX_FLAGS "${VECTOR_CXX_FLAGS} ${VECTOR_CXX_VERBOSE}")
endif (CMAKE_VECTOR_VERBOSE)

target_include_directories(globalsums PRIVATE vectorclass)
set_source_files_properties(globalsums.c PROPERTIES COMPILE_FLAGS ${VECTOR_C_FLAGS})
set_source_files_properties(do_full_qdsum.c PROPERTIES COMPILE_FLAGS ${VECTOR_C_FLAGS})
set_source_files_properties(do_full_qdsum_wtrunc.c PROPERTIES COMPILE_FLAGS ${VECTOR_C_FLAGS})
set_source_files_properties(do_kahan_sum.c PROPERTIES COMPILE_FLAGS ${VECTOR_C_FLAGS})
set_source_files_properties(serial_sum_novec.c PROPERTIES COMPILE_FLAGS ${VECTOR_NOVEC_C_FLAGS})
set_source_files_properties(do_kahan_sum_omp.c PROPERTIES COMPILE_FLAGS ${VECTOR_C_FLAGS})
set_source_files_properties(do_kahan_sum_omp_wbittrunc.c PROPERTIES COMPILE_FLAGS ${VECTOR_C_FLAGS})
set_source_files_properties(do_knuth_sum.c PROPERTIES COMPILE_FLAGS ${VECTOR_C_FLAGS})
set_source_files_properties(do_ldsum.c PROPERTIES COMPILE_FLAGS ${VECTOR_C_FLAGS})
set_source_files_properties(do_ldsum_wbittrunc.c PROPERTIES COMPILE_FLAGS ${VECTOR_C_FLAGS})
set_source_files_properties(do_ldsum_wdigittrunc.c PROPERTIES COMPILE_FLAGS ${VECTOR_C_FLAGS})
set_source_files_properties(do_pair_sum.c PROPERTIES COMPILE_FLAGS ${VECTOR_C_FLAGS})
set_source_files_properties(do_qdsum.c PROPERTIES COMPILE_FLAGS ${VECTOR_C_FLAGS})
set_source_files_properties(do_qdsum_wtrunc.c PROPERTIES COMPILE_FLAGS ${VECTOR_C_FLAGS})
set_source_files_properties(do_sum.c PROPERTIES COMPILE_FLAGS ${VECTOR_C_FLAGS})
set_source_files_properties(do_sum_omp.c PROPERTIES COMPILE_FLAGS ${VECTOR_C_FLAGS})
set_source_files_properties(do_sum_omp_wbittrunc.c PROPERTIES COMPILE_FLAGS ${VECTOR_C_FLAGS})
set_source_files_properties(do_sum_wbittrunc.c PROPERTIES COMPILE_FLAGS ${VECTOR_C_FLAGS})
set_source_files_properties(do_sum_wdigittrunc.c PROPERTIES COMPILE_FLAGS ${VECTOR_C_FLAGS})
set_source_files_properties(kahan_fog_vector.cpp PROPERTIES COMPILE_FLAGS ${VECTOR_CXX_FLAGS})
set_source_files_properties(kahan_fog_vector8.cpp PROPERTIES COMPILE_FLAGS ${VECTOR_CXX_FLAGS})
set_source_files_properties(kahan_gcc_vector.c PROPERTIES COMPILE_FLAGS ${VECTOR_C_FLAGS})
set_source_files_properties(kahan_gcc_vector8.c PROPERTIES COMPILE_FLAGS ${VECTOR_C_FLAGS})
set_source_files_properties(kahan_intel_vector.c PROPERTIES COMPILE_FLAGS ${VECTOR_C_FLAGS})
set_source_files_properties(kahan_intel_vector8.c PROPERTIES COMPILE_FLAGS ${VECTOR_C_FLAGS})
set_source_files_properties(knuth_fog_vector.cpp PROPERTIES COMPILE_FLAGS ${VECTOR_CXX_FLAGS})
set_source_files_properties(knuth_fog_vector8.cpp PROPERTIES COMPILE_FLAGS ${VECTOR_CXX_FLAGS})
set_source_files_properties(knuth_gcc_vector.c PROPERTIES COMPILE_FLAGS ${VECTOR_C_FLAGS})
set_source_files_properties(knuth_gcc_vector8.c PROPERTIES COMPILE_FLAGS ${VECTOR_C_FLAGS})
set_source_files_properties(knuth_intel_vector.c PROPERTIES COMPILE_FLAGS ${VECTOR_C_FLAGS})
set_source_files_properties(knuth_intel_vector8.c PROPERTIES COMPILE_FLAGS ${VECTOR_C_FLAGS})
if(OPENMP_FOUND)
   set_target_properties(globalsums PROPERTIES COMPILE_FLAGS ${OpenMP_C_FLAGS})
   set_target_properties(globalsums PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
   set_target_properties(globalsums PROPERTIES LINK_FLAGS "${OpenMP_C_FLAGS}")
endif(OPENMP_FOUND)
if (HAVE_QUADMATH)
target_link_libraries(globalsums quadmath)
endif (HAVE_QUADMATH)

target_link_libraries(globalsums m)

# Cleanup
add_custom_target(distclean COMMAND rm -rf CMakeCache.txt CMakeFiles
                  Makefile cmake_install.cmake globalsums.dSYM ipo_out.optrpt)
