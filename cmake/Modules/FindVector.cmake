
# Set vectorization flags for a few compilers
if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang") # using Clang
    set(VECTOR_C_BASE_FLAGS "${VECTOR_C_BASE_FLAGS} -fstrict-aliasing -march=native -mtune=native")

    set(VECTOR_C_FLAGS "${VECTOR_C_FLAGS} ${VECTOR_C_BASE_FLAGS} -fvectorize -fopenmp-simd")
    set(VECTOR_NOVEC_C_FLAGS "${VECTOR_NOVEC_C_FLAGS} ${VECTOR_C_BASE_FLAGS} -fno-vectorize")
    set(VECTOR_C_VERBOSE "${VECTOR_C_VERBOSE} -Rpass=loop-vectorize -Rpass-missed=loop-vectorize -Rpass-analysis=loop-vectorize")

elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU") # using GCC
    set(VECTOR_C_BASE_FLAGS "${VECTOR_C_BASE_FLAGS} -fstrict-aliasing")
    if ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
        set(VECTOR_C_BASE_FLAGS "${VECTOR_C_BASE_FLAGS} -march=native -mtune=native")
    elseif ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "ppc64le")
        set(VECTOR_C_BASE_FLAGS "${VECTOR_C_BASE_FLAGS} -mcpu=powerpc64le")
    elseif ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "aarch64")
        set(VECTOR_C_BASE_FLAGS "${VECTOR_C_BASE_FLAGS} -march=native -mtune=native")
    endif ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")

    set(VECTOR_C_FLAGS "${VECTOR_C_FLAGS} ${VECTOR_C_BASE_FLAGS} -ftree-vectorize -fopenmp-simd")
    if ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
        if ("${CMAKE_C_COMPILER_VERSION}" VERSION_GREATER "7.4.0")
            set(VECTOR_C_FLAGS "${VECTOR_C_FLAGS} -mprefer-vector-width=512")
        endif ("${CMAKE_C_COMPILER_VERSION}" VERSION_GREATER "7.4.0")
    endif ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")

    set(VECTOR_NOVEC_C_FLAGS "${VECTOR_NOVEC_C_FLAGS} ${VECTOR_C_BASE_FLAGS} -fno-tree-vectorize")
    set(VECTOR_C_VERBOSE "${VECTOR_C_VERBOSE} -fopt-info-vec-optimized -fopt-info-vec-missed -fopt-info-loop-optimized -fopt-info-loop-missed")

elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel") # using Intel C
    set(VECTOR_C_BASE_FLAGS "${VECTOR_C_BASE_FLAGS} -ansi-alias -fp-model:precise")

    set(VECTOR_C_FLAGS "${VECTOR_C_FLAGS} ${VECTOR_C_BASE_FLAGS} -xHOST -qopenmp-simd")
    if ("${CMAKE_C_COMPILER_VERSION}" VERSION_GREATER "17.0.4")
        set(VECTOR_C_FLAGS "${VECTOR_C_FLAGS} -qopt-zmm-usage=high")
    endif ("${CMAKE_C_COMPILER_VERSION}" VERSION_GREATER "17.0.4")
    set(VECTOR_NOVEC_C_FLAGS "${VECTOR_NOVEC_C_FLAGS} ${VECTOR_C_BASE_FLAGS} -no-vec")
    set(VECTOR_C_VERBOSE "${VECTOR_C_VERBOSE} -qopt-report=5 -qopt-report-phase=openmp,loop,vec")

elseif (CMAKE_C_COMPILER_ID MATCHES "PGI")
    set(VECTOR_C_BASE_FLAGS "${VECTOR_C_BASE_FLAGS} -alias=ansi")
    set(VECTOR_C_FLAGS "${VECTOR_C_FLAGS} ${VECTOR_C_BASE_FLAGS} -Mvect=simd")

    set(VECTOR_NOVEC_C_FLAGS "${VECTOR_NOVEC_C_FLAGS} ${VECTOR_C_BASE_FLAGS} -Mnovect ")
    set(VECTOR_C_VERBOSE "${VECTOR_C_VERBOSE} -Minfo=loop,inline,vect")

elseif (CMAKE_C_COMPILER_ID MATCHES "MSVC")
    set(VECTOR_C_FLAGS "${VECTOR_C_FLAGS}" " ")

    set(VECTOR_NOVEC_C_FLAGS "${VECTOR_NOVEC_C_FLAGS}" " ")
    set(VECTOR_C_VERBOSE "${VECTOR_C_VERBOSE} -Qvec-report:2")

elseif (CMAKE_C_COMPILER_ID MATCHES "XL")
    set(VECTOR_C_BASE_FLAGSS "${VECTOR_C_BASE_FLAGS} -qalias=restrict -qstrict -qhot -qarch=auto -qtune=auto")
    set(VECTOR_NOVEC_C_FLAGS "${VECTOR_NOVEC_C_FLAGS} ${VECTOR_C_BASE_FLAGS} -qsimd=noauto")
    # "long vector" optimizations
    #set(VECTOR_NOVEC_C_FLAGS "${VECTOR_NOVEC_C_FLAGS} -qhot=novector")
    set(CMAKE_VEC_FLAGS "${CMAKE_VEC_FLAGS} ${VECTOR_C_BASE_FLAGS} -qsimd=auto")

    set(VECTOR_C_VERBOSE "${VECTOR_C_VERBOSE} -qreport")

elseif (CMAKE_C_COMPILER_ID MATCHES "Cray")
    set(VECTOR_C_BASE_FLAGS "${VECTOR_C_BASE_FLAGS} -h restrict=a")
    set(VECTOR_C_FLAGS "${VECTOR_C_FLAGS} ${VECTOR_C_BASE_FLAGS}-h vector=3")
  
   set(VECTOR_NOVEC_C_FLAGS "${VECTOR_NOVEC_C_FLAGS} ${VECTOR_C_BASE_FLAGS} -h vector=0")
   set(VECTOR_C_VERBOSE "${VECTOR_C_VERBOSE} -h msgs -h negmsgs -h list=a")

endif()


# Set vectorization flags for a few compilers
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") # using Clang
    set(VECTOR_CXX_BASE_FLAGS "${VECTOR_CXX_BASE_FLAGS} -fstrict-aliasing -march=native -mtune=native")

    set(VECTOR_CXX_FLAGS "${VECTOR_CXX_FLAGS} ${VECTOR_CXX_BASE_FLAGS} -fvectorize -fopenmp-simd")
    set(VECTOR_NOVEC_CXX_FLAGS "${VECTOR_NOVEC_CXX_FLAGS} ${VECTOR_CXX_BASE_FLAGS} -fno-vectorize")
    set(VECTOR_CXX_VERBOSE "${VECTOR_CXX_VERBOSE} -Rpass=loop-vectorize -Rpass-missed=loop-vectorize -Rpass-analysis=loop-vectorize")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") # using GCC
    set(VECTOR_CXX_BASE_FLAGS "${VECTOR_CXX_BASE_FLAGS} -fstrict-aliasing")
    if ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
        set(VECTOR_CXX_BASE_FLAGS "${VECTOR_CXX_BASE_FLAGS} -march=native -mtune=native")
    elseif ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "ppc64le")
        set(VECTOR_CXX_BASE_FLAGS "${VECTOR_CXX_BASE_FLAGS} -mcpu=powerpc64le")
    elseif ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "aarch64")
        set(VECTOR_CXX_BASE_FLAGS "${VECTOR_CXX_BASE_FLAGS} -march=native -mtune=native")
    endif ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")

    set(VECTOR_CXX_FLAGS "${VECTOR_CXX_FLAGS} ${VECTOR_CXX_BASE_FLAGS} -ftree-vectorize -fopenmp-simd")
    if ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
        if ("${CMAKE_CXX_COMPILER_VERSION}" VERSION_GREATER "7.4.0")
            set(VECTOR_CXX_FLAGS "${VECTOR_CXX_FLAGS} -mprefer-vector-width=512")
        endif ("${CMAKE_CXX_COMPILER_VERSION}" VERSION_GREATER "7.4.0")
    endif ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")

    set(VECTOR_NOVEC_CXX_FLAGS "${VECTOR_NOVEC_CXX_FLAGS} ${VECTOR_CXX_BASE_FLAGS} -fno-tree-vectorize")
    set(VECTOR_CXX_VERBOSE "${VECTOR_CXX_VERBOSE} -fopt-info-vec-optimized -fopt-info-vec-missed -fopt-info-loop-optimized -fopt-info-loop-missed")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel") # using Intel C
    set(VECTOR_CXX_BASE_FLAGS "${VECTOR_CXX_BASE_FLAGS} -ansi-alias -fp-model:precise")

    set(VECTOR_CXX_FLAGS "${VECTOR_CXX_FLAGS} ${VECTOR_CXX_BASE_FLAGS} -xHOST -qopenmp-simd")
    if ("${CMAKE_CXX_COMPILER_VERSION}" VERSION_GREATER "17.0.4")
        set(VECTOR_CXX_FLAGS "${VECTOR_CXX_FLAGS} -qopt-zmm-usage=high")
    endif ("${CMAKE_CXX_COMPILER_VERSION}" VERSION_GREATER "17.0.4")
    set(VECTOR_NOVEC_CXX_FLAGS "${VECTOR_NOVEC_CXX_FLAGS} ${VECTOR_CXX_BASE_FLAGS} -no-vec")
    set(VECTOR_CXX_VERBOSE "${VECTOR_CXX_VERBOSE} -qopt-report=5 -qopt-report-phase=openmp,loop,vec")

elseif (CMAKE_CXX_COMPILER_ID MATCHES "PGI")
    set(VECTOR_CXX_BASE_FLAGS "${VECTOR_CXX_BASE_FLAGS} -alias=ansi")
    set(VECTOR_CXX_FLAGS "${VECTOR_CXX_FLAGS} ${VECTOR_CXX_BASE_FLAGS} -Mvect=simd")

    set(VECTOR_NOVEC_CXX_FLAGS "${VECTOR_NOVEC_CXX_FLAGS} ${VECTOR_CXX_BASE_FLAGS} -Mnovect ")
    set(VECTOR_CXX_VERBOSE "${VECTOR_CXX_VERBOSE} -Minfo=loop,inline,vect")

elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    set(VECTOR_CXX_FLAGS "${VECTOR_CXX_FLAGS}" " ")

    set(VECTOR_NOVEC_CXX_FLAGS "${VECTOR_NOVEC_CXX_FLAGS}" " ")
    set(VECTOR_CXX_VERBOSE "${VECTOR_CXX_VERBOSE} -Qvec-report:2")

elseif (CMAKE_CXX_COMPILER_ID MATCHES "XL")
    set(VECTOR_CXX_BASE_FLAGSS "${VECTOR_CXX_BASE_FLAGS} -qalias=restrict -qstrict -qhot -qarch=auto -qtune=auto")
    set(VECTOR_NOVEC_CXX_FLAGS "${VECTOR_NOVEC_CXX_FLAGS} ${VECTOR_CXX_BASE_FLAGS} -qsimd=noauto")
    # "long vector" optimizations
    #set(VECTOR_NOVEC_CXX_FLAGS "${VECTOR_NOVEC_CXX_FLAGS} -qhot=novector")
    set(CMAKE_VEC_FLAGS "${CMAKE_VEC_FLAGS} ${VECTOR_CXX_BASE_FLAGS} -qsimd=auto")

    set(VECTOR_CXX_VERBOSE "${VECTOR_CXX_VERBOSE} -qreport")

elseif (CMAKE_CXX_COMPILER_ID MATCHES "Cray")
    set(VECTOR_CXX_BASE_FLAGS "${VECTOR_CXX_BASE_FLAGS} -h restrict=a")
    set(VECTOR_CXX_FLAGS "${VECTOR_CXX_FLAGS} ${VECTOR_CXX_BASE_FLAGS}-h vector=3")
  
   set(VECTOR_NOVEC_CXX_FLAGS "${VECTOR_NOVEC_CXX_FLAGS} ${VECTOR_CXX_BASE_FLAGS} -h vector=0")
   set(VECTOR_CXX_VERBOSE "${VECTOR_CXX_VERBOSE} -h msgs -h negmsgs -h list=a")

endif()

