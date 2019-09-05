# Version 0.2 Increment by 0.1 every change

# Set vectorization flags for a few compilers
if(CMAKE_C_COMPILER_LOADED)
    if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang") # using Clang
        set(VECTOR_ALIASING_C_FLAGS "${VECTOR_ALIASING_C_FLAGS} -fstrict-aliasing")
        set(VECTOR_ARCH_C_FLAGS "${VECTOR_ARCH_C_FLAGS} -march=native -mtune=native")

        set(VECTOR_OPENMP_SIMD_C_FLAGS "${VECTOR_OPENMP_SIMD_C_FLAGS} -fopenmp-simd")
        set(VECTOR_C_FLAGS "${VECTOR_C_FLAGS} -fvectorize")
        set(VECTOR_NOVEC_C_FLAGS "${VECTOR_NOVEC_C_FLAGS} -fno-vectorize")
        set(VECTOR_C_VERBOSE "${VECTOR_C_VERBOSE} -Rpass=loop-vectorize -Rpass-missed=loop-vectorize -Rpass-analysis=loop-vectorize")

    elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU") # using GCC
        set(VECTOR_ALIASING_C_FLAGS "${VECTOR_ALIASING_C_FLAGS} -fstrict-aliasing")
        if ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
            set(VECTOR_ARCH_C_FLAGS "${VECTOR_ARCH_C_FLAGS} -march=native -mtune=native")
        elseif ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "ppc64le")
            set(VECTOR_ARCH_C_FLAGS "${VECTOR_ARCH_C_FLAGS} -mcpu=powerpc64le")
        elseif ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "aarch64")
            set(VECTOR_ARCH_C_FLAGS "${VECTOR_ARCH_C_FLAGS} -march=native -mtune=native")
        endif ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")

        set(VECTOR_OPENMP_SIMD_C_FLAGS "${VECTOR_OPENMP_SIMD_C_FLAGS} -fopenmp-simd")
        set(VECTOR_C_FLAGS "${VECTOR_C_FLAGS} -ftree-vectorize")
        if ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
            if ("${CMAKE_C_COMPILER_VERSION}" VERSION_GREATER "7.4.0")
                set(VECTOR_C_FLAGS "${VECTOR_C_FLAGS} -mprefer-vector-width=512")
            endif ("${CMAKE_C_COMPILER_VERSION}" VERSION_GREATER "7.4.0")
        endif ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")

        set(VECTOR_NOVEC_C_FLAGS "${VECTOR_NOVEC_C_FLAGS} -fno-tree-vectorize")
        set(VECTOR_C_VERBOSE "${VECTOR_C_VERBOSE} -fopt-info-vec-optimized -fopt-info-vec-missed -fopt-info-loop-optimized -fopt-info-loop-missed")

    elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel") # using Intel C
        set(VECTOR_ALIASING_C_FLAGS "${VECTOR_ALIASING_C_FLAGS} -ansi-alias")
        set(VECTOR_FPMODEL_C_FLAGS "${VECTOR_FPMODEL_C_FLAGS} -fp-model:precise")

        set(VECTOR_OPENMP_SIMD_C_FLAGS "${VECTOR_OPENMP_SIMD_C_FLAGS} -qopenmp-simd")
        set(VECTOR_C_FLAGS "${VECTOR_C_FLAGS} -xHOST")
        if ("${CMAKE_C_COMPILER_VERSION}" VERSION_GREATER "17.0.4")
            set(VECTOR_C_FLAGS "${VECTOR_C_FLAGS} -qopt-zmm-usage=high")
        endif ("${CMAKE_C_COMPILER_VERSION}" VERSION_GREATER "17.0.4")
        set(VECTOR_NOVEC_C_FLAGS "${VECTOR_NOVEC_C_FLAGS} -no-vec")
        set(VECTOR_C_VERBOSE "${VECTOR_C_VERBOSE} -qopt-report=5 -qopt-report-phase=openmp,loop,vec")

    elseif (CMAKE_C_COMPILER_ID MATCHES "PGI")
        set(VECTOR_ALIASING_C_FLAGS "${VECTOR_ALIASING_C_FLAGS} -alias=ansi")
        set(VECTOR_OPENMP_SIMD_C_FLAGS "${VECTOR_OPENMP_SIMD_C_FLAGS} -Mvect=simd")

        set(VECTOR_NOVEC_C_FLAGS "${VECTOR_NOVEC_C_FLAGS} -Mnovect ")
        set(VECTOR_C_VERBOSE "${VECTOR_C_VERBOSE} -Minfo=loop,inline,vect")

    elseif (CMAKE_C_COMPILER_ID MATCHES "MSVC")
        set(VECTOR_C_FLAGS "${VECTOR_C_FLAGS}" " ")

        set(VECTOR_NOVEC_C_FLAGS "${VECTOR_NOVEC_C_FLAGS}" " ")
        set(VECTOR_C_VERBOSE "${VECTOR_C_VERBOSE} -Qvec-report:2")

    elseif (CMAKE_C_COMPILER_ID MATCHES "XL")
        set(VECTOR_ALIASING_C_FLAGSS "${VECTOR_ALIASING_C_FLAGS} -qalias=restrict")
        set(VECTOR_FPMODEL_C_FLAGSS "${VECTOR_FPMODEL_C_FLAGS} -qstrict")
        set(VECTOR_ARCH_C_FLAGSS "${VECTOR_ARCH_C_FLAGS} -qhot -qarch=auto -qtune=auto")

        set(CMAKE_VEC_C_FLAGS "${CMAKE_VEC_FLAGS} -qsimd=auto")
        set(VECTOR_NOVEC_C_FLAGS "${VECTOR_NOVEC_C_FLAGS} -qsimd=noauto")
        # "long vector" optimizations
        #set(VECTOR_NOVEC_C_FLAGS "${VECTOR_NOVEC_C_FLAGS} -qhot=novector")
        set(VECTOR_C_VERBOSE "${VECTOR_C_VERBOSE} -qreport")

    elseif (CMAKE_C_COMPILER_ID MATCHES "Cray")
        set(VECTOR_ALIASING_C_FLAGS "${VECTOR_ALIASING_C_FLAGS} -h restrict=a")
        set(VECTOR_C_FLAGS "${VECTOR_C_FLAGS} -h vector=3")
  
       set(VECTOR_NOVEC_C_FLAGS "${VECTOR_NOVEC_C_FLAGS} -h vector=0")
       set(VECTOR_C_VERBOSE "${VECTOR_C_VERBOSE} -h msgs -h negmsgs -h list=a")

    endif()

    set(VECTOR_BASE_C_FLAGS "${VECTOR_ALIASING_C_FLAGS} ${VECTOR_ARCH_C_FLAGS} ${VECTOR_FPMODEL_C_FLAGS}")
    set(VECTOR_NOVEC_C_FLAGS "${VECTOR_BASE_C_FLAGS} ${VECTOR_NOVEC_C_FLAGS}")
    set(VECTOR_C_FLAGS "${VECTOR_BASE_C_FLAGS} ${VECTOR_C_FLAGS} ${VECTOR_OPENMP_SIMD_C_FLAGS}")

    mark_as_advanced(VECTOR_NOVEC_C_FLAGS VECTOR_C_FLAGS VECTOR_ALIASING_C_FLAGS
                    VECTOR_ARCH_C_FLAGS VECTOR_FPMODEL_C_FLAGS VECTOR_OPENMP_SIMD_C_FLAGS)

endif(CMAKE_C_COMPILER_LOADED)

#Start CXX Flags
if(CMAKE_CXX_COMPILER_LOADED)
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") # using Clang
        set(VECTOR_ALIASING_CXX_FLAGS "${VECTOR_ALIASING_CXX_FLAGS} -fstrict-aliasing")
        set(VECTOR_ARCH_CXX_FLAGS "${VECTOR_ARCH_CXX_FLAGS} -march=native -mtune=native")

        set(VECTOR_OPENMP_SIMD_CXX_FLAGS "${VECTOR_OPENMP_SIMD_CXX_FLAGS} -fopenmp-simd")
        set(VECTOR_CXX_FLAGS "${VECTOR_CXX_FLAGS} -fvectorize")
        set(VECTOR_NOVEC_CXX_FLAGS "${VECTOR_NOVEC_CXX_FLAGS} -fno-vectorize")
        set(VECTOR_CXX_VERBOSE "${VECTOR_CXX_VERBOSE} -Rpass=loop-vectorize -Rpass-missed=loop-vectorize -Rpass-analysis=loop-vectorize")

    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") # using GCC
        set(VECTOR_ALIASING_CXX_FLAGS "${VECTOR_ALIASING_CXX_FLAGS} -fstrict-aliasing")
        if ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
            set(VECTOR_ARCH_CXX_FLAGS "${VECTOR_ARCH_CXX_FLAGS} -march=native -mtune=native")
        elseif ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "ppc64le")
            set(VECTOR_ARCH_CXX_FLAGS "${VECTOR_ARCH_CXX_FLAGS} -mcpu=powerpc64le")
        elseif ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "aarch64")
            set(VECTOR_ARCH_CXX_FLAGS "${VECTOR_ARCH_CXX_FLAGS} -march=native -mtune=native")
        endif ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")

        set(VECTOR_OPENMP_SIMD_CXX_FLAGS "${VECTOR_OPENMP_SIMD_CXX_FLAGS} -fopenmp-simd")
        set(VECTOR_CXX_FLAGS "${VECTOR_CXX_FLAGS} -ftree-vectorize")
        if ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
            if ("${CMAKE_CXX_COMPILER_VERSION}" VERSION_GREATER "7.4.0")
                set(VECTOR_CXX_FLAGS "${VECTOR_CXX_FLAGS} -mprefer-vector-width=512")
            endif ("${CMAKE_CXX_COMPILER_VERSION}" VERSION_GREATER "7.4.0")
        endif ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")

        set(VECTOR_NOVEC_CXX_FLAGS "${VECTOR_NOVEC_CXX_FLAGS} -fno-tree-vectorize")
        set(VECTOR_CXX_VERBOSE "${VECTOR_CXX_VERBOSE} -fopt-info-vec-optimized -fopt-info-vec-missed -fopt-info-loop-optimized -fopt-info-loop-missed")

    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel") # using Intel C
        set(VECTOR_ALIASING_CXX_FLAGS "${VECTOR_ALIASING_CXX_FLAGS} -ansi-alias")
        set(VECTOR_FPMODEL_CXX_FLAGS "${VECTOR_FPMODEL_CXX_FLAGS} -fp-model:precise")

        set(VECTOR_OPENMP_SIMD_CXX_FLAGS "${VECTOR_OPENMP_SIMD_CXX_FLAGS} -qopenmp-simd")
        set(VECTOR_CXX_FLAGS "${VECTOR_CXX_FLAGS} -xHOST")
        if ("${CMAKE_CXX_COMPILER_VERSION}" VERSION_GREATER "17.0.4")
            set(VECTOR_CXX_FLAGS "${VECTOR_CXX_FLAGS} -qopt-zmm-usage=high")
        endif ("${CMAKE_CXX_COMPILER_VERSION}" VERSION_GREATER "17.0.4")
        set(VECTOR_NOVEC_CXX_FLAGS "${VECTOR_NOVEC_CXX_FLAGS} -no-vec")
        set(VECTOR_CXX_VERBOSE "${VECTOR_CXX_VERBOSE} -qopt-report=5 -qopt-report-phase=openmp,loop,vec")

    elseif (CMAKE_CXX_COMPILER_ID MATCHES "PGI")
        set(VECTOR_ALIASING_CXX_FLAGS "${VECTOR_ALIASING_CXX_FLAGS} -alias=ansi")
        set(VECTOR_OPENMP_SIMD_CXX_FLAGS "${VECTOR_OPENMP_SIMD_CXX_FLAGS} -Mvect=simd")

        set(VECTOR_NOVEC_CXX_FLAGS "${VECTOR_NOVEC_CXX_FLAGS} -Mnovect ")
        set(VECTOR_CXX_VERBOSE "${VECTOR_CXX_VERBOSE} -Minfo=loop,inline,vect")

    elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
        set(VECTOR_CXX_FLAGS "${VECTOR_CXX_FLAGS}" " ")

        set(VECTOR_NOVEC_CXX_FLAGS "${VECTOR_NOVEC_CXX_FLAGS}" " ")
        set(VECTOR_CXX_VERBOSE "${VECTOR_CXX_VERBOSE} -Qvec-report:2")

    elseif (CMAKE_CXX_COMPILER_ID MATCHES "XL")
        set(VECTOR_ALIASING_CXX_FLAGSS "${VECTOR_ALIASING_CXX_FLAGS} -qalias=restrict")
        set(VECTOR_FPMODEL_CXX_FLAGSS "${VECTOR_FPMODEL_CXX_FLAGS} -qstrict")
        set(VECTOR_ARCH_CXX_FLAGSS "${VECTOR_ARCH_CXX_FLAGS} -qhot -qarch=auto -qtune=auto")

        set(CMAKE_VEC_CXX_FLAGS "${CMAKE_VEC_FLAGS} -qsimd=auto")
        set(VECTOR_NOVEC_CXX_FLAGS "${VECTOR_NOVEC_CXX_FLAGS} -qsimd=noauto")
        # "long vector" optimizations
        #set(VECTOR_NOVEC_CXX_FLAGS "${VECTOR_NOVEC_CXX_FLAGS} -qhot=novector")
        set(VECTOR_CXX_VERBOSE "${VECTOR_CXX_VERBOSE} -qreport")

    elseif (CMAKE_CXX_COMPILER_ID MATCHES "Cray")
        set(VECTOR_ALIASING_CXX_FLAGS "${VECTOR_ALIASING_CXX_FLAGS} -h restrict=a")
        set(VECTOR_CXX_FLAGS "${VECTOR_CXX_FLAGS} -h vector=3")
  
        set(VECTOR_NOVEC_CXX_FLAGS "${VECTOR_NOVEC_CXX_FLAGS} -h vector=0")
        set(VECTOR_CXX_VERBOSE "${VECTOR_CXX_VERBOSE} -h msgs -h negmsgs -h list=a")

    endif()

    set(VECTOR_BASE_CXX_FLAGS "${VECTOR_ALIASING_CXX_FLAGS} ${VECTOR_ARCH_CXX_FLAGS} ${VECTOR_FPMODEL_CXX_FLAGS}")
    set(VECTOR_NOVEC_CXX_FLAGS "${VECTOR_BASE_CXX_FLAGS} ${VECTOR_NOVEC_CXX_FLAGS}")
    set(VECTOR_CXX_FLAGS "${VECTOR_BASE_CXX_FLAGS} ${VECTOR_CXX_FLAGS} ${VECTOR_OPENMP_SIMD_CXX_FLAGS}")

    mark_as_advanced(VECTOR_NOVEC_CXX_FLAGS VECTOR_CXX_FLAGS VECTOR_ALIASING_CXX_FLAGS
                     VECTOR_ARCH_CXX_FLAGS VECTOR_FPMODEL_CXX_FLAGS VECTOR_OPENMP_SIMD_CXX_FLAGS)
endif(CMAKE_CXX_COMPILER_LOADED)

# Start Fortran flags
if(CMAKE_Fortran_COMPILER_LOADED)
    if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Clang") # using Clang
        set(VECTOR_ALIASING_Fortran_FLAGS "${VECTOR_ALIASING_Fortran_FLAGS} -fstrict-aliasing")
        set(VECTOR_ARCH_Fortran_FLAGS "${VECTOR_ARCH_Fortran_FLAGS} -march=native -mtune=native")

        set(VECTOR_OPENMP_SIMD_Fortran_FLAGS "${VECTOR_OPENMP_SIMD_Fortran_FLAGS} -fopenmp-simd")
        set(VECTOR_Fortran_FLAGS "${VECTOR_Fortran_FLAGS} -fvectorize")
        set(VECTOR_NOVEC_Fortran_FLAGS "${VECTOR_NOVEC_Fortran_FLAGS} -fno-vectorize")
        set(VECTOR_Fortran_VERBOSE "${VECTOR_Fortran_VERBOSE} -Rpass=loop-vectorize -Rpass-missed=loop-vectorize -Rpass-analysis=loop-vectorize")

    elseif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU") # using GCC
        set(VECTOR_ALIASING_Fortran_FLAGS "${VECTOR_ALIASING_Fortran_FLAGS} -fstrict-aliasing")
        if ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
            set(VECTOR_ARCH_Fortran_FLAGS "${VECTOR_ARCH_Fortran_FLAGS} -march=native -mtune=native")
        elseif ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "ppc64le")
            set(VECTOR_ARCH_Fortran_FLAGS "${VECTOR_ARCH_Fortran_FLAGS} -mcpu=powerpc64le")
        elseif ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "aarch64")
            set(VECTOR_ARCH_Fortran_FLAGS "${VECTOR_ARCH_Fortran_FLAGS} -march=native -mtune=native")
        endif ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")

        set(VECTOR_OPENMP_SIMD_Fortran_FLAGS "${VECTOR_OPENMP_SIMD_Fortran_FLAGS} -fopenmp-simd")
        set(VECTOR_Fortran_FLAGS "${VECTOR_Fortran_FLAGS} -ftree-vectorize")
        if ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
            if ("${CMAKE_Fortran_COMPILER_VERSION}" VERSION_GREATER "7.4.0")
                set(VECTOR_Fortran_FLAGS "${VECTOR_Fortran_FLAGS} -mprefer-vector-width=512")
            endif ("${CMAKE_Fortran_COMPILER_VERSION}" VERSION_GREATER "7.4.0")
        endif ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")

        set(VECTOR_NOVEC_Fortran_FLAGS "${VECTOR_NOVEC_Fortran_FLAGS} -fno-tree-vectorize")
        set(VECTOR_Fortran_VERBOSE "${VECTOR_Fortran_VERBOSE} -fopt-info-vec-optimized -fopt-info-vec-missed -fopt-info-loop-optimized -fopt-info-loop-missed")

    elseif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel") # using Intel C
        set(VECTOR_ALIASING_Fortran_FLAGS "${VECTOR_ALIASING_Fortran_FLAGS} -ansi-alias")
        set(VECTOR_FPMODEL_Fortran_FLAGS "${VECTOR_FPMODEL_Fortran_FLAGS} -fp-model:precise")

        set(VECTOR_OPENMP_SIMD_Fortran_FLAGS "${VECTOR_OPENMP_SIMD_Fortran_FLAGS} -qopenmp-simd")
        set(VECTOR_Fortran_FLAGS "${VECTOR_Fortran_FLAGS} -xHOST")
        if ("${CMAKE_Fortran_COMPILER_VERSION}" VERSION_GREATER "17.0.4")
            set(VECTOR_Fortran_FLAGS "${VECTOR_Fortran_FLAGS} -qopt-zmm-usage=high")
        endif ("${CMAKE_Fortran_COMPILER_VERSION}" VERSION_GREATER "17.0.4")
        set(VECTOR_NOVEC_Fortran_FLAGS "${VECTOR_NOVEC_Fortran_FLAGS} -no-vec")
        set(VECTOR_Fortran_VERBOSE "${VECTOR_Fortran_VERBOSE} -qopt-report=5 -qopt-report-phase=openmp,loop,vec")

    elseif (CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
        set(VECTOR_ALIASING_Fortran_FLAGS "${VECTOR_ALIASING_Fortran_FLAGS} -alias=ansi")
        set(VECTOR_OPENMP_SIMD_Fortran_FLAGS "${VECTOR_OPENMP_SIMD_Fortran_FLAGS} -Mvect=simd")

        set(VECTOR_NOVEC_Fortran_FLAGS "${VECTOR_NOVEC_Fortran_FLAGS} -Mnovect ")
        set(VECTOR_Fortran_VERBOSE "${VECTOR_Fortran_VERBOSE} -Minfo=loop,inline,vect")

    elseif (CMAKE_Fortran_COMPILER_ID MATCHES "MSVC")
        set(VECTOR_Fortran_FLAGS "${VECTOR_Fortran_FLAGS}" " ")

        set(VECTOR_NOVEC_Fortran_FLAGS "${VECTOR_NOVEC_Fortran_FLAGS}" " ")
        set(VECTOR_Fortran_VERBOSE "${VECTOR_Fortran_VERBOSE} -Qvec-report:2")

    elseif (CMAKE_Fortran_COMPILER_ID MATCHES "XL")
        set(VECTOR_ALIASING_Fortran_FLAGSS "${VECTOR_ALIASING_Fortran_FLAGS} -qalias=restrict")
        set(VECTOR_FPMODEL_Fortran_FLAGSS "${VECTOR_FPMODEL_Fortran_FLAGS} -qstrict")
        set(VECTOR_ARCH_Fortran_FLAGSS "${VECTOR_ARCH_Fortran_FLAGS} -qhot -qarch=auto -qtune=auto")

        set(CMAKE_VEC_Fortran_FLAGS "${CMAKE_VEC_FLAGS} -qsimd=auto")
        set(VECTOR_NOVEC_Fortran_FLAGS "${VECTOR_NOVEC_Fortran_FLAGS} -qsimd=noauto")
        # "long vector" optimizations
        #set(VECTOR_NOVEC_Fortran_FLAGS "${VECTOR_NOVEC_Fortran_FLAGS} -qhot=novector")
        set(VECTOR_Fortran_VERBOSE "${VECTOR_Fortran_VERBOSE} -qreport")

    elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
        set(VECTOR_ALIASING_Fortran_FLAGS "${VECTOR_ALIASING_Fortran_FLAGS} -h restrict=a")
        set(VECTOR_Fortran_FLAGS "${VECTOR_Fortran_FLAGS} -h vector=3")
  
       set(VECTOR_NOVEC_Fortran_FLAGS "${VECTOR_NOVEC_Fortran_FLAGS} -h vector=0")
       set(VECTOR_Fortran_VERBOSE "${VECTOR_Fortran_VERBOSE} -h msgs -h negmsgs -h list=a")

    endif()


    set(VECTOR_BASE_Fortran_FLAGS "${VECTOR_ALIASING_Fortran_FLAGS} ${VECTOR_ARCH_Fortran_FLAGS} ${VECTOR_FPMODEL_Fortran_FLAGS}")
    set(VECTOR_NOVEC_Fortran_FLAGS "${VECTOR_BASE_Fortran_FLAGS} ${VECTOR_NOVEC_Fortran_FLAGS}")
    set(VECTOR_Fortran_FLAGS "${VECTOR_BASE_Fortran_FLAGS} ${VECTOR_Fortran_FLAGS} ${VECTOR_OPENMP_SIMD_Fortran_FLAGS}")

    mark_as_advanced(VECTOR_NOVEC_Fortran_FLAGS VECTOR_Fortran_FLAGS VECTOR_ALIASING_Fortran_FLAGS
                     VECTOR_ARCH_Fortran_FLAGS VECTOR_FPMODEL_Fortran_FLAGS VECTOR_OPENMP_SIMD_Fortran_FLAGS)
endif(CMAKE_Fortran_COMPILER_LOADED)

message(STATUS  "Setting Vector C Vector flags to -- ${VECTOR_C_FLAGS}")
message(STATUS  "Setting Vector C No-Vector flags to -- ${VECTOR_NOVEC_C_FLAGS}")
message(STATUS  "Setting Vector C Verbose flags to -- ${VECTOR_C_VERBOSE}")
