cmake_minimum_required(VERSION 3.14)
project(Hubbard-2leg-dynamics)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")
set(CMAKE_SOURCE_DIR "${PROJECT_SOURCE_DIR}/src/")

# Build type.
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

#if(NOT PHYSICAL_MODEL)
#    set(PHYSICAL_MODEL Fermion)
#endif()

if (${CMAKE_BUILD_TYPE} MATCHES Debug)
    set(CMAKE_VERBOSE_MAKEFILE ON)
endif ()

option(QLTEN_TIMING_MODE "Timing mode." OFF)
if (QLTEN_TIMING_MODE)
    add_definitions(-DQLTEN_TIMING_MODE)
endif ()

option(QLMPS_TIMING_MODE "Timing mode." OFF)
if (QLMPS_TIMING_MODE)
    add_definitions(-DQLMPS_TIMING_MODE)
endif ()

option(QLTEN_MPI_TIMING_MODE "Timing mode." OFF)
if (QLTEN_MPI_TIMING_MODE)
    add_definitions(-DQLTEN_MPI_TIMING_MODE)
endif ()

option(QLMPS_MPI_TIMING_MODE "Timing mode." OFF)
if (QLMPS_MPI_TIMING_MODE)
    add_definitions(-DQLMPS_MPI_TIMING_MODE)
endif ()

option(COMPILE_FIX_CODE "compile the fixing code" OFF)

option(HIGHER_SYMMETRY "1 or 2 particle conservation." OFF)
if (HIGHER_SYMMETRY)
    add_definitions(-DSYMMETRY_LEVLE=1)
else ()
    add_definitions(-DSYMMETRY_LEVLE=0)
endif ()
# Compilation and linking control.
set(CMAKE_CXX_STANDARD 17)
add_definitions(-Wall -g)


if (APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
    if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        set(MATH_LIB_COMPILE_FLAGS "-I$ENV{MKLROOT}/include")
        #Need test
        set(MATH_LIB_LINK_FLAGS $ENV{MKLROOT}/lib/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/libmkl_intel_thread.a $ENV{MKLROOT}/lib/libmkl_core.a -liomp5 -lpthread -lm -ldl)
    endif ()
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        set(MATH_LIB_COMPILE_FLAGS -m64 -I$ENV{MKLROOT}/include)
        # May not work
        set(MATH_LIB_LINK_FLAGS $ENV{MKLROOT}/lib/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/libmkl_intel_thread.a $ENV{MKLROOT}/lib/libmkl_core.a -L$ENV{MKLROOT}/lib -L$ENV{CMPLR_ROOT}/mac/compiler/lib/ -liomp5 -lpthread -lm -ldl)
    endif ()
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(MATH_LIB_COMPILE_FLAGS -m64 -I$ENV{MKLROOT}/include)
        # Note as of Intel oneAPI 2021.2, "source /opt/intel/oneapi/compiler/latest/env/vars.sh"
        set(MATH_LIB_LINK_FLAGS $ENV{MKLROOT}/lib/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/libmkl_intel_thread.a $ENV{MKLROOT}/lib/libmkl_core.a -L$ENV{MKLROOT}/lib -L$ENV{CMPLR_ROOT}/mac/compiler/lib/ -Wl, -rpath $ENV{CMPLR_ROOT}/mac/compiler/lib/libiomp5.dylib -liomp5 -lpthread -lm -ldl)
    endif ()
elseif (APPLE)
    #ARM
    add_definitions(-DUSE_OPENBLAS)
    set(BLA_VENDOR OpenBLAS)
    #        FIND_PACKAGE(BLAS REQUIRED)
    #        FIND_PACKAGE(LAPACK REQUIRED)
    set(OpenBLAS_ROOT "/opt/homebrew/opt/openblas/")
    set(Lapack_ROOT "/opt/homebrew/opt/lapack")
    set(OpenBLAS_INCLUDE_DIRS "${OpenBLAS_ROOT}/include")
    set(OpenBLAS_LIBRARIES "${OpenBLAS_ROOT}/lib/libblas.dylib")
    set(MATH_LIB_COMPILE_FLAGS -I${OpenBLAS_INCLUDE_DIRS} -pthread)
    set(MATH_LIB_LINK_FLAGS ${OpenBLAS_LIBRARIES} ${OpenBLAS_ROOT}/lib/liblapack.dylib -lm -lpthread -ldl -fopenmp -lclapack)
elseif (UNIX)
    #UNIX include APPLE, but we except it here
    if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
        set(MATH_LIB_COMPILE_FLAGS "-I$ENV{MKLROOT}/include")
        set(MATH_LIB_LINK_FLAGS -Wl,--start-group $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/intel64/libmkl_intel_thread.a $ENV{MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl)
    endif ()
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        set(MATH_LIB_COMPILE_FLAGS -m64 -I$ENV{MKLROOT}/include)
        # Link the Intel's OpenMP library to avoid performance issue when the library calls the MKL's gesdd function.
        set(MATH_LIB_LINK_FLAGS -Wl,--start-group $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/intel64/libmkl_intel_thread.a $ENV{MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -L$ENV{MKLROOT}/lib/intel64 -liomp5 -lpthread -lm -ldl)
    endif ()
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(MATH_LIB_COMPILE_FLAGS -m64 -I$ENV{MKLROOT}/include)
        set(MATH_LIB_LINK_FLAGS -Wl,--start-group $ENV{MKLROOT}/lib/intel64/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/intel64/libmkl_intel_thread.a $ENV{MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -L$ENV{MKLROOT}/lib/intel64 -liomp5 -lpthread -lm -ldl)
    endif ()
endif ()


find_path(TENSOR_HEADER_PATH "qlten")
find_path(MPS2_HEADER_PATH "qlmps")
find_package(hptt)
find_package(Boost
        1.74.0
        REQUIRED
        COMPONENTS serialization mpi)
find_package(MPI REQUIRED)
#set omp flag
if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(OMP_FLAGS -qopenmp)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(OMP_FLAGS -fopenmp)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(OMP_FLAGS -fopenmp)
endif ()

#find openmp?
add_compile_options(${MATH_LIB_COMPILE_FLAGS} ${OMP_FLAGS})

include_directories(${MPS2_HEADER_PATH}
        ${TENSOR_HEADER_PATH}
        ${hptt_INCLUDE_DIR}
        ${Boost_INCLUDE_DIR}
        ${MPI_CXX_HEADER_DIR}
)

link_libraries(
        ${hptt_LIBRARY}
        ${MATH_LIB_LINK_FLAGS}
        Boost::serialization
        Boost::mpi
        ${MPI_CXX_LINK_FLAGS}
        ${MPI_mpi_LIBRARY}
)

add_executable(mpogen src/mpogen.cpp src/operators.cpp)
add_executable(dmrg src/dmrg.cpp src/myutil.cpp src/operators.cpp)
add_executable(tdvp src/tdvp.cpp src/myutil.cpp src/operators.cpp)
add_executable(static_measure1 src/static_measure1.cpp src/operators.cpp)
#add_executable(measure2 src/measure2.cpp src/myutil.cpp src/operators.cpp)
