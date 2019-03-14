cmake_minimum_required(VERSION 3.5)
project(Coursework)

if(APPLE)
    SET(CMAKE_C_COMPILER /usr/local/bin/mpicc)
    SET(CMAKE_CXX_COMPILER /usr/local/bin/mpicxx)
else()
    SET(CMAKE_C_COMPILER /usr/bin/mpicc)
    SET(CMAKE_CXX_COMPILER /usr/bin/mpicxx)
endif()

add_compile_options(-Wall -Wextra -Wpedantic)

#SET(GCC_COVERAGE_LINK_FLAGS    "-llapack -lblas")
#SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")

SET(GCC_COVERAGE_COMPILE_FLAGS "-O3 -ffast-math -funroll-loops -march=native -ftree-vectorize -fopenmp")
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

set(CMAKE_CXX_STANDARD 14)

add_executable(burgers main.cpp Model.cpp Burgers.cpp)