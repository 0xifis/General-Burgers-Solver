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

add_executable(burgers main.cpp Model.cpp Burgers.cpp MyMPI.cpp MyMPI.h)

add_custom_target(burg
        COMMAND burgers OMP_NUM_THREADS=10 mpirun -np 1 burgers  --physics 1 0.5 1 0.02 --steps 2001 2001 4000 -x 1 1
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Running burgers in ${CMAKE_CURRENT_SOURCE_DIR}"
        SOURCES ${burgers_SOURCES}
        )

add_custom_target(advx
        COMMAND OMP_NUM_THREADS=10 mpirun -np 1 burgers  --physics 1 0 0 0 --steps 2001 2001 4000 -x 1 1
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Running burgers in ${CMAKE_BINARY_DIR}"
        SOURCES ${burgers_SOURCES}
        )

add_custom_target(advy
        COMMAND OMP_NUM_THREADS=10 mpirun -np 1 burgers  --physics 0 1 0 0 --steps 2001 2001 4000 -x 1 1
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Running burgers in ${CMAKE_BINARY_DIR}"
        SOURCES ${burgers_SOURCES}
        )

add_custom_target(diff
        COMMAND OMP_NUM_THREADS=10 mpirun -np 1 burgers  --physics 0 0 0 1 --steps 201 201 4000 -x 1 1
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Running burgers in ${CMAKE_BINARY_DIR}"
        SOURCES ${burgers_SOURCES}
        )


add_custom_target(burgp
        COMMAND OMP_NUM_THREADS=8 mpirun -np 2 burgers  --physics 1 0.5 1 0.02 --steps 2001 2001 4000 -x 2 1
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Running burgers in ${CMAKE_BINARY_DIR}"
        SOURCES ${burgers_SOURCES}
        )

add_custom_target(advxp
        COMMAND OMP_NUM_THREADS=8 mpirun -np 2 burgers  --physics 1 0 0 0 --steps 2001 2001 4000 -x 2 1
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Running burgers in ${CMAKE_BINARY_DIR}"
        SOURCES ${burgers_SOURCES}
        )

add_custom_target(advyp
        COMMAND OMP_NUM_THREADS=8 mpirun -np 2 burgers  --physics 0 1 0 0 --steps 2001 2001 4000 -x 2 1
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Running burgers in ${CMAKE_BINARY_DIR}"
        SOURCES ${burgers_SOURCES}
        )

add_custom_target(diffp
        COMMAND OMP_NUM_THREADS=8 mpirun -np 2 burgers  --physics 0 0 0 1 --steps 201 201 4000 -x 2 1
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Running burgers in ${CMAKE_BINARY_DIR}"
        SOURCES ${burgers_SOURCES}
        )


add_dependencies(burgp burgers)
add_dependencies(advxp burgers)
add_dependencies(advyp burgers)
add_dependencies(diffp burgers)

add_dependencies(burg burgers)
add_dependencies(advy burgers)
add_dependencies(advx burgers)
add_dependencies(diff burgers)


