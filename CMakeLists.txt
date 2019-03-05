cmake_minimum_required(VERSION 3.13)
project(Coursework)

add_compile_options(-Wall -Wextra -Wpedantic)

set(CMAKE_CXX_STANDARD 14)

find_package( blas )

find_package( blaze )
if( blaze_FOUND )
    add_library( blaze_target INTERFACE )
    target_link_libraries( blaze_target INTERFACE blaze::blaze )
endif()

add_executable(burgers main.cpp Model.cpp Burgers.cpp)