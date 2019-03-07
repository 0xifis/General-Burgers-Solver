cmake_minimum_required(VERSION 3.5)
project(Coursework)

add_compile_options(-Wall -Wextra -Wpedantic)

set(CMAKE_CXX_STANDARD 14)

add_executable(burgers main.cpp Model.cpp Burgers.cpp)

find_package( BLAS REQUIRED )
target_link_libraries(burgers ${BLAS_LIBRARIES})

include_directories(./includes)
find_package( blaze REQUIRED )
add_library( blaze_target INTERFACE )
target_link_libraries( blaze_target INTERFACE blaze::blaze )

target_link_libraries( burgers ${BLAZE_LIBRARIES} )