cmake_minimum_required(VERSION 3.12)
project(xrotor_nographics_c)

set(CMAKE_CXX_STANDARD 14)
#find_package(Armadillo REQUIRED)
#include_directories(${ARMADILLO_INCLUDE_DIRS})
include_directories("include")

add_executable(xrotor_nographics_c main.cpp
        src/spline.cpp
        src/userio.cpp
        src/vortex.cpp
        src/xaero.cpp
        src/xio.cpp
        src/xrotor.cpp
        src/xutils.cpp)
#target_link_libraries(xrotor_nographics_c ${ARMADILLO_LIBRARIES})
