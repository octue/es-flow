cmake_minimum_required(VERSION 3.6)
project(es_flow)

# Add cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

# Set make options
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_VERBOSE_MAKEFILE ON)


# Variables used for building mex files, platform independent
# TODO shift to using find_package for MATLAB MEX https://cmake.org/cmake/help/v3.3/module/FindMatlab.html
set(MATLAB_LIBS
        -lemlrt
        -lcovrt
        -lut
        -lmx
        -lmex
        -lmat
        -lmwmathutil)
if(${CMAKE_SHARED_LIBRARY_SUFFIX} STREQUAL ".dylib")
    set(ARCH maci64)
    set(MATLAB_ROOT /Applications/MATLAB_R2016a.app)
    set(MATLAB_MEXEXT .mexmaci64)
    set(MATLAB_INCLUDE_DIRS ${MATLAB_ROOT}/extern/include)
    set(MATLAB_LIB_DIRS
            ${MATLAB_ROOT}/bin/maci64
            ${MATLAB_ROOT}/sys/os/maci64
            )
elseif(${CMAKE_SHARED_LIBRARY_SUFFIX} STREQUAL ".so")
    set(ARCH linux)
    # TODO set MATLABROOT
    set(MATLAB_MEXEXT .mexa64)
    # TODO set LIBS location
    set(MATLAB_LIB_DIRS
            ${MATLAB_ROOT}/bin/x86_64
            ${MATLAB_ROOT}/sys/os/x86_64
            )
elseif(${CMAKE_SHARED_LIBRARY_SUFFIX} STREQUAL ".dll")
    set(ARCH win64)
    # TODO set MATLABROOT
    set(MATLAB_MEXEXT .mexw64)
    # TODO set LIBS location - useful make file as top answer at http://stackoverflow.com/questions/8552580/using-gcc-mingw-as-matlabs-mex-compiler/9451392#9451392
    set(MATLAB_LIB_DIRS
            ${MATLAB_ROOT}/bin/x86_64
            ${MATLAB_ROOT}/sys/os/x86_64
            )
else()
    message(ERROR "Cannot identify architecture using system library extension when defining MATLAB libraries")
endif()


# Add include directories so the test routines can pick up the library contents
include_directories(source)

# Add the souce and test directories
add_subdirectory(source)
add_subdirectory(test)
