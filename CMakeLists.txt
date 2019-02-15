cmake_minimum_required(VERSION 3.6)
project(es_flow)


# Add cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")


# Set make options
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
set(CMAKE_VERBOSE_MAKEFILE ON)


# Add directories where third party and other tools exist
set(THIRD_PARTY_DIR ${CMAKE_SOURCE_DIR}/../thirdparty)
message("-- Defined THIRD_PARTY_DIR: ${THIRD_PARTY_DIR}")

# Add argument parsing tool
find_package(CXXOPTS REQUIRED)
include_directories(${CXXOPTS_INCLUDE_DIRS})

# Add plotting tool
find_package(CPPLOT REQUIRED)
include_directories(${CPPLOT_INCLUDE_DIRS})

# Add ceres-solver dependency
find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})

# Add numerical integration for eigen
find_package(NUMINT REQUIRED)
include_directories(${NUMINT_INCLUDE_DIRS})

# Add Intel MKL dependency
set(MKL_MULTI_THREADED "1")
find_package(MKL REQUIRED)
include_directories(${MKL_INCLUDE_DIRS})

# Add matio dependency
find_package(MATIO REQUIRED)
include_directories(${MATIO_INCLUDE_DIRS})

#Add libboost (mainly for string and path manipulation)
find_package(Boost REQUIRED)

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
MESSAGE("-- Using dynamic library extension (platform dependent): ${CMAKE_SHARED_LIBRARY_SUFFIX}")
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


# Build the documentation, if we can find sphinx
set(DOCS_INSTALL_DIR "share/doc/es-flow" CACHE STRING "es-flow docs install path relative to CMAKE_INSTALL_PREFIX")
find_package(Sphinx QUIET)
if (NOT SPHINX_FOUND)
    message(WARNING "-- Failed to find Sphinx, disabling build of documentation.")
else()
    # Generate the User's Guide (html).
    # The corresponding target is es_flow_docs, but is included in ALL.
    message("-- Found Sphinx")
    add_subdirectory(docs)
endif()


# Add source to include_directories so the test routines can pick up the library contents
include_directories(source)


# Add the souce and test directories
add_subdirectory(source)
add_subdirectory(test)


# Report included directories
get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
message(STATUS "Include directories:")
foreach(dir ${dirs})
    message(STATUS "    '${dir}'")
endforeach()
