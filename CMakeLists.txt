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

# Add Intel TBB dependency
find_package(TBB REQUIRED)
set(USE_TBB "1")
include_directories(${TBB_INCLUDE_DIRS})

# Add Intel MKL dependency
set(MKL_MULTI_THREADED "1")
set(MKL_THREADING "TBB")
find_package(MKL REQUIRED)
include_directories(${MKL_INCLUDE_DIR})

# Ensure any eigen operations use MKL
add_definitions(-DEIGEN_USE_MKL_ALL)

# Add matio dependency
find_package(MATIO REQUIRED)
include_directories(${MATIO_INCLUDE_DIRS})

# Add libboost (mainly for string and path manipulation)
find_package(Boost REQUIRED)


# Build the documentation
add_subdirectory(docs)


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
