
set(SOURCE_FILES
        profiles/Bins.cpp
        profiles/Profile.cpp
        profiles/SpectrumProfile.cpp
        profiles/StressProfile.cpp
        )

add_library(flow STATIC ${SOURCE_FILES})

#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/codegen/mex/biotsavart ${CMAKE_CURRENT_SOURCE_DIR}/codegen/mex/biotsavart/interface ${MATLAB_INCLUDE_DIRS})
#
#link_directories(${CMAKE_CURRENT_SOURCE_DIR}/codegen/mex/biotsavart ${CMAKE_CURRENT_SOURCE_DIR}/codegen/mex/biotsavart/interface ${MATLAB_LIB_DIRS})
