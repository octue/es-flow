
include(FindPackageHandleStandardArgs)

# Add the Eigen NumericalIntegration library
#include(ExternalProject)
#ExternalProject_Add (
#        numint
#        PREFIX ${CMAKE_BINARY_DIR}/source/numint
#        GIT_REPOSITORY "https://github.com/tbs1980/NumericalIntegration"
#        GIT_TAG "master"
#        CONFIGURE_COMMAND ""
#        BUILD_COMMAND ""
#        INSTALL_COMMAND "" )
#
#ExternalProject_Get_Property(numint source_dir)
#find_path(NUMINT_INCLUDE_DIRS NumericalIntegration.h HINTS ${source_dir})

#find_package(NUMINT REQUIRED)
#include_directories(${NUMINT_INCLUDE_DIRS})

find_path(NUMINT_INCLUDE_DIRS NumericalIntegration.h HINTS ${THIRD_PARTY_DIR}/NumericalIntegration)

# Handle the QUIETLY and REQUIRED arguments and set MATIO_FOUND to TRUE if all listed variables are TRUE
find_package_handle_standard_args(NUMINT DEFAULT_MSG NUMINT_INCLUDE_DIRS)

if(NUMINT_FOUND)
    MESSAGE("-- Found NumericalIntegration package (NUMINT)")
endif()

MESSAGE("-- Defined NUMINT_INCLUDE_DIRS: ${NUMINT_INCLUDE_DIRS}")
