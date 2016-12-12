# Find matio library

include(FindPackageHandleStandardArgs)

find_path(MATIO_INCLUDE_DIRS matio.h)

set(MATIO_LIBRARY_DIRS "${MATIO_INCLUDE_DIRS}/../lib/")

find_library(MATIO_LIBRARIES matio PATHS ${MATIO_LIBRARY_DIRS})

# Handle the QUIETLY and REQUIRED arguments and set MATIO_FOUND to TRUE if all listed variables are TRUE
find_package_handle_standard_args(MATIO DEFAULT_MSG MATIO_INCLUDE_DIRS MATIO_LIBRARY_DIRS MATIO_LIBRARIES)

if(MATIO_FOUND)
    MESSAGE("-- Found matio package (MATIO)")
endif()

MESSAGE("-- Defined MATIO_INCLUDE_DIRS: ${MATIO_INCLUDE_DIRS}")
MESSAGE("-- Defined MATIO LIBRARY_DIRS: ${MATIO_LIBRARY_DIRS}")
MESSAGE("-- Defined MATIO_LIBRARIES: ${MATIO_LIBRARIES}")
