
include(FindPackageHandleStandardArgs)

find_path(CXXOPTS_INCLUDE_DIRS cxxopts.hpp HINTS ${THIRD_PARTY_DIR}/cxxopts/include)
find_path(TBB_INCLUDE_DIRS tbb/tbb.h
        HINTS ${TBB_INCLUDE_DIR} ${TBB_SEARCH_DIR}
        PATHS ${TBB_DEFAULT_SEARCH_DIR}
        PATH_SUFFIXES include)

find_package_handle_standard_args(CXXOPTS DEFAULT_MSG CXXOPTS_INCLUDE_DIRS)

if(CXXOPTS_FOUND)
    MESSAGE("-- Found cxxopts package (CXXOPTS)")
endif()

MESSAGE("-- Defined CXXOPTS_INCLUDE_DIRS: ${CXXOPTS_INCLUDE_DIRS}")
