
include(FindPackageHandleStandardArgs)

find_path(CXXOPTS_INCLUDE_DIRS cxxopts.hpp PATHS ${THIRD_PARTY_DIR}/cxxopts/include)

# Handle the QUIETLY and REQUIRED arguments and set MATIO_FOUND to TRUE if all listed variables are TRUE
find_package_handle_standard_args(CXXOPTS DEFAULT_MSG CXXOPTS_INCLUDE_DIRS)

if(CXXOPTS_FOUND)
    MESSAGE("-- Found cxxopts package (CXXOPTS)")
endif()

MESSAGE("-- Defined CXXOPTS_INCLUDE_DIRS: ${CXXOPTS_INCLUDE_DIRS}")
