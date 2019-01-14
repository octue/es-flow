
include(FindPackageHandleStandardArgs)

find_path(NUMINT_INCLUDE_DIRS NumericalIntegration.h HINTS ${THIRD_PARTY_DIR}/NumericalIntegration/include)

# Handle the QUIETLY and REQUIRED arguments and set MATIO_FOUND to TRUE if all listed variables are TRUE
find_package_handle_standard_args(NUMINT DEFAULT_MSG NUMINT_INCLUDE_DIRS)

if(NUMINT_FOUND)
    MESSAGE("-- Found NumericalIntegration package (NUMINT)")
endif()

MESSAGE("-- Defined NUMINT_INCLUDE_DIRS: ${NUMINT_INCLUDE_DIRS}")
