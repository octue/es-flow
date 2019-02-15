
include(FindPackageHandleStandardArgs)
MESSAGE(${THIRD_PARTY_DIR}../../thclark/cpplot/source)
find_path(CPPLOT_INCLUDE_DIRS cpplot.h  HINTS ${THIRD_PARTY_DIR}/../../thclark/cpplot/source/)

# Handle the QUIETLY and REQUIRED arguments and set <NAME>_FOUND to TRUE if all listed variables are TRUE
find_package_handle_standard_args(CPPLOT DEFAULT_MSG CPPLOT_INCLUDE_DIRS)

if(CPPLOT_FOUND)
    MESSAGE("-- Found cpplot package (CPPLOT)")
endif()

MESSAGE("-- Defined CPPLOT_INCLUDE_DIRS: ${CPPLOT_INCLUDE_DIRS}")
