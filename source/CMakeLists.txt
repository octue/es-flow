
# Add argument parsing tool
find_package(CXXOPTS REQUIRED)
include_directories(${CXXOPTS_INCLUDE_DIRS})

# Add ceres-solver dependency
find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})

# Add argument parsing tool
find_package(NumericalIntegration REQUIRED)
include_directories(${NUMINT_INCLUDE_DIRS})

# Add Intel MKL dependency
set(MKL_MULTI_THREADED "1")
find_package(MKL REQUIRED)
include_directories(${MKL_INCLUDE_DIRS})

# Add matio dependency
find_package(MATIO REQUIRED)
include_directories(${MATIO_INCLUDE_DIRS})

# Add include directories so the test routines can pick up the library contents
include_directories(relations)

# Add the souce and test directories
add_subdirectory(relations)

# Add the flow library
set(SOURCE_FILES
        adem_analysis.cpp
        adem_analysis.h
        adem_signature.cpp
        adem_signature.h
        data_types.h
        fit_profiles.cpp
        fit_profiles.h
        instrument_analysis.cpp
        instrument_analysis.h
        relations/velocity.cpp
        relations/velocity.h
        profile.cpp profile.h
        readers.cpp readers.h
        variable_readers.h
        window.cpp window.h utilities/cumtrapz.h)
add_library(flow STATIC ${SOURCE_FILES} utilities/cumtrapz.h)

# The hello world function which gives an executable target for the library
add_executable(es_flow main.cpp utilities/cumtrapz.h)
target_link_libraries(es_flow ${CERES_LIBRARIES} ${MKL_LIBRARIES} ${MATIO_LIBRARIES})
