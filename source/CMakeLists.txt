
# Add include directories so the test routines can pick up the library contents
include_directories(relations)

# Add the souce and test directories
add_subdirectory(relations)

# Add the flow library
set(SOURCE_FILES
        adem/analysis.cpp
        adem/analysis.h
        adem/signature.cpp
        adem/signature.h
        data_types.h
        fit_profiles.cpp
        fit_profiles.h
        instrument_analysis.cpp
        instrument_analysis.h
        relations/velocity.cpp
        relations/velocity.h
        relations/stress.h
        profile.cpp profile.h
        readers.cpp readers.h
        variable_readers.h
        window.cpp window.h
        utilities/cumtrapz.h)
add_library(flow STATIC ${SOURCE_FILES})


# The hello world function which gives an executable target for the library
add_executable(es_flow main.cpp)
target_link_libraries(es_flow ${CERES_LIBRARIES} ${MKL_LIBRARIES} ${MATIO_LIBRARIES})

