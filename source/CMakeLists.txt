
# Add ceres-solver dependency
find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})

# Add Intel MKL dependency
set(MKL_MULTI_THREADED "1")
find_package(MKL REQUIRED)
include_directories(${MKL_INCLUDE_DIRS})

# Add matio dependency
find_package(MATIO REQUIRED)
include_directories(${MATIO_INCLUDE_DIRS})


# Add the flow library
set(SOURCE_FILES
        adem_analysis.cpp
        adem_analysis.h
        adem_signature.cpp
        adem_signature.h
        fit_profiles.cpp
        fit_profiles.h
        instrument_analysis.cpp
        instrument_analysis.h
        profiles_velocity.cpp
        profiles_velocity.h
        profiles_spectrum.cpp
        profiles_spectrum.h
        profiles_stress.cpp
        profiles_stress.h
        window.cpp window.h
        profile.cpp profile.h)
add_library(flow STATIC ${SOURCE_FILES})

# The hello world function which gives an executable target for the library
add_executable(es_flow main.cpp)
target_link_libraries(es_flow ${CERES_LIBRARIES} ${MKL_LIBRARIES} ${MATIO_LIBRARIES})
