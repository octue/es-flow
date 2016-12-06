
# List the unit test files
set(UNIT_TEST_SOURCE_FILES
        unit/adem_analysis_test.cpp
        unit/fit_profiles_test.cpp
        unit/profiles_velocity_test.cpp
        unit/adem_signature_test.cpp
        unit/profiles_spectrum_test.cpp
        unit/profiles_stress_test.cpp
        unit/window_test.cpp)

# Include google test
enable_testing()
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# The unit test executable
add_executable(es_flow_test ${UNIT_TEST_SOURCE_FILES})

# Include the es_flow library headers
set(ES_INCLUDE_DIRS
        source/adem
        source/fit
        source/instrument
        source/profiles)
include_directories(${ES_INCLUDE_DIRS})

# Link the google test libraries
target_link_libraries(es_flow_test ${GTEST_BOTH_LIBRARIES})

# Link the library to test
target_link_libraries(es_flow_test flow)
#add_test(AllTestsInFoo foo)