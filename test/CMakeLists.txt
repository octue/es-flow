
# Add ceres-solver dependency
find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS})

# List the unit test files
set(UNIT_TEST_SOURCE_FILES
        unit/adem_analysis_test.cpp
        unit/adem_signature_test.cpp
        unit/fit_profiles_test.cpp
        unit/instrument_analysis_test.cpp
        unit/profiles_test.cpp
        unit/window_test.cpp
        unit/readers_test.cpp ../source/utilities/cumtrapz.h)

# Include google test
enable_testing()
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# The unit test executable
add_executable(es_flow_test ${UNIT_TEST_SOURCE_FILES} ../source/utilities/cumtrapz.h)

# Link the google test libraries
target_link_libraries(es_flow_test ${GTEST_BOTH_LIBRARIES})

# Link the third party libraries
target_link_libraries(es_flow_test ${CERES_LIBRARIES} ${MKL_LIBRARIES} ${MATIO_LIBRARIES})

# Link the library to test
target_link_libraries(es_flow_test flow)
