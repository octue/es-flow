
# Include google test
enable_testing()
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

# The unit test executables (there need to be more than one due to duplicate symbol conflicts)
add_executable(es_flow_test
        main_test.cpp
        unit/relations_test.cpp
        unit/readers_test.cpp
        unit/utilities_cumtrapz_test.cpp
        unit/utilities_trapz_test.cpp
        )

add_executable(es_flow_test_adem
        main_test.cpp
        unit/adem_test.cpp
        )

add_executable(es_flow_test_fit
        main_test.cpp
        unit/fit_test.cpp
        )

add_executable(es_flow_test_profile
        main_test.cpp
        unit/profile_test.cpp
        )

# Link the google test libraries
target_link_libraries(es_flow_test ${GTEST_BOTH_LIBRARIES})
target_link_libraries(es_flow_test_adem ${GTEST_BOTH_LIBRARIES})
target_link_libraries(es_flow_test_fit ${GTEST_BOTH_LIBRARIES})
target_link_libraries(es_flow_test_profile ${GTEST_BOTH_LIBRARIES})

# Link the third party libraries
target_link_libraries(es_flow_test ${CERES_LIBRARIES} ${MKL_LIBRARIES} ${MATIO_LIBRARIES})
target_link_libraries(es_flow_test_adem ${CERES_LIBRARIES} ${MKL_LIBRARIES} ${MATIO_LIBRARIES})
target_link_libraries(es_flow_test_fit ${CERES_LIBRARIES} ${MKL_LIBRARIES} ${MATIO_LIBRARIES})
target_link_libraries(es_flow_test_profile ${CERES_LIBRARIES} ${MKL_LIBRARIES} ${MATIO_LIBRARIES})

# Link the library to test
target_link_libraries(es_flow_test flow)
target_link_libraries(es_flow_test_adem flow)
target_link_libraries(es_flow_test_fit flow)
target_link_libraries(es_flow_test_profile flow)
