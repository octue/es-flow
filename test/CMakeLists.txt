
# Include google test
enable_testing()

# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" . RESULT_VARIABLE result WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build . RESULT_VARIABLE result WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src ${CMAKE_CURRENT_BINARY_DIR}/googletest-build EXCLUDE_FROM_ALL)
include_directories("${gtest_SOURCE_DIR}/include")

# The unit test executables (there need to be more than one due to duplicate symbol conflicts)
add_executable(es_flow_test
        main_test.cpp
        unit/relations_test.cpp
        unit/readers_test.cpp
        unit/utilities_conv_test.cpp
        unit/utilities_cumtrapz_test.cpp
        unit/utilities_cumsum_test.cpp
        unit/utilities_integration_test.cpp
        unit/utilities_interp_test.cpp
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
target_link_libraries(es_flow_test gtest gtest_main)
target_link_libraries(es_flow_test_adem gtest gtest_main)
target_link_libraries(es_flow_test_fit gtest gtest_main)
target_link_libraries(es_flow_test_profile gtest gtest_main)

# Link the third party libraries
target_link_libraries(es_flow_test ${CERES_LIBRARIES} ${MKL_LIBRARIES} ${TBB_LIBRARIES} ${MATIO_LIBRARIES})
target_link_libraries(es_flow_test_adem ${CERES_LIBRARIES} ${MKL_LIBRARIES} ${TBB_LIBRARIES} ${MATIO_LIBRARIES})
target_link_libraries(es_flow_test_fit ${CERES_LIBRARIES} ${MKL_LIBRARIES} ${TBB_LIBRARIES} ${MATIO_LIBRARIES})
target_link_libraries(es_flow_test_profile ${CERES_LIBRARIES} ${MKL_LIBRARIES} ${TBB_LIBRARIES} ${MATIO_LIBRARIES})

# Link the library to test
target_link_libraries(es_flow_test flow)
target_link_libraries(es_flow_test_adem flow)
target_link_libraries(es_flow_test_fit flow)
target_link_libraries(es_flow_test_profile flow)
