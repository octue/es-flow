/*
 * READERS_TEST.CPP Test fixtures for file readers
 *
 * References:
 *
 *   [1]
 *
 * Future Improvements:
 *
 *   [1] Possibly useful to parameterise tests. To see how, check out minute 9 onward in
 *   https://www.youtube.com/watch?v=16FI1-d2P4E&t=20s
 *
 * Author:                   T. Clark
 * Work address:             Ocean Array Systems Ltd
 *                           Hauser Forum
 *                           3 Charles Babbage Road
 *                           Cambridge
 *                           CB3 0GT
 * Email:                    tom.clark@oceanarraysystems.com
 * Website:                  www.oceanarraysystems.com
 *
 * Copyright (c) 2016-17 Ocean Array Systems, All Rights Reserved.
 *
 */
#include <unistd.h>
#include <stdio.h>
#include <errno.h>
#include "gtest/gtest.h"
#include "readers.h"
#include "data_types.h"

using namespace es;

class ReaderTest : public ::testing::Test {
public:

    std::string data_path;

protected:

    virtual void SetUp() {
        // Code here will be called immediately after the constructor (right
        // before each test).

        // Get the test data directory path from the environmnet variable
        if(const char* env_p = std::getenv("TEST_DATA_DIR")) {
            std::cout << std::endl << "Setting up ReaderTest" <<std::endl;
            std::cout << "TEST_DATA_DIR = " << env_p << std::endl;
            data_path = env_p;
        } else {
            throw std::invalid_argument("Invalid environment variable 'TEST_DATA_DIR'");
        }
    }

    virtual void TearDown() {
        // Code here will be called immediately after each test (right
        // before the destructor).
        std::cout << "Tearing down ReaderTest()..." << std::endl << std::endl;
    }

};

TEST_F(ReaderTest, test_lidar_basic){

    // Get lidar_basic test file
    std::string file = data_path + std::string("/fake_lidar_basic.mat");

    // Construct a lidar data reader
    Reader<BasicLidar> lr(file);

    // Check the file is of a standard type
    lr.checkFileType(OAS_STANDARD);

    // Read the test file into a BasicLidar data type object
    auto data = lr.read();

    // Apply windowing settings for timeseries processing (default 10 minute windows, 50% overlap)
    lr.checkTimeseries(STRICTLY_MONOTONIC);
    double window_duration = 10*60;
    double overlap = 0.5;
    lr.setWindowDuration(window_duration);
//    lr.setWindowOverlap(window_overlap);


    // Ensure that the printing operator does not error
    std::cout << lr << std::endl;


}

class AnalysisTest : public ::testing::Test {

protected:

    virtual void SetUp() {
        std::cout << std::endl << "Setting up AnalysisTest()..." << std::endl;
    }

    virtual void TearDown() {
        std::cout << "Tearing down AnalysisTest()..." << std::endl << std::endl;
    }

};

// Unit tests for the profile class
TEST_F(AnalysisTest, test_construct_double_profile){




}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
