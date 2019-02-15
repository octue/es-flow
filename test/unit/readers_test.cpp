/*
 * READERS_TEST.CPP Test fixtures for file readers
 *
 * Author:                   Tom Clark (thclark @ github)
 *
 * Copyright (c) 2016-9 Octue Ltd. All Rights Reserved.
 *
 */

#include <errno.h>
#include "gtest/gtest.h"
#include <stdio.h>
#include <unistd.h>

#include "data_types.h"
#include "readers.h"


using namespace es;


class ReaderTest : public ::testing::Test {
public:

    std::string data_path;

protected:

    virtual void SetUp() {

        // Get the test data directory path from the environment variable
        if(const char* env_p = std::getenv("TEST_DATA_DIR")) {
            std::cout << "TEST_DATA_DIR = " << env_p << std::endl;
            data_path = env_p;
        } else {
            throw std::invalid_argument("Invalid environment variable 'TEST_DATA_DIR'");
        }
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
