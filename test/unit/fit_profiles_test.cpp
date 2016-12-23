/*
 * FIT_PROFILES_TEST.CPP Test fixtures for bins, profiles and profile fitting methods.
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

#include "gtest/gtest.h"
#include "profile.h"

using namespace es;

// Test fixture for the Bin class
class BinTest : public ::testing::Test {

protected:

    virtual void SetUp() {
        // Code here will be called immediately after the constructor (right
        // before each test).
        std::cout << std::endl << "Setting up BinTest()..." << std::endl;
    }

    virtual void TearDown() {
        // Code here will be called immediately after each test (right
        // before the destructor).
        std::cout << "Tearing down BinTest()..." << std::endl;
    }

};

// Unit tests for the Bin class
TEST_F(BinTest, test_construct_bins){

    // Construct 10 bins spaced equally 1 m apart from 1m above reference using the base constructor
    std::vector<double> z;
    for (int i=0; i<10; i++) z.push_back(double(i+1));
    Bins b1 = Bins(z);

    // Ensure that the printing operator does not error
    std::cout << b1 << std::endl;

    // Check value
    EXPECT_EQ(z[9], double(10));

}

// Test fixture for the Profile class
class ProfileTest : public ::testing::Test {

protected:

    virtual void SetUp() {
        std::cout << std::endl << "Setting up BinTest()..." << std::endl;
    }

    virtual void TearDown() {
        std::cout << "Tearing down BinTest()..." << std::endl;
    }

};

// Unit tests for the profile class
TEST_F(ProfileTest, test_construct_double_profile){

    // Construct profiles with 100 bins spaced equally 1 m apart from 1m above reference using the base constructor
    std::vector<double> z;
    for (int i=0; i<100; i++) z.push_back(double(i+1));
    Bins bins = Bins(z);

    Profile<double> p1(bins);
    Profile<double> p2(bins, 10.1, 0.2, 0.1);

    // Ensure that the printing operator does not error (before values are assigned)
    std::cout << p1 << std::endl;

    // Assign values (init as zeros, the value doesn't matter here)
    std::vector<double> vals(100);
    p1.setValues(vals);

    // This should crap out due to incorrect number of bins
    vals.push_back(0.0);
    try {
        p1.setValues(vals);
        FAIL() << "Expected std::out_of_range exception";
    }
    catch(std::out_of_range const & err) {
        EXPECT_EQ(err.what(), std::string("size of vector 'values' does not equal the number of bins for this profile"));
    }
    catch(...) {
        FAIL() << "Expected std::out_of_range - exception of another type found";
    }

}

// Test fixture for the fitting process
class FitTest : public ::testing::Test {

protected:

    // TODO remove any or all of the following functions if bodies are empty

    FitTest() {
        // You can do set-up work for each test here.
    }

    virtual ~FitTest() {
        // You can do clean-up work that doesn't throw exceptions here.
    }

    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    virtual void SetUp() {
        // Code here will be called immediately after the constructor (right
        // before each test).
        std::cout << std::endl << "Setting up FitTest()..." << std::endl;
    }

    virtual void TearDown() {
        // Code here will be called immediately after each test (right
        // before the destructor).
        std::cout << std::endl << "Tearing down FitTest()..." << std::endl;
    }

    // Objects declared here can be used by all tests
    int a= 11;
};
