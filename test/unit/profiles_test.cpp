/*
 * PROFILES_TEST.CPP Test fixtures for analytical profiles
 *
 * References:
 *
 *   [1]
 *
 * Future Improvements:
 *
 *   [1]
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

// Test fixture for generating analytical profiles
class AnalyticalProfileTest : public ::testing::Test {

protected:

    virtual void SetUp() {
        std::cout << std::endl << "Setting up AnalyticalProfileTest()..." << std::endl;
    }

    virtual void TearDown() {
        std::cout << "Tearing down AnalyticalProfileTest()..." << std::endl << std::endl;
    }

};

// Unit tests for the profile class
TEST_F(AnalyticalProfileTest, test_adem_profile){

    /*
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
*/

}
