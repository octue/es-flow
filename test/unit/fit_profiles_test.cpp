/*
 * Fit_test.cpp
 *
 *  Created on: 29 Nov 2016
 *      Author: thc29
 */

#include "gtest/gtest.h"
#include "fit_profiles.h"

using namespace es;

TEST(basic_check, test_eq) {
    Fit my_fit;
    int my_check_value = 0;

    EXPECT_EQ(my_check_value, my_fit.Fitness());
}
