/*
 * utilities_trapz_test.cpp Test fixtures for the trapz function
 *
 * Author:                   Tom Clark (thclark @ github)
 *
 * Copyright (c) 2016-9 Octue Ltd. All Rights Reserved.
 *
 */

#include "gtest/gtest.h"
#include <Eigen/Dense>

#include "utilities/cumsum.h"


using namespace utilities;


class CumsumTest : public ::testing::Test {};


TEST_F(CumsumTest, test_cumsum) {

    // Test that zeros are returned when trying to colwise integrate an array with 1 row
    Eigen::VectorXd in(4);
    in << 1, 2, 6, 9;
    Eigen::VectorXd out_correct(4);
    out_correct << 1, 3, 9, 18;

    Eigen::VectorXd out;
    cumsum(out, in);
    EXPECT_EQ(out.matrix(), out_correct.matrix());

}
