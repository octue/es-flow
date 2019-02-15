/*
 * utilities_cumtrapz_test.cpp Test fixtures for the cumtrapz function
 *
 * Author:                   Tom Clark (thclark @ github)
 *
 * Copyright (c) 2016-9 Octue Ltd. All Rights Reserved.
 *
 */

#include "gtest/gtest.h"
#include <Eigen/Dense>

#include "utilities/cumtrapz.h"


using namespace es;


class CumtrapzTest : public ::testing::Test {};


TEST_F(CumtrapzTest, test_colwise_cumtrapz_1_row) {

    // Test that zeros are returned when trying to colwise integrate an array with 1 row
    Eigen::ArrayXXf integrand(1, 4);
    integrand << 1, 2, 6, 9;
    Eigen::ArrayXXf integral_correct(1, 4);
    integral_correct << 0, 0, 0, 0;
    Eigen::ArrayXXf integral = cumtrapz(integrand);
    EXPECT_EQ(integral.matrix(), integral_correct.matrix());

}


TEST_F(CumtrapzTest, test_colwise_cumtrapz_2_rows) {

    // Test that uniform spaced integration works on an array with 2 rows
    Eigen::ArrayXXf integrand(2,6);
    integrand << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
        7.0, 8.0, 9.0, 10.0, 11.0, 12.0;
    Eigen::ArrayXXf integral_correct(2,6);
    integral_correct << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;
    Eigen::ArrayXXf integral = cumtrapz(integrand);
    EXPECT_EQ(integral.matrix(), integral_correct.matrix());

}


TEST_F(CumtrapzTest, test_colwise_cumtrapz_6_rows) {

    Eigen::ArrayXXf integrand(6,1);
    integrand << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;

    Eigen::ArrayXXf integral_correct(6,1);
    integral_correct<< 0.0, 1.5, 4.0, 7.5, 12.0, 17.5;

    Eigen::ArrayXXf integral = cumtrapz(integrand);
    EXPECT_EQ(integral.matrix(), integral_correct.matrix());

}


TEST_F(CumtrapzTest, test_colwise_nonuniform_cumtrapz) {

    // Test that non-uniform spacing works with an n x 1 spacing array
    Eigen::ArrayXXf spacing(3, 1);
    spacing << 2,
        5,
        7;
    Eigen::ArrayXXf integrand(3,4);
    integrand << 1, 2, 6, 9,
        3, 1, 7, 2,
        4, 8, 3, 1;
    Eigen::ArrayXXf integral_correct(3,4);
    integral_correct << 0.0, 0.0, 0.0, 0.0,
                        6.0, 4.5, 19.5, 16.5,
                        13.0, 13.5, 29.5, 19.5;
    Eigen::ArrayXXf integral(3,4);
    cumtrapz(integral, spacing, integrand);
    EXPECT_EQ(integral.matrix(), integral_correct.matrix());

}
