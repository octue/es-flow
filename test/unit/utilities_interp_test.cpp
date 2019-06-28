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

#include "utilities/interp.h"


using namespace utilities;


class InterpTest : public ::testing::Test {};


TEST_F(InterpTest, test_linear_interp_double) {

    Eigen::VectorXd x(5);
    Eigen::VectorXd y(5);

    x << 0, 1, 2, 3, 4;
    y << 0, 2, 4, 6, 8;
    double xi = 1.5;

    utilities::LinearInterpolant s(x, y);
    double yi = s(xi);

    ASSERT_NEAR(yi, 3.0, 0.00001);

}

TEST_F(InterpTest, test_linear_interp_at_nodes) {

    Eigen::VectorXd x(5);
    Eigen::VectorXd y(5);
    x << 0, 1, 2, 3, 4;
    y << 0, 2, 4, 6, 8;
    utilities::LinearInterpolant s(x, y);

    // First node (endpoint)
    double xi = 0.0;
    double yi = s(xi);
    ASSERT_NEAR(yi, 0.0, 0.00001);

    // Any inter node
    xi = 3.0;
    yi = s(xi);
    ASSERT_NEAR(yi, 6.0, 0.00001);

    // Last node (endpoint)
    xi = 4.0;
    yi = s(xi);
    ASSERT_NEAR(yi, 8.0, 0.00001);

}


TEST_F(InterpTest, test_linear_interp_beyond_endpoints) {

    Eigen::VectorXd x(5);
    Eigen::VectorXd y(5);
    x << 0, 1, 2, 3, 4;
    y << 0, 2, 4, 6, 8;
    utilities::LinearInterpolant s(x, y);

    double xi = -1.0;
    double yi = s(xi);
    ASSERT_NEAR(yi, 0.0, 0.00001);

    xi = 5.0;
    yi = s(xi);
    ASSERT_NEAR(yi, 8.0, 0.00001);

}


TEST_F(InterpTest, test_linear_interp_vector) {

    Eigen::ArrayXd x(5);
    Eigen::ArrayXd xi(2);
    Eigen::ArrayXd y(5);
    Eigen::ArrayXd yi(2);
    Eigen::ArrayXd yi_correct(2);

    x << 0, 1, 2, 3, 4;
    y << 0, 2, 4, 6, 8;
    xi << 0.5, 1.5;
    yi_correct << 1.0, 3.0;

    utilities::LinearInterpolant s(x, y);
    yi = s(xi);

    ASSERT_TRUE(yi.isApprox(yi_correct));

}


TEST_F(InterpTest, test_cubic_interp_double) {

    Eigen::VectorXd x(5);
    Eigen::VectorXd y(5);

    x << 0, 1, 2, 3, 4;
    y << 0, 1, 4, 9, 16;
    double xi = 1.5;

    utilities::CubicSplineInterpolant s(x, y);
    double yi = s(xi);

    ASSERT_NEAR(yi, 2.25, 0.00001);

}


TEST_F(InterpTest, test_cubic_interp_vector) {

    Eigen::VectorXd x(5);
    Eigen::VectorXd xi(2);
    Eigen::VectorXd y(5);
    Eigen::VectorXd yi(2);
    Eigen::VectorXd yi_correct(2);

    x << 0, 1, 2, 3, 4;
    y << 0, 1, 4, 9, 16;
    xi << 0.5, 1.5;
    yi_correct << 0.25, 2.25;

    utilities::CubicSplineInterpolant s(x, y);
    yi = s(xi);

    ASSERT_TRUE(yi.isApprox(yi_correct));

}
