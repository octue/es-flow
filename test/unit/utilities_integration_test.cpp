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

#include "utilities/integration.h"


using namespace utilities;

template<typename T_scalar>
class PowFunctor {
public:
    PowFunctor() {};
    T_scalar operator() (T_scalar eta) const{
        return pow(eta, 2.0);
    };
};

class IntegrationTest : public ::testing::Test {};


TEST_F(IntegrationTest, test_colwise_cumulative_integrate_1_row) {

    // Test that zeros are returned when trying to colwise integrate an array with 1 row
    Eigen::ArrayXd integrand(1);
    integrand << 1;
    Eigen::ArrayXd integral_correct(1);
    integral_correct << 0;
    PowFunctor<double> functor;
    Eigen::ArrayXd integral = cumulative_integrate(integrand, functor);
    EXPECT_EQ(integral.matrix(), integral_correct.matrix());

}

TEST_F(IntegrationTest, test_colwise_cumulative_integrate) {
    Eigen::ArrayXd integrand(3);
    Eigen::ArrayXd integral_correct(3);
    integrand << 0, 2, 4;
    integral_correct << 0, 8.0/3.0, 64.0/3.0;
    PowFunctor<double> pow_functor;
    Eigen::ArrayXd integral = cumulative_integrate(integrand, pow_functor);
    ASSERT_TRUE(integral.isApprox(integral_correct));
}
