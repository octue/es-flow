/*
 * utilities_trapz_test.cpp Test fixtures for the trapz function
 *
 * Author:                   T. Clark
 * Email:                    tom@octue.com
 * Website:                  www.octue.com
 *
 * Copyright (c) 2019 Octue Ltd
 *
 */

#include "gtest/gtest.h"
#include <Eigen/Dense>
#include "utilities/trapz.h"


using namespace es;


class TrapzTest : public ::testing::Test {

protected:

    virtual void SetUp() {
        std::cout << std::endl << "Setting up TrapzTest()..." << std::endl;
    }

    virtual void TearDown() {
        std::cout << "Tearing down TrapzTest()..." << std::endl << std::endl;
    }

};

TEST_F(TrapzTest, test_colwise_trapz_1_row) {

    // Test that zeros are returned when trying to colwise integrate an array with 1 row
    Eigen::ArrayXXf integrand(1,4);
    integrand << 1, 2, 6, 9;
    Eigen::ArrayXXf integral_correct(1,4);
    integral_correct<< 0, 0, 0, 0;
    Eigen::ArrayXXf integral = trapz(integrand);
    EXPECT_EQ(integral.matrix(), integral_correct.matrix());

};

TEST_F(TrapzTest, test_colwise_trapz_2_rows) {

    // Test that uniform spaced integration works on an array with 2 rows
    Eigen::ArrayXXf integrand(2,4);
    integrand << 1, 2, 6, 9,
        3, 1, 7, 2;
    Eigen::ArrayXXf integral_correct(1,4);
    integral_correct << 2, 1.5, 6.5, 5.5;
    Eigen::ArrayXXf integral = trapz(integrand);
    EXPECT_EQ(integral.matrix(), integral_correct.matrix());

};

TEST_F(TrapzTest, test_colwise_trapz_3_rows) {

    // Test that uniform spaced integration works on an array with 3 rows
    Eigen::ArrayXXf integrand(3,4);
    integrand << 1, 2, 6, 9,
                 3, 1, 7, 2,
                 4, 8, 3, 1;
    Eigen::ArrayXXf integral_correct(1,4);
    integral_correct << 5.5, 6, 11.5, 7;
    Eigen::ArrayXXf integral = trapz(integrand);
    EXPECT_EQ(integral.matrix(), integral_correct.matrix());

};

TEST_F(TrapzTest, test_colwise_nonuniform_trapz) {

    // Test that non-uniform spacing works with an n x 1 spacing array
    Eigen::ArrayXXf spacing(3, 1);
    spacing << 2,
               5,
               7;
    Eigen::ArrayXXf integrand(3,4);
    integrand << 1, 2, 6, 9,
                 3, 1, 7, 2,
                 4, 8, 3, 1;
    Eigen::ArrayXXf integral_correct(1,4);
    integral_correct << 13, 13.5, 29.5, 19.5;
    Eigen::ArrayXXf integral(1,4);
    trapz(integral, spacing, integrand);
    EXPECT_EQ(integral.matrix(), integral_correct.matrix());

};
