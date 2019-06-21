/*
 * utilities_convolution_matrix_test.cpp Test fixtures for the convolution_matrix function
 *
 * Author:                   Tom Clark (thclark @ github)
 *
 * Copyright (c) 2016-9 Octue Ltd. All Rights Reserved.
 *
 */

#include "gtest/gtest.h"
#include <Eigen/Dense>

#include "utilities/conv.h"


using namespace utilities;


class ConvTest : public ::testing::Test {};


TEST_F(ConvTest, test_loaded_diagonal_deconv) {

    Eigen::VectorXd input(9);
    Eigen::VectorXd kernel(5);
    input << 1, 1, 2, 3, 2, 1, 2, 3, 2;
    kernel << 1, 2, 3, 2, 1;
    double alpha = 0.01;
    Eigen::VectorXd out = diagonal_loading_deconv(input, kernel, alpha);
    Eigen::VectorXd correct(9);
    correct << 0.6364,
              -0.1880,
               0.6558,
               0.7658,
              -1.7118,
               1.0871,
               2.6934,
              -2.9723,
              -0.5919;
    ASSERT_TRUE(out.isApprox(correct, 0.01));
}


TEST_F(ConvTest, test_convolution_matrix_slim) {
    Eigen::VectorXd k(5);
    Eigen::Index n = 3;
    k << 1, 2, 3, 4, 5;
    Eigen::MatrixXd c = utilities::convolution_matrix(k, n);
    Eigen::MatrixXd correct(7,3);
    correct << 1, 0, 0,
        2, 1, 0,
        3, 2, 1,
        4, 3, 2,
        5, 4, 3,
        0, 5, 4,
        0, 0, 5;
    EXPECT_EQ(c, correct);
}


TEST_F(ConvTest, test_convolution_matrix_fat) {
    Eigen::VectorXd k(3);
    Eigen::Index n = 5;
    k << 1, 2, 3;
    Eigen::MatrixXd c = utilities::convolution_matrix(k, n);
    Eigen::MatrixXd correct(7,5);
    correct << 1, 0, 0, 0, 0,
        2, 1, 0, 0, 0,
        3, 2, 1, 0, 0,
        0, 3, 2, 1, 0,
        0, 0, 3, 2, 1,
        0, 0, 0, 3, 2,
        0, 0, 0, 0, 3;
    EXPECT_EQ(c, correct);
}


TEST_F(ConvTest, test_convolution_matrix_square) {
    Eigen::VectorXd k(5);
    Eigen::Index n = 5;
    k << 1, 2, 3, 4, 5;
    Eigen::MatrixXd c = utilities::convolution_matrix(k, n);
    Eigen::MatrixXd correct(9,5);
    correct << 1, 0, 0, 0, 0,
        2, 1, 0, 0, 0,
        3, 2, 1, 0, 0,
        4, 3, 2, 1, 0,
        5, 4, 3, 2, 1,
        0, 5, 4, 3, 2,
        0, 0, 5, 4, 3,
        0, 0, 0, 5, 4,
        0, 0, 0, 0, 5;
    EXPECT_EQ(c, correct);
}


TEST_F(ConvTest, test_convolution_matrix_one_square) {
    Eigen::VectorXd k(1);
    Eigen::Index n = 1;
    k << 1;
    Eigen::MatrixXd c = utilities::convolution_matrix(k, n);
    Eigen::MatrixXd correct(1,1);
    correct << 1;
    EXPECT_EQ(c, correct);
}
