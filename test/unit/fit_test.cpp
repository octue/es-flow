/*
 * fit_test.cpp Test fixtures for fitting velocity profiles
 *
 * Author:                   Tom Clark (thclark @ github)
 *
 * Copyright (c) 2016-9 Octue Ltd. All Rights Reserved.
 *
 */

#include <Eigen/Dense>
#include <Eigen/Core>
#include "gtest/gtest.h"
#include <unsupported/Eigen/AutoDiff>

#include "fit.h"
#include "relations/velocity.h"

#include "cpplot.h"

using namespace es;
using namespace Eigen;


// Add the test environment
extern ::testing::Environment* const environment;


// Test fixture for the fitting process
class FitTest : public ::testing::Test {
};

// Unit tests for fitting routines
TEST_F(FitTest, test_fit){

    // Basic profile parameters
    double z_ref = 60;
    double u_ref = 10;
    double alpha_true = 0.3;

    // Create a power law distribution with random noise added
    Eigen::ArrayXd z = Eigen::ArrayXd::LinSpaced(20, 1, 200);
    Eigen::ArrayXd u_original(20);
    Eigen::ArrayXd u_noisy(20);
    Eigen::ArrayXd u_fitted(20);
    u_original = power_law_speed(z, u_ref, z_ref, alpha_true);
    u_noisy = u_original + ArrayXd::Random(20);

    // Fit to find the value of alpha
    double alpha_fitted = fit_power_law_speed(z, u_noisy, z_ref, u_ref);
    u_fitted = power_law_speed(z, u_ref, z_ref, alpha_fitted);

    // Display original, noisy and fitted profiles on scatter plot
    cpplot::Figure fig = cpplot::Figure();
    cpplot::ScatterPlot p_original = cpplot::ScatterPlot();
    cpplot::ScatterPlot p_noisy = cpplot::ScatterPlot();
    cpplot::ScatterPlot p_fitted = cpplot::ScatterPlot();

    p_original.y = z.matrix();
    p_noisy.y = z.matrix();
    p_fitted.y = z.matrix();

    p_original.x = u_original;
    p_noisy.x = u_noisy;
    p_fitted.x = u_fitted;

    fig.add(p_original);
    fig.add(p_noisy);
    fig.add(p_fitted);

    fig.write("test_fit_power_law_speed.json");

    // Check that the fit worked
    // TODO use a known seed for the random noise which is added to avoid occasional failure here
    ASSERT_NEAR(alpha_fitted,  alpha_true, 0.15);
}

