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
    double alpha = 0.3;

    // Create a power law distribution with random noise added
    Eigen::ArrayXd z = Eigen::ArrayXd::LinSpaced(20, 1, 20);
    Eigen::ArrayXd u(20);
    u = power_law_speed(z, u_ref, z_ref, alpha) + ArrayXd::Random(20) - 0.5;

    std::cout << "basic power law profile (U = " << u << " m/s)" << std::endl;

//    // Check that it works for a VectorXd input (vertically spaced z)
//    double low = 1;
//    double high = 100;
//    size_t n_bins = 100;
//    VectorXd z = VectorXd::LinSpaced(n_bins, low, high);
//    VectorXd speed = power_law_speed(z, u_ref, z_ref, alpha);
//    std::cout << "checked VectorXd operation" << std::endl;
//
//    // Check that it works for an AutoDiffScalar
//    // Also provides minimal example of how to get the derivative through the profile
//    typedef Eigen::AutoDiffScalar<Eigen::VectorXd> ADScalar;
//    ADScalar ads_z;
//    ADScalar ads_speed;
//    VectorXd dspeed_dz;
//    dspeed_dz.setZero(n_bins);
//    for (int k = 0; k < n_bins; k++) {
//    ads_z.value() = z[k];
//    ads_z.derivatives() = Eigen::VectorXd::Unit(1, 0);  // Also works once outside the loop without resetting the
//    // derivative guess each step
//    ads_speed = power_law_speed(ads_z, u_ref, z_ref, alpha);
//    dspeed_dz[k] = ads_speed.derivatives()[0];
//    }
//    std::cout << "checked AutoDiffScalar operation" << std::endl;
//
//    // Print useful diagnostics values
//    std::cout << "speed = ["     << speed.transpose()     << "];" << std::endl;
//    std::cout << "dspeed_dz = [" << dspeed_dz.transpose() << "];" << std::endl;

}

