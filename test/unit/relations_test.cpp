/*
 * VELOCITY_RELATIONS_TEST.CPP Test fixtures for analytical profiles
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

#include "definitions.h"
#include "relations/stress.h"
#include "relations/veer.h"
#include "relations/velocity.h"


using namespace es;
using namespace Eigen;


// Test fixture for generating analytical profiles
class VelocityRelationsTest : public ::testing::Test {};


TEST_F(VelocityRelationsTest, test_power_law_profile) {

    // Basic test parameters
    double z_ref = 60;
    double u_ref = 10;
    double alpha = 0.3;

    // Check that it works for a z value of type double
    double z_doub = 1.;
    double speed2 = power_law_speed(z_doub, u_ref, z_ref, alpha);
    std::cout << "checked scalar double operation (U = " << speed2 << " m/s)" << std::endl;

    // Check that it works for a VectorXd input (vertically spaced z)
    double low = 1;
    double high = 100;
    size_t n_bins = 100;
    VectorXd z = VectorXd::LinSpaced(n_bins, low, high);
    VectorXd speed = power_law_speed(z, u_ref, z_ref, alpha);
    std::cout << "checked VectorXd operation" << std::endl;

    // Check that it works for an AutoDiffScalar
    // Also provides minimal example of how to get the derivative through the profile
    typedef Eigen::AutoDiffScalar<Eigen::VectorXd> ADScalar;
    ADScalar ads_z;
    ADScalar ads_speed;
    VectorXd dspeed_dz;
    dspeed_dz.setZero(n_bins);

    typedef Eigen::AutoDiffScalar<Eigen::VectorXd> ADScalar;
    for (int k = 0; k < n_bins; k++) {
        ads_z.value() = z[k];
        ads_z.derivatives() = Eigen::VectorXd::Unit(1, 0);  // Also works once outside the loop without resetting the
                                                            // derivative guess each step
        ads_speed = power_law_speed(ads_z, u_ref, z_ref, alpha);
        dspeed_dz[k] = ads_speed.derivatives()[0];
    }
    std::cout << "checked AutoDiffScalar operation" << std::endl;

    // Print useful diagnostics values
    std::cout << "speed = ["     << speed.transpose()     << "];" << std::endl;
    std::cout << "dspeed_dz = [" << dspeed_dz.transpose() << "];" << std::endl;
}


//TEST_F(VelocityRelationsTest, test_most_profile) {
//
//    // von karman constant
//    double kappa = 0.41;
//    // zero plane offset distance (e.g. for forest canopies)
//    double d = 0;
//    // roughness length
//    double z0 = 0.5;
//    // Monin-Obukhov length
//    double L = 20.;
//
//    // Check that it works for a z value of type double
//    double z_doub = 1.;
//    double speed2 = most_law_speed(z_doub, kappa, d, z0, L);
//    std::cout << "checked scalar double operation (U = " << speed2 << " m/s)" << std::endl;
//
//    // Check that it works for a VectorXd input (vertically spaced z)
//    double low = 1;
//    double high = 100;
//    size_t n_bins = 100;
//    VectorXd z = VectorXd::LinSpaced(n_bins, low, high);
//    VectorXd speed = most_law_speed(z, kappa, d, z0, L);
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
//        ads_z.value() = z[k];
//        ads_z.derivatives() = Eigen::VectorXd::Unit(1, 0);  // Also works once outside the loop without resetting the
//        // derivative guess each step
//        ads_speed = most_law_speed(ads_z, kappa, d, z0, L);
//        dspeed_dz[k] = ads_speed.derivatives()[0];
//    }
//    std::cout << "checked AutoDiffScalar operation" << std::endl;
//
//    // Print useful diagnostics values
//    std::cout << "speed = ["     << speed.transpose()     << "];" << std::endl;
//    std::cout << "dspeed_dz = [" << dspeed_dz.transpose() << "];" << std::endl;
//}


TEST_F(VelocityRelationsTest, test_marusic_jones_profile) {

    // Basic test parameters
    double pi_j = 0.42;
    double kappa = 0.41;
    double delta = 1000.0;
    double u_inf = 20.0;
    double s = 23.6;
    double u_tau = u_inf / s;
    double z_0 = 0.0;

    // Check that it works for a z value of type double
    double z_doub = 1.0;
    double speed1 = marusic_jones_speed(z_doub, pi_j, kappa, z_0, delta, u_inf, u_tau);
    std::cout << "checked scalar double operation (U = " << speed1 << " m/s)" << std::endl;

    // Check that it works for a VectorXd input (vertically spaced z)
    double low = 1;
    double high = 100;
    size_t n_bins = 100;
    VectorXd z = VectorXd::LinSpaced(n_bins, low, high);
    VectorXd speed = marusic_jones_speed(z, pi_j, kappa, z_0, delta, u_inf, u_tau);
    std::cout << "checked VectorXd operation" << std::endl;

    // Check that it works for an AutoDiffScalar in z
    // Also provides minimal example of how to get the derivative through the profile
    typedef Eigen::AutoDiffScalar<Eigen::VectorXd> ADScalar;
    ADScalar ads_z;
    ADScalar ads_speed;
    VectorXd dspeed_dz;
    dspeed_dz.setZero(n_bins);
    for (int k = 0; k < n_bins; k++) {
        ads_z.value() = z[k];
        ads_z.derivatives() = Eigen::VectorXd::Unit(1, 0);  // Also works once outside the loop without resetting the derivative guess each step
        ads_speed = marusic_jones_speed(ads_z, pi_j, kappa, z_0, delta, u_inf, u_tau);
//        std::cout << ads_speed.value() << std::endl;
        dspeed_dz[k] = ads_speed.derivatives()[0];
    }
    std::cout << "checked AutoDiffScalar operation" << std::endl;

    // Print useful diagnostics values
    std::cout << "speed = ["     << speed.transpose()     << "];" << std::endl;
    std::cout << "dspeed_dz = [" << dspeed_dz.transpose() << "];" << std::endl;

    // Check it works for an AutoDiffScalar in z and pi_j
    VectorXd partial_dspeed_dz, partial_dspeed_dpi, partial_dspeed_dpi_check;
    ADScalar ads_pi_j, ads_speed_check1, ads_speed_check2;
    partial_dspeed_dz.setZero(n_bins);
    partial_dspeed_dpi.setZero(n_bins);
    partial_dspeed_dpi_check.setZero(n_bins);
    for (int k = 0; k < n_bins; k++) {

        // Set the current values
        ads_z.value() = z[k];
        ads_pi_j.value() = pi_j;

        // Initialise both derivative terms to unit vectors, telling us which term relates to which derivative
        ads_z.derivatives() = Eigen::VectorXd::Unit(2,0);
        ads_pi_j.derivatives() = Eigen::VectorXd::Unit(2,1);
        ads_speed = marusic_jones_speed(ads_z, ads_pi_j, kappa, z_0, delta, u_inf, u_tau);
        partial_dspeed_dz[k] = ads_speed.derivatives()[0];
        partial_dspeed_dpi[k] = ads_speed.derivatives()[1];

        // Run a check with a crude central differencer
        ads_speed_check1 = marusic_jones_speed(ads_z, pi_j-0.01, kappa, z_0, delta, u_inf, u_tau);
        ads_speed_check2 = marusic_jones_speed(ads_z, pi_j+0.01, kappa, z_0, delta, u_inf, u_tau);
        ads_speed_check2 = (ads_speed_check2 - ads_speed_check1) / 0.02;
        partial_dspeed_dpi_check[k] = ads_speed_check2.value();

    }
    std::cout << "partial_dspeed_dpi       = [" << partial_dspeed_dpi.transpose() << "];" << std::endl;
    std::cout << "partial_dspeed_dpi_check = [" << partial_dspeed_dpi_check.transpose() << "];" << std::endl;
}


TEST_F(VelocityRelationsTest, test_lewkowicz_profile) {

    // Basic test parameters
    double pi_coles = 0.42;
    double kappa = 0.41;
    double delta = 1000.0;
    double u_inf = 20.0;
    double shear_ratio = 23.6;
    double z_0 = 0.0;

    // Check that it works for a z value of type double
    double z_doub = 1.0;
    double speed1 = lewkowicz_speed(z_doub, pi_coles, kappa, u_inf, shear_ratio, delta);
    std::cout << "checked scalar double operation (U = " << speed1 << " m/s)" << std::endl;

    // Check that it works for a ArrayXd input (vertically spaced z)
    double low = 1;
    double high = 100;
    int n_bins = 10;
    VectorXd z = VectorXd::LinSpaced(n_bins, low, high);
    VectorXd speed = lewkowicz_speed(z, pi_coles, kappa, u_inf, shear_ratio, delta);
    std::cout << "checked VectorXd operation: "<< speed.transpose() << std::endl;

    // Check that it works for an AutoDiffScalar
    // Also provides minimal example of how to get the derivative through the profile
    typedef Eigen::AutoDiffScalar<Eigen::VectorXd> ADScalar;
    ADScalar ads_z;
    ADScalar ads_speed;
    VectorXd dspeed_dz;
    dspeed_dz.setZero(n_bins);
    for (int k = 0; k < n_bins; k++) {
        ads_z.value() = z[k];
        ads_z.derivatives() = Eigen::VectorXd::Unit(1, 0);  // Also works once outside the loop without resetting the derivative guess each step
        ads_speed = lewkowicz_speed(ads_z, pi_coles, kappa, u_inf, shear_ratio, delta);
        dspeed_dz[k] = ads_speed.derivatives()[0];
    }
    std::cout << "checked AutoDiffScalar<VectorXd> operation" << std::endl;

    // Print useful diagnostics values
    std::cout << "speed = ["     << speed.transpose()     << "];" << std::endl;
    std::cout << "dspeed_dz = [" << dspeed_dz.transpose() << "];" << std::endl;

    // Print variables to plot comparison with MATLAB based equivalent calculation
    //    std::cout << "pi_j = " << pi_j << ";" << std::endl;
    //    std::cout << "kappa = " << kappa << ";" << std::endl;
    //    std::cout << "delta = " << delta << ";" << std::endl;
    //    std::cout << "s = " << s << ";" << std::endl;
    //    std::cout << "u_inf = " << u_inf << ";" << std::endl;
    //    std::cout << "z_0 = " << z_0 << ";" << std::endl;
    //    std::cout << "z = [" << z << "];" << std::endl;

}


class VeerRelationsTest : public ::testing::Test {};


TEST_F(VeerRelationsTest, test_veer_profile) {

//    // Elevation of site in degrees latitude
//    double phi_latitude = 52;
//
//    Vector3d k;
//    k << 0, 0, 2*omega_world*sind(phi_latitude);
//
//    VectorXd v_bar;
//    v_bar << 10., 3.;

    // Assume homogeneous, isotropic turbulence such that R13 = R23. Not valid, but OK for the sake of the unit test.

}
