/*
 * PROFILES_TEST.CPP Test fixtures for analytical profiles
 *
 * References:
 *
 *   [1]
 *
 * Future Improvements:
 *
 *   [1]
 *
 * Author:                   T. Clark
 * Work address:             Ocean Array Systems Ltd
 *                           Hauser Forum
 *                           3 Charles Babbage Road
 *                           Cambridge
 *                           CB3 0GT
 * Email:                    tom.clark@oceanarraysystems.com
 * Website:                  www.oceanarraysystems.com
 *
 * Copyright (c) 2016-17 Ocean Array Systems, All Rights Reserved.
 *
 */

#include "gtest/gtest.h"
#include "profile.h"
#include "relations/velocity.h"
#include "relations/stress.h"
#include "constants.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <unsupported/Eigen/AutoDiff>

using namespace es;
using namespace Eigen;

// Test fixture for generating analytical profiles
class AnalyticalProfileTest : public ::testing::Test {

protected:

    virtual void SetUp() {
        std::cout << std::endl << "Setting up AnalyticalProfileTest()..." << std::endl;
    }

    virtual void TearDown() {
        std::cout << "Tearing down AnalyticalProfileTest()..." << std::endl << std::endl;
    }

};

TEST_F(AnalyticalProfileTest, test_power_law_profile) {
    // Get analytical values for velocity using power law profile
    std::cout << "test_power_law_profile" << std::endl;

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

//TEST_F(AnalyticalProfileTest, test_most_profile) {
//    // Get analytical values for velocity using log law profile and psi function
//    std::cout << "test_most_profile" << std::endl;
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

TEST_F(AnalyticalProfileTest, test_marusic_jones_profile) {
    // Get analytical values for velocity using law of wall and wake
    std::cout << "test_marusic_jones_profile" << std::endl;

    // Basic test parameters
    double pi_j = 0.42;
    double kappa = 0.41;
    double delta = 1000.0;
    double u_inf = 20.0;
    double s = 23.6;
    double u_tau = u_inf / s;
    double z_0 = 0.0;

    // Check that it works for a z value of type double
    double z_doub = 1.;
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









    // Print variables to plot comparison with MATLAB based equivalent calculation
    //    std::cout << "pi_j = " << pi_j << ";" << std::endl;
    //    std::cout << "kappa = " << kappa << ";" << std::endl;
    //    std::cout << "delta = " << delta << ";" << std::endl;
    //    std::cout << "s = " << s << ";" << std::endl;
    //    std::cout << "u_inf = " << u_inf << ";" << std::endl;
    //    std::cout << "z_0 = " << z_0 << ";" << std::endl;
    //    std::cout << "z = [" << z << "];" << std::endl;
}

TEST_F(AnalyticalProfileTest, test_r13_profile) {
    // Get an R13 profile from a basic parameter set

    // Test a basic integrator out
    double arg = 0.0;
    double res = do_something(arg);

    std::cout << "test_r13_profile" << std::endl;
}

TEST_F(AnalyticalProfileTest, test_lewkowicz_profile) {
    // Get analytical values for velocity using lewkowicz law of wall and wake
    std::cout << "test_lewkowicz_profile" << std::endl;

    // Basic test parameters
    double pi_coles = 0.42;
    double kappa = 0.41;
    double delta = 1000.0;
    double u_inf = 20.0;
    double s = 23.6;
    double u_tau = u_inf / s;
    double z_0 = 0.0;

    // Check that it works for a z value of type double
    double eta_doub = 1.0/delta;
    double speed1 = lewkowicz_speed(eta_doub, pi_coles, kappa, u_inf, u_tau);
    std::cout << "checked scalar double operation (U = " << speed1 << " m/s)" << std::endl;

    // Check that it works for a VectorXd input (vertically spaced z)
    double low = 1;
    double high = 100;
    size_t n_bins = 100;
    VectorXd z = VectorXd::LinSpaced(n_bins, low, high);
    VectorXd eta = z.array() / delta;
    VectorXd speed = lewkowicz_speed(eta, pi_coles, kappa, u_inf, u_tau);
    std::cout << "checked VectorXd operation" << std::endl;

    // Check that it works for an AutoDiffScalar
    // Also provides minimal example of how to get the derivative through the profile
    typedef Eigen::AutoDiffScalar<Eigen::VectorXd> ADScalar;
    ADScalar ads_eta;
    ADScalar ads_speed;
    VectorXd dspeed_dz;
    dspeed_dz.setZero(n_bins);
    for (int k = 0; k < n_bins; k++) {
        ads_eta.value() = eta[k];
        ads_eta.derivatives() = Eigen::VectorXd::Unit(1, 0);  // Also works once outside the loop without resetting the derivative guess each step
        ads_speed = lewkowicz_speed(ads_eta, pi_coles, kappa, u_inf, u_tau);
        dspeed_dz[k] = ads_speed.derivatives()[0] / delta;
    }
    std::cout << "checked AutoDiffScalar operation" << std::endl;

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


TEST_F(AnalyticalProfileTest, test_veer_profile) {
    std::cout << "test_veer_profile" << std::endl;

//    // Elevation of site in degrees latitude
//    double phi_latitude = 52;
//
//    Vector3d k;
//    k << 0, 0, 2*omega_world*sind(phi_latitude);
//
//    VectorXd v_bar;
//    v_bar << 10., 3.;


    // Assume homogeneous, isotropic turbulence such that R13 = R23. Not valid, but OK for the sake of the unit test.



//        size_t dim_x = 28, dim_y = 126;
//    Eigen::FFT<float> fft;
//    Eigen::MatrixXf in = Eigen::MatrixXf::Random(dim_x, dim_y);
//    Eigen::MatrixXcf out;
//    out.setZero(dim_x, dim_y);
//
//    for (int k = 0; k < in.rows(); k++) {
//        Eigen::VectorXcf tmpOut(dim_x);
//        fft.fwd(tmpOut, in.row(k));
//        out.row(k) = tmpOut;
//    }
//
//    for (int k = 0; k < in.cols(); k++) {
//        Eigen::VectorXcf tmpOut(dim_y);
//        fft.fwd(tmpOut, out.col(k));
//        out.col(k) = tmpOut;
//    }

    /*

    Profile<double> p1(bins);
    Profile<double> p2(bins, 10.1, 0.2, 0.1);

// Ensure that the printing operator does not error (before values are assigned)
    std::cout << p1 << std::endl;

// Assign values (init as zeros, the value doesn't matter here)
    std::vector<double> vals(100);
    p1.setValues(vals);

// This should crap out due to incorrect number of bins
    vals.push_back(0.0);
    try {
        p1.setValues(vals);
        FAIL() << "Expected std::out_of_range exception";
    }
    catch(std::out_of_range const & err) {
        EXPECT_EQ(err.what(), std::string("size of vector 'values' does not equal the number of bins for this profile"));
    }
    catch(...) {
        FAIL() << "Expected std::out_of_range - exception of another type found";
    }
*/

}
