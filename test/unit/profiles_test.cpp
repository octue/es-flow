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
#include "profiles_velocity.h"
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
    std::cout << "speed:       " << speed.transpose() << std::endl;
    std::cout << "derivative: " << dspeed_dz.transpose() << std::endl;
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
//    std::cout << "speed:       " << speed.transpose() << std::endl;
//    std::cout << "derivative: " << dspeed_dz.transpose() << std::endl;
//}

TEST_F(AnalyticalProfileTest, test_marusic_jones_profile) {
    // Get analytical values for velocity using law of wall and wake
    std::cout << "test_marusic_jones_profile" << std::endl;

    // Vertical spacing
    double low = 1;
    double high = 100;
    size_t n_bins = 100;
    VectorXd z = VectorXd::LinSpaced(n_bins, low, high);

    // Jones' modification of the Coles wake factor
    double pi_j = 0.42;
    // von karman constant
    double kappa = 0.41;
    // boundary layer thickness
    double delta = 1000.0;
    // free stream velocity
    double u_inf = 20.0;
    // shear / skin friction velocity (s = u_inf / u_tau)
    double s = 23.6;
    double u_tau = u_inf / s;
    // distance of hypothetical smooth wall from actual rough wall z0 = 0.25k_s
    double z_0 = 0.0;

    VectorXd eta = (z.array() + z_0) / (delta + z_0);
    VectorXd eta_cubed = eta.array().cube();
    VectorXd term1 = eta.array().log() / kappa;
    VectorXd term2 = (eta_cubed.array() - 1.0) / (3.0 * kappa);
    VectorXd term3 = 2.0 * pi_j * (1.0 - eta.array().square() * 3.0 + eta_cubed.array() * 2.0) / kappa;
    VectorXd u_deficit = term2 - term1 + term3;
    VectorXd u_bar = u_inf - u_deficit.array() * u_tau;

    // Print variables to plot comparison with MATLAB based equivalent calculation
    //    std::cout << "pi_j = " << pi_j << ";" << std::endl;
    //    std::cout << "kappa = " << kappa << ";" << std::endl;
    //    std::cout << "delta = " << delta << ";" << std::endl;
    //    std::cout << "s = " << s << ";" << std::endl;
    //    std::cout << "u_inf = " << u_inf << ";" << std::endl;
    //    std::cout << "z_0 = " << z_0 << ";" << std::endl;
    //    std::cout << "z = [" << z << "];" << std::endl;
    //    std::cout << "u_bar = [" <<u_bar << "];" << std::endl;
}

TEST_F(AnalyticalProfileTest, test_r13_profile) {
    // Get an R13 profile from a basic parameter set
    std::cout << "test_r13_profile" << std::endl;
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
