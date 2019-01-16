/*
 * adem_analysis_test.cpp Test fixtures for running adem analysis
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
 * Email:                    tom@octue.com
 * Website:                  www.octue.com
 *
 * Copyright (c) 2019 Octue Ltd, All Rights Reserved.
 *
 */

#include "gtest/gtest.h"
#include "adem/adem.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <unsupported/Eigen/AutoDiff>

using namespace es;
using namespace Eigen;


class AdemAnalysisTest : public ::testing::Test {
public:

    std::string data_path;

protected:

    virtual void SetUp() {
        std::cout << std::endl << "Setting up AdemAnalysisTest..." << std::endl;

        // Get the test data directory path from the environment variable
        if(const char* env_p = std::getenv("TEST_DATA_DIR")) {
            std::cout << "TEST_DATA_DIR = " << env_p << std::endl;
            data_path = env_p;
        } else {
            throw std::invalid_argument("Invalid environment variable 'TEST_DATA_DIR'");
        }
    }

    virtual void TearDown() {
        std::cout << "Tearing down AdemAnalysisTest..." << std::endl << std::endl;
    }

};


TEST_F(AdemAnalysisTest, test_analysis) {

    // Construct an eddy signature object
    EddySignature sig = EddySignature();

    // Load a signature file
    std::string file_name = data_path + std::string("/signatures_A.mat");
    sig.load(file_name);

//    std::cout << "test_analysis" << std::endl;
//
//    // Basic test parameters
//    double pi_coles = 0.42;
//    double kappa = 0.41;
//    double delta = 1000.0;
//    double u_inf = 20.0;
//    double s = 23.6;
//    double u_tau = u_inf / s;
//    double z_0 = 0.0;
//    adem(delta_c, u_inf, pi_coles, s, beta, zeta)
//    // Check that it works for a z value of type double
//    double eta_doub = 1.0/delta;
//    double speed1 = lewkowicz_speed(eta_doub, pi_coles, kappa, u_inf, u_tau);
//    std::cout << "checked scalar double operation (U = " << speed1 << " m/s)" << std::endl;
//
//    // Check that it works for a VectorXd input (vertically spaced z)
//    double low = 1;
//    double high = 100;
//    size_t n_bins = 100;
//    VectorXd z = VectorXd::LinSpaced(n_bins, low, high);
//    VectorXd eta = z.array() / delta;
//    VectorXd speed = lewkowicz_speed(eta, pi_coles, kappa, u_inf, u_tau);
//    std::cout << "checked VectorXd operation" << std::endl;
//
//    // Check that it works for an AutoDiffScalar
//    // Also provides minimal example of how to get the derivative through the profile
//    typedef Eigen::AutoDiffScalar<Eigen::VectorXd> ADScalar;
//    ADScalar ads_eta;
//    ADScalar ads_speed;
//    VectorXd dspeed_dz;
//    dspeed_dz.setZero(n_bins);
//    for (int k = 0; k < n_bins; k++) {
//    ads_eta.value() = eta[k];
//    ads_eta.derivatives() = Eigen::VectorXd::Unit(1, 0);  // Also works once outside the loop without resetting the derivative guess each step
//    ads_speed = lewkowicz_speed(ads_eta, pi_coles, kappa, u_inf, u_tau);
//    dspeed_dz[k] = ads_speed.derivatives()[0] / delta;
//    }
//    std::cout << "checked AutoDiffScalar operation" << std::endl;
//
//    // Print useful diagnostics values
//    std::cout << "speed = ["     << speed.transpose()     << "];" << std::endl;
//    std::cout << "dspeed_dz = [" << dspeed_dz.transpose() << "];" << std::endl;

    // Print variables to plot comparison with MATLAB based equivalent calculation
    //    std::cout << "pi_j = " << pi_j << ";" << std::endl;
    //    std::cout << "kappa = " << kappa << ";" << std::endl;
    //    std::cout << "delta = " << delta << ";" << std::endl;
    //    std::cout << "s = " << s << ";" << std::endl;
    //    std::cout << "u_inf = " << u_inf << ";" << std::endl;
    //    std::cout << "z_0 = " << z_0 << ";" << std::endl;
    //    std::cout << "z = [" << z << "];" << std::endl;
}

