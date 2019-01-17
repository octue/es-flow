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

    // Ensemble signatures as we load them
    EddySignature signature_a = EddySignature();
    signature_a.load(data_path + std::string("/signatures_A.mat"));
    EddySignature signature_b = EddySignature();
    signature_b.load(data_path + std::string("/signatures_B1.mat"));
    EddySignature signature_bx = EddySignature();
    signature_bx.load(data_path + std::string("/signatures_B2.mat"));
    signature_b = signature_b + signature_bx;
    signature_bx.load(data_path + std::string("/signatures_B3.mat"));
    signature_b = signature_b + signature_bx;
    signature_bx.load(data_path + std::string("/signatures_B4.mat"));
    signature_b = signature_b + signature_bx;

    // Basic test parameters
    double pi_coles = 0.42;
    double kappa = 0.41;
    double delta_c = 1000.0;
    double u_inf = 20.0;
    double shear_ratio = 23.6;
    double u_tau = u_inf / shear_ratio;
    double z_0 = 0.0;
    double beta = 0.0;
    double zeta = 0.0;

    // Run ADEM for these parameters to get full spectra and stuff
    AdemData data = adem(beta, delta_c, kappa, pi_coles, shear_ratio, u_inf, zeta, signature_a, signature_b);

    // Load in the verification data (computed with MATLAB)
//    AdemData check_data = AdemData();
//    check_data.load(data_path + std::string("/adem_data_check.mat"));



    // Print useful diagnostics values
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

