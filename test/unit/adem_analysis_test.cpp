/*
 * adem_analysis_test.cpp Test fixtures for running adem analysis
 *
 * Author:                   Tom Clark (thclark @ github)
 *
 * Copyright (c) 2016-9 Octue Ltd. All Rights Reserved.
 *
 */

#include <Eigen/Dense>
#include <Eigen/Core>
#include <unsupported/Eigen/AutoDiff>

#include "gtest/gtest.h"

#include "adem/adem.h"
#include "relations/stress.h"
#include "utilities/filter.h"


using namespace es;
using namespace Eigen;
using namespace utilities;

// Add the test environment
extern ::testing::Environment* const environment;


class AdemAnalysisTest : public ::testing::Test {
public:

    std::string data_path;

protected:

    virtual void SetUp() {

        // Get the test data directory path from the environment variable
        if(const char* env_p = std::getenv("TEST_DATA_DIR")) {
            std::cout << "TEST_DATA_DIR = " << env_p << std::endl;
            data_path = env_p;
        } else {
            throw std::invalid_argument("Invalid environment variable 'TEST_DATA_DIR'");
        }
    }

};


TEST_F(AdemAnalysisTest, test_analysis) {

    // Ensemble signatures as we load them, divide to average the B signatures.
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

    // TODO This divide-by-4 should be in the code somewhere!!!!
    signature_b = signature_b / 4.0;

    // Basic test parameters
    double beta = 0.0;
    double delta_c = 1000.0;
    double kappa = 0.41;
    double pi_coles = 0.42;
    double shear_ratio = 23.6;
    double u_inf = 2.2;
    double u_tau = u_inf / shear_ratio;
    double z_0 = 0.0;
    double zeta = 0.0;

    // Run ADEM for these parameters to get full spectra and stuff
    AdemData data = adem(beta, delta_c, kappa, pi_coles, shear_ratio, u_inf, zeta, signature_a, signature_b);

    // Verification data
    Eigen::Tensor<double, 3> psi;
    psi = data.psi_a + data.psi_b;

    // Verify against data computed with MATLAB
    ASSERT_NEAR(psi(0, 10, 3), 0, 1e-6);
    ASSERT_NEAR(psi(99, 10, 3), 0.0086, 1e-4);
    ASSERT_NEAR(psi(9953, 10, 3),  0.0142, 1e-4);

    ASSERT_NEAR(psi(0, 4, 5),  0, 1e-6);
    ASSERT_NEAR(psi(99, 4, 5), 0.0038, 1e-4);
    ASSERT_NEAR(psi(9953, 4, 5),  0.0031, 1e-4);

    ASSERT_NEAR(psi(0, 48, 2), 0, 1e-6);
    ASSERT_NEAR(psi(99, 48, 2), -1.3117e-04, 1e-05);
    ASSERT_NEAR(psi(9953, 48, 2),  -6.9840e-05, 1e-06);

    // Print useful diagnostics values
    // std::cout << "speed = ["     << speed.transpose()     << "];" << std::endl;
    // std::cout << "dspeed_dz = [" << dspeed_dz.transpose() << "];" << std::endl;

    // Print variables to plot comparison with MATLAB based equivalent calculation
    //    std::cout << "pi_j = " << pi_j << ";" << std::endl;
    //    std::cout << "kappa = " << kappa << ";" << std::endl;
    //    std::cout << "delta = " << delta << ";" << std::endl;
    //    std::cout << "s = " << s << ";" << std::endl;
    //    std::cout << "u_inf = " << u_inf << ";" << std::endl;
    //    std::cout << "z_0 = " << z_0 << ";" << std::endl;
    //    std::cout << "z = [" << z << "];" << std::endl;
}


TEST_F(AdemAnalysisTest, test_get_reynolds_stress_13) {

    // Results obtained and validated against MATLAB implementation:
    /* eta = [0.001, 0.1, 0.3, 0.6, 1.0];
     * pi_coles = 0.42;
     * zeta = 0;
     * beta = 0;
     * shear_ratio = 23.6;
     * [r13_a, r13_b] = getReynoldsStress13(eta, pi_coles, shear_ratio, zeta, beta)
     */
    Eigen::ArrayXd eta;
    Eigen::ArrayXd r13_a;
    Eigen::ArrayXd r13_b;
    Eigen::ArrayXd r13_a_correct;
    Eigen::ArrayXd r13_b_correct;
    double beta = 0.0;
    double kappa = 0.41;
    double pi_coles = 0.42;
    double shear_ratio = 23.6;
    double zeta = 0.0;

    eta.setZero(5);
    r13_a_correct.setZero(5);
    r13_b_correct.setZero(5);

    eta << 0.001, 0.1, 0.3, 0.6, 1.0;
    r13_a_correct << -0.991958214632306, -0.663877109187046, -0.323638466476833, -0.072136229566514, 0;
    r13_b_correct << -0.015941694827107, -0.292636189577533, -0.515677923020367, -0.430112433486726, 0.000000000000000;

    reynolds_stress_13(r13_a, r13_b, beta, eta, kappa, pi_coles, shear_ratio, zeta);

    std::cout << "r13a expected: " << r13_a_correct.transpose() << std::endl;
    std::cout << "r13a returned: " << r13_a.transpose() << std::endl;

    std::cout << "r13b expected: " << r13_b_correct.transpose() << std::endl;
    std::cout << "r13b returned: " << r13_b.transpose() << std::endl;

    ASSERT_TRUE(r13_a.isApprox(r13_a_correct));
    ASSERT_TRUE(r13_b.isApprox(r13_b_correct));

}


TEST_F(AdemAnalysisTest, test_filter_and_deconv) {

    Eigen::ArrayXd a(5);
    Eigen::ArrayXd b(6);
    Eigen::ArrayXd x(10);
    Eigen::ArrayXd y;
    Eigen::ArrayXd y_correct(10);
    Eigen::ArrayXd z;
    Eigen::ArrayXd z_correct(2);

    // Set up arbitrary filter inputs
    a << 1,    0.4,  0.2,  0.5,  0.1;
    b << 0.45, 0.34, 0.65, 0.32, 0.46, 0.98;
    x << 0.814723686393179, 0.905791937075619, 0.126986816293506, 0.913375856139019, 0.632359246225410,
         0.097540404999410, 0.278498218867048, 0.546881519204984, 0.957506835434298, 0.964888535199277;

    // Correct outputs for filter and deconvolution determined with MATLAB
    y_correct << 0.366625658876931, 0.537962161506937, 0.606173725715194, 0.770296239521390, 0.607280310491784,
                 1.354464466281476, 0.698884413086931, 0.220025836201273, 1.049300675330637, 0.920311139464538;
    z_correct << 0.450000000000000, 0.160000000000000;

    filter(y, b, a, x);
    deconv(z, b, a);

    std::cout << "y expected: " << y_correct.transpose() << std::endl;
    std::cout << "y returned: " << y.transpose() << std::endl;

    std::cout << "z expected: " << z_correct.transpose() << std::endl;
    std::cout << "z returned: " << z.transpose() << std::endl;

    ASSERT_TRUE(y.isApprox(y_correct));
    ASSERT_TRUE(z.isApprox(z_correct));

}
