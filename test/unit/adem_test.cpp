/*
 * adem_test.cpp Test fixtures for running adem analysis
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

#include "cpplot.h"

#include "adem/adem.h"
#include "utilities/filter.h"
#include "utilities/tensors.h"


using namespace es;
using namespace Eigen;
using namespace utilities;
using namespace cpplot;

// Add the test environment
extern ::testing::Environment* const environment;


class AdemTest : public ::testing::Test {
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


TEST_F(AdemTest, test_analysis) {

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
    double zeta = 0.0;

    // Run ADEM for these parameters to get full spectra
    AdemData data = adem(beta, delta_c, kappa, pi_coles, shear_ratio, u_inf, zeta, signature_a, signature_b);

    // Print adem data for reference
    std::cout << data << std::endl;

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

    // Initialise dimensions and spectral tensor arrays
    Eigen::array<Eigen::Index, 2> dims = {psi.dimensions()[0], psi.dimensions()[1]};
    auto n_elems_per_slice = dims[0] * dims[1];
    std::string term_str;

    Eigen::ArrayXXd s(dims[0], dims[1]);
    Eigen::ArrayXXd s11(dims[0], dims[1]);
    Eigen::ArrayXXd s12(dims[0], dims[1]);
    Eigen::ArrayXXd s13(dims[0], dims[1]);
    Eigen::ArrayXXd s22(dims[0], dims[1]);
    Eigen::ArrayXXd s23(dims[0], dims[1]);
    Eigen::ArrayXXd s33(dims[0], dims[1]);

    // For each of the 6 auto / cross spectra terms, map a slice out of the psi tensor
    for (Eigen::Index j = 0; j < 6; j++) {

        std::cout << "j: " << j << std::endl;
        Eigen::Tensor<double, 2> slice_tens = psi.chip(j, 2);
        Eigen::ArrayXXd slice = tensor_to_array(slice_tens, dims[0], dims[1]);

        // Premultiply the spectrum with the wavenumber, which gives the spectral density
        slice = slice;

        // TODO use std vec of pointers and labels duh
        switch (j) {
            case 0: {
                s11 << slice;
                term_str = "s11";
                break;
            }
            case 1: {
                s12 = slice;
                term_str = "s12";
                break;
            }
            case 2: {
                s13 = slice;
                term_str = "s13";
                break;
            }
            case 3: {
                s22 = slice;
                term_str = "s22";
                break;
            }
            case 4: {
                s23 = slice;
                term_str = "s23";
                break;
            }
            case 5: {
                s33 = slice;
                term_str = "s33";
                break;
            }
            default: {}
        };

        // Create a plot to show the t2w terms
        Figure figt = Figure();
        ScatterPlot pt = ScatterPlot();
        auto n = data.t2wa.rows();
        pt.x = Eigen::ArrayXd::LinSpaced(n, 1, n);
        pt.y = data.t2wa;
        pt.name = "t2wa";
        figt.add(pt);
        ScatterPlot pt2 = ScatterPlot();
        pt2.x = Eigen::ArrayXd::LinSpaced(n, 1, n);
        pt2.y = data.t2wb;
        pt2.name = "t2wb";
        figt.add(pt2);
        figt.write("test_t2w_plot.json");

        // Create a plot to show the z distribution
        Figure figz = Figure();
        ScatterPlot pz = ScatterPlot();
        n = data.z.rows();
        pz.x = Eigen::ArrayXd::LinSpaced(n, 1, n);
        pz.y = data.z;
        pz.name = "z";
        figz.add(pz);
        ScatterPlot pz2 = ScatterPlot();
        pz2.x = Eigen::ArrayXd::LinSpaced(n, 1, n);
        pz2.y = data.eta;
        pz2.name = "eta";
        figz.add(pz2);
        ScatterPlot pz3 = ScatterPlot();
        pz3.x = Eigen::ArrayXd::LinSpaced(n, 1, n);
        pz3.y = data.lambda_e;
        pz3.name = "lambda_e";
        figz.add(pz3);
        figz.write("test_z_plot.json");


        // Create a surface plot to show the spectrum term
        Figure fig = Figure();
        SurfacePlot p = SurfacePlot();

        // x direction is frequency
//        p.x = Eigen::RowVectorXd::LinSpaced(s13.cols(), 1, s13.cols()).replicate(s13.rows(), 1).array();
        p.x = data.k1z.transpose().leftCols(400);
        std::cout << "size x " << p.x.rows() << " " << p.x.cols() << std::endl;

        // y direction is vertical height z
        p.y = data.z.transpose().replicate(p.x.rows(), 1).array().leftCols(400);
        std::cout << "size y " << p.y.rows() << " " << p.y.cols() << std::endl;

        // Copy spectrum amplitude from mapped tensor
        p.z = slice.transpose().leftCols(400);
        std::cout << "size z " << p.z.rows() << " " << p.z.cols() << std::endl;

        fig.add(p);
        fig.write("test_spectrum_surface_plot_" + term_str + ".json");
    }
}


TEST_F(AdemTest, test_get_reynolds_stress_13) {

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


TEST_F(AdemTest, test_filter_and_deconv) {

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


TEST_F(AdemTest, test_signature_k1z) {

    // Load a signature
    EddySignature signature_a = EddySignature();
    signature_a.load(data_path + std::string("/signatures_A.mat"));

    //
    ArrayXd eta = ArrayXd::LinSpaced(9, 0.1, 1);
    ArrayXXd k1z = signature_a.k1z(eta);

}
