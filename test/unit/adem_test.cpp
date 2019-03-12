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


// Helper to plot eddy intensities for different eddy types
void plotEddyIntensities(const EddySignature &sig, const std::string &type) {

    // Plot the intensity functions as figures
    Figure fig_i = Figure();
    Figure fig_j = Figure();

    // Plot the J and I functions
    std::vector<std::string> labels = {"11", "12", "13", "22", "23", "33"};
    for (auto col=0; col<6; col++) {
        ScatterPlot p_i = ScatterPlot();
        ScatterPlot p_j = ScatterPlot();
        p_i.x = sig.eta;
        p_i.y = sig.j.col(col);
        p_i.name = "$I_{"+labels[col]+"}$";
        fig_i.add(p_i);
        p_j.x = sig.lambda;
        p_j.y = sig.j.col(col);
        p_j.name = "$J_{"+labels[col]+"}$";
        fig_j.add(p_j);
    }

    // Add axis labeling
    Layout lay_i = Layout("Type " + type + " Eddy Intensity");
    lay_i.xTitle("$\\eta$");
    lay_i.yTitle("$I_{ij}$");
    fig_i.setLayout(lay_i);
    Layout lay_j = Layout("Type " + type + " Eddy Intensity");
    lay_j.xTitle("$\\lambda$");
    lay_j.yTitle("$J_{ij}$");
    fig_j.setLayout(lay_j);

    // Write figures
    fig_i.write("check_type_" + type + "_iij_from_cpp.json");
    fig_j.write("check_type_" + type + "_jij_from_cpp.json");
};


TEST_F(AdemTest, test_get_type_a_eddy_intensity) {
    int n_lambda = 200;
    double dx = 0.01;
    EddySignature sig = EddySignature();
    sig.getEddyIntensity("A", n_lambda, dx);
    plotEddyIntensities(sig, "A");
};


TEST_F(AdemTest, test_get_type_b1_eddy_intensity) {
    int n_lambda = 200;
    double dx = 0.01;
    EddySignature sig = EddySignature();
    sig.getEddyIntensity("B1", n_lambda, dx);
    plotEddyIntensities(sig, "B1");
};


TEST_F(AdemTest, test_get_type_b2_eddy_intensity) {
    int n_lambda = 200;
    double dx = 0.01;
    EddySignature sig = EddySignature();
    sig.getEddyIntensity("B2", n_lambda, dx);
    plotEddyIntensities(sig, "B2");
};


TEST_F(AdemTest, test_get_type_b3_eddy_intensity) {
    int n_lambda = 200;
    double dx = 0.01;
    EddySignature sig = EddySignature();
    sig.getEddyIntensity("B3", n_lambda, dx);
    plotEddyIntensities(sig, "B3");
};


TEST_F(AdemTest, test_get_type_b4_eddy_intensity) {
    int n_lambda = 200;
    double dx = 0.01;
    EddySignature sig = EddySignature();
    sig.getEddyIntensity("B4", n_lambda, dx);
    plotEddyIntensities(sig, "B4");
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
    double zeta = 0.15;

    // Run ADEM for these parameters to get full spectra
    AdemData data = adem(beta, delta_c, kappa, pi_coles, shear_ratio, u_inf, zeta, signature_a, signature_b);

    // Print adem data for reference
    std::cout << data << std::endl;

    // Verification data
    Eigen::Tensor<double, 3> psi;
    psi = data.psi_a + data.psi_b;


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

    // Create a plot to show Reynolds Stress variations with eta
    Figure fig_rs = Figure();
    std::cout << "sizea " << data.reynolds_stress_a.rows() << " " << data.reynolds_stress_a.cols() << std::endl;
    std::cout << "sizeb " << data.reynolds_stress_b.rows() << " " << data.reynolds_stress_b.cols() << std::endl;
    std::vector<std::string> term_labels = { "11", "12", "13", "22", "23", "33" };
    for (Eigen::Index k = 0; k < 6; k++) {
        ScatterPlot p_rs = ScatterPlot();
        p_rs.x = data.eta;
        p_rs.y = data.reynolds_stress_a.col(k);
        p_rs.name = "$R_{" + term_labels[k] + "a}$";
        fig_rs.add(p_rs);
    }
    for (Eigen::Index k = 0; k < 6; k++) {
        ScatterPlot p_rs = ScatterPlot();
        p_rs.x = data.eta;
        p_rs.y = data.reynolds_stress_b.col(k);
        p_rs.name = "$R_{" + term_labels[k] + "b}$";
        fig_rs.add(p_rs);
    }
    Layout lay_rs = Layout();
    lay_rs.xTitle("$\\eta$");
    lay_rs.yTitle("$R_{ij}$");
    fig_rs.setLayout(lay_rs);
    fig_rs.write("test_rs_plot.json");


    // Initialise dimensions and spectral tensor arrays
    Eigen::array<Eigen::Index, 2> dims = {psi.dimensions()[0], psi.dimensions()[1]};
    std::vector<Eigen::ArrayXXd> spectra;
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
        slice = slice.transpose().eval();// * data.k1z;

        // Add to the vector of arrays (much easier to deal with than tensors for later plotting)
        spectra.push_back(slice);

        // Create a surface plot to show the spectrum term
        Figure fig = Figure();
        Layout lay = Layout();
        lay.xTitle("k1z");
        lay.yTitle("z");
        lay.zTitle("Psi");
        lay.xLog();
        fig.setLayout(lay);
        SurfacePlot p = SurfacePlot();

        // x direction is frequency
        p.x = data.k1z.rightCols(600);
        std::cout << "size x " << p.x.rows() << " " << p.x.cols() << std::endl;

        // y direction is vertical height z
        p.y = data.z.transpose().replicate(p.x.rows(), 1).array().rightCols(600);
        std::cout << "size y " << p.y.rows() << " " << p.y.cols() << std::endl;

        // Copy spectrum amplitude from mapped tensor
        p.z = slice.rightCols(600);
        std::cout << "size z " << p.z.rows() << " " << p.z.cols() << std::endl;

        // Trim the zero wavenumber value
        p.x = p.x.topRows(150);
        p.y = p.y.topRows(150);
        p.z = p.z.topRows(150);
        fig.add(p);
        fig.write("test_spectrum_surface_plot_" + term_labels[j] + ".json");
    }

    // Create a plot to show spectra at particular heights. Use nearest above, instead of interpolating
    std::vector<double> values = {0.10, 0.17, 0.27, 0.39, 0.54, 0.72, 0.93};
    std::vector<Eigen::Index> inds;
    for (auto &value : values) {
        for (Eigen::Index i=0; i<data.eta.size(); ++i) {
            if (data.eta(i) >= value) {
                std::cout << "Pushing back index " << i << std::endl;
                inds.push_back(i);
                break;
            };
        };
    }

    // For each of the 6 auto / cross spectra terms, map a slice out of the psi tensor
    for (Eigen::Index spectrum_ind = 0; spectrum_ind < 6; spectrum_ind++) {
        Figure fig_ls = Figure();
        for (auto &ind : inds) {
            std::cout << "Adding spectrum trace for index " << ind << std::endl;
            ScatterPlot p = ScatterPlot();
            std::cout << "data k1z size " << data.k1z.rows() << " x " << data.k1z.cols() << std::endl;
            p.x = data.k1z.col(ind);
            std::cout << "spectra[spectrum_ind] size " << spectra[spectrum_ind].rows() << " x "
                      << spectra[spectrum_ind].cols() << std::endl;
            p.y = spectra[spectrum_ind].col(ind);
            p.name = "S_" + term_labels[spectrum_ind] + " at \\eta = " + std::to_string(data.eta(ind));
            fig_ls.add(p);
        }
        Layout lay = Layout();
        lay.xTitle("k1z");
        lay.yTitle("psi");
        lay.xLog();
        fig_ls.setLayout(lay);
        fig_ls.write("test_spectrum_line_plot_" + term_labels[spectrum_ind] + ".json");
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


// Helper function to keep code dry whist plotting data
ScatterPlot getTrace(const ArrayXXd &data, const std::string &color, const std::string &dash, std::string name="") {
    ScatterPlot p = ScatterPlot();
    p.x = data.row(0);
    p.y = data.row(1);
    p.name = name;
    p.setDash(dash);
    p.setColor(color);
    return p;
}


// Helper function to run ADEM for Reynolds Stress comparisons
void getStressTraces(
    ScatterPlot &r11, ScatterPlot &r22, ScatterPlot &r33,
    const EddySignature &signature_a, const EddySignature &signature_b,
    const double pi_coles, const double shear_ratio, const double zeta, const double beta,
    const std::string marker, const std::string name, const std::string color
    ) {

    // TODO should use the K_tau values given in M&P with the shear ratio, to fix at least the ratio of delta_c to U1
    // Use these same parameters for all. The outputs are nondimensional so it shouldn't matter
    double delta_c = 1000.0;
    double kappa = 0.41;
    double u_inf = 2.2;

    AdemData data = adem(beta, delta_c, kappa, pi_coles, shear_ratio, u_inf, zeta, signature_a, signature_b, false);
    std::cout << "Computed stresses using ADEM for case " << name << std::endl;
    std::cout << data << std::endl;

    r11.x = data.eta;
    r11.y = data.reynolds_stress.col(0);
    r11.name = name;
    r11.setDash("solid");
    r11.setColor(color);

    std::cout << "trace1 " << name << std::endl;

    r22.x = data.eta;
    r22.y = data.reynolds_stress.col(3);
    r22.name = name;
    r22.setDash("solid");
    r22.setColor(color);

    std::cout << "trace2 " << name << std::endl;

    r33.x = data.eta;
    r33.y = data.reynolds_stress.col(5);
    r33.name = name;
    r33.setDash("solid");
    r33.setColor(color);
    std::cout << "trace3 " << name << std::endl;


}


TEST_F(AdemTest, test_validate_stresses_perry_marusic) {

    // Load the eddy signatures. We only want to do this once.
//    EddySignature signature_a = EddySignature();
//    signature_a.load(data_path + std::string("/signatures_A.mat"));
//    EddySignature signature_b = EddySignature();
//    signature_b.load(data_path + std::string("/signatures_B1.mat"));
//    EddySignature signature_bx = EddySignature();
//    signature_bx.load(data_path + std::string("/signatures_B2.mat"));
//    signature_b = signature_b + signature_bx;
//    signature_bx.load(data_path + std::string("/signatures_B3.mat"));
//    signature_b = signature_b + signature_bx;
//    signature_bx.load(data_path + std::string("/signatures_B4.mat"));
//    signature_b = signature_b + signature_bx;
//    signature_b = signature_b / 4.0;
    int n_lambda = 200;
    double dx = 0.01;
    EddySignature signature_a = EddySignature();
    signature_a.getEddyIntensity("A", n_lambda, dx);

    EddySignature signature_b = EddySignature();
    signature_b.getEddyIntensity("B1", n_lambda, dx);
    EddySignature sig = EddySignature();
    sig.getEddyIntensity("B2", n_lambda, dx);
    signature_b = signature_b + sig;
    sig.getEddyIntensity("B3", n_lambda, dx);
    signature_b = signature_b + sig;
    sig.getEddyIntensity("B4", n_lambda, dx);
    signature_b = signature_b + sig;
    Figure fig_i = Figure();

    // Plot the J and I functions
    std::vector<std::string> labels = {"11", "12", "13", "22", "23", "33"};

    std::vector<ScatterPlot> p;
    for (auto i=0; i<6; i++) {
        ScatterPlot pij = ScatterPlot();
        pij.x = signature_a.eta;
        pij.y = signature_a.j.col(i);
        pij.name = "$I_{"+labels[i]+"}$";
        p[i] = pij;
        fig_i.add(p[i]);
    }
    Layout lay = Layout();
    lay.xTitle("$\\eta$");
    lay.yTitle("$Type A Eddy Intensity I_{ij}");
    fig_i.setLayout(lay);
    fig_i.write("check_type_a_iij_from_cpp.json");

    Figure fig_j = Figure();
    std::vector<ScatterPlot> pa;
    for (auto i=0; i<6; i++) {
        ScatterPlot pij = ScatterPlot();

        pij.x = signature_a.lambda;
        pij.y = signature_a.j.col(i);
        pij.name = "$J_{"+labels[i]+"}$";
        pa[i] = pij;
        fig_j.add(pa[i]);
    }
    Layout lay_j = Layout();
    lay_j.xTitle("$\\lambda$");
    lay_j.yTitle("$Eddy Intensity J_{ij}");
    fig_j.setLayout(lay);
    fig_j.write("check_Jij_from_cpp.json");
    
    // Run ADEM for all measurement locations in the "10APG" flow case in Perry and Marusic 1995
//    ScatterPlot r11_trace_1 = ScatterPlot();
//    ScatterPlot r22_trace_1 = ScatterPlot();
//    ScatterPlot r33_trace_1 = ScatterPlot();
//    getStressTraces(r11_trace_1, r22_trace_1, r33_trace_1, signature_a, signature_b, 3.23, 38.4, 15.32, 7.16, "diamond", "es-flow, $R_\\theta = 7257$", "#1f77b4");
//
//    ScatterPlot r11_trace_2 = ScatterPlot();
//    ScatterPlot r22_trace_2 = ScatterPlot();
//    ScatterPlot r33_trace_2 = ScatterPlot();
//    getStressTraces(r11_trace_2, r22_trace_2, r33_trace_2, signature_a, signature_b, 2.46, 34.5, 8.01, 4.48, "lefttriangle", "es-flow, $R_\\theta = 6395$", "#ff7f0e");

    ScatterPlot r11_trace_3 = ScatterPlot();
    ScatterPlot r22_trace_3 = ScatterPlot();
    ScatterPlot r33_trace_3 = ScatterPlot();
    getStressTraces(r11_trace_3, r22_trace_3, r33_trace_3, signature_a, signature_b, 1.87, 31.5, 4.64, 2.90, "uptriangle", "es-flow, $R_\\theta = 5395$", "#2ca02c");

//    ScatterPlot r11_trace_4 = ScatterPlot();
//    ScatterPlot r22_trace_4 = ScatterPlot();
//    ScatterPlot r33_trace_4 = ScatterPlot();
//    getStressTraces(r11_trace_4, r22_trace_4, r33_trace_4, signature_a, signature_b, 1.19, 28.1, 2.18, 1.45, "square", "es-flow, $R_\\theta = 4155$", "#d62728");
//
//    ScatterPlot r11_trace_5 = ScatterPlot();
//    ScatterPlot r22_trace_5 = ScatterPlot();
//    ScatterPlot r33_trace_5 = ScatterPlot();
//    getStressTraces(r11_trace_5, r22_trace_5, r33_trace_5, signature_a, signature_b, 0.68, 25.4, 0.94, 0.65, "downtriangle", "es-flow, $R_\\theta = 3153$", "#9467bd");
//
//    ScatterPlot r11_trace_6 = ScatterPlot();
//    ScatterPlot r22_trace_6 = ScatterPlot();
//    ScatterPlot r33_trace_6 = ScatterPlot();
//    getStressTraces(r11_trace_6, r22_trace_6, r33_trace_6, signature_a, signature_b, 0.42, 23.6, 0.15, 0.0, "circle", "es-flow, $R_\\theta = 2206$", "#e377c2");

    
    // Digitally extracted from Perry and Marusic 1995 Part 1, Figure 6a:
    ArrayXXd fig_6a_trace_1 = Eigen::ArrayXXd(2, 32);
    fig_6a_trace_1 << 0.0100961,   0.012242,   0.0152986,   0.0194686,   0.0246286,   0.0330971,   0.0450213,   0.0562821,   0.0716523,   0.0879855,   0.108698,   0.138418,   0.167934,   0.206221,   0.24866,   0.285682,   0.32615,   0.376796,   0.422329,   0.456479,   0.505256,   0.542659,   0.579328,   0.610986,   0.644355,   0.679589,   0.712374,   0.769516,   0.811594,   0.871647,   0.936379,   1,
        9.54832,   9.38548,   9.25332,   9.1837,   9.20878,   9.38977,   9.66497,   9.87945,   10.2195,   10.6235,   11.0589,   11.6511,   12.0555,   12.491,   12.7695,   12.9233,   12.7936,   12.4745,   12.0302,   11.4293,   10.4812,   9.56531,   8.64963,   7.60823,   6.53532,   5.52543,   4.42117,   3.22155,   2.21166,   1.2643,   0.569042,   0;

    ArrayXXd fig_6a_trace_2 = Eigen::ArrayXXd(2, 31);
    fig_6a_trace_2 << 0.0100243, 0.0117941, 0.0144739, 0.0178717, 0.0224683, 0.0261233, 0.0314893, 0.039358, 0.0471635, 0.0614952, 0.0741318, 0.0904717, 0.111074, 0.133119, 0.160509, 0.190042, 0.227753, 0.274562, 0.332888, 0.389232, 0.452207, 0.500637, 0.547509, 0.584617, 0.62039, 0.662395, 0.7072, 0.755083, 0.796598, 0.850668, 0.993959,
        8.44554, 8.28353, 8.08884, 7.98853, 7.8247, 7.85208, 7.81546, 7.84086, 7.93045, 8.04924, 8.07564, 8.35381, 8.56879, 8.78443, 9.03141, 9.15284, 9.33697, 9.39488, 9.1375, 8.78657, 8.05766, 7.33006, 6.41372, 5.68712, 4.80311, 4.01348, 3.16083, 2.3712, 1.64492, 1.01286,     -0.0313476;

    ArrayXXd fig_6a_trace_3 = Eigen::ArrayXXd(2, 25);
    fig_6a_trace_3 << 0.0100779, 0.0120732, 0.0149951, 0.0189659, 0.0238439, 0.0308948, 0.0410114, 0.0544391, 0.0705484, 0.0881828, 0.113601, 0.14286, 0.182939, 0.228639, 0.285736, 0.346414, 0.414853, 0.476197, 0.539977, 0.59775, 0.657688, 0.727964, 0.786704, 0.860631, 1,
        7.7521, 7.55808, 7.30003, 7.13604, 6.97221, 6.83906, 6.79979, 6.72901, 6.75343, 6.84187, 6.96099, 7.08077, 7.23157, 7.19395, 7.09331, 6.7729, 6.23224, 5.50365, 4.61783, 3.7957, 2.91071, 1.96253, 1.20408, 0.602864, 0;

    ArrayXXd fig_6a_trace_4 = Eigen::ArrayXXd(2, 25);
    fig_6a_trace_4 << 0.0100713, 0.0119195, 0.0148937, 0.0188364, 0.0246986, 0.0304938, 0.0376499, 0.0487848, 0.0624535, 0.0789988, 0.0975501, 0.122652, 0.150549, 0.185896, 0.218729, 0.273326, 0.341475, 0.401661, 0.483846, 0.572313, 0.649069, 0.736166, 0.839997, 0.924808, 1,
        7.09034, 6.8021, 6.5439, 6.31687, 6.02584, 5.83099, 5.66765, 5.56602, 5.4332, 5.36374, 5.32646, 5.25716, 5.25155, 5.18275, 5.08376, 4.88858, 4.47281, 4.05869, 3.3603, 2.50484, 1.77658, 1.11136, 0.445961, 0.19122, 0;

    ArrayXXd fig_6a_trace_5 = Eigen::ArrayXXd(2, 22);
    fig_6a_trace_5 << 0.0100665, 0.0124273, 0.0158106, 0.0204827, 0.0268564, 0.0362919, 0.0467368, 0.0598279, 0.0765858, 0.0992298, 0.129358, 0.168614, 0.211995, 0.266495, 0.325051, 0.408513, 0.48922, 0.582315, 0.660475, 0.727023, 0.814982, 1,
        6.61765, 6.32826, 6.00653, 5.71583, 5.39328, 5.10141, 4.87389, 4.67805, 4.48221, 4.31756, 4.24727, 4.05094, 3.95013, 3.69176, 3.4027, 2.89223, 2.35157, 1.74804, 1.11433, 0.670509, 0.352248, 0;

    ArrayXXd fig_6a_trace_6 = Eigen::ArrayXXd(2, 21);
    fig_6a_trace_6 << 0.0100627, 0.0122729, 0.0151511, 0.0196277, 0.0249714, 0.0331358, 0.0434482, 0.0612395, 0.0803008, 0.109829, 0.152043, 0.193449, 0.247626, 0.320801, 0.396022, 0.471516, 0.561313, 0.656237, 0.776591, 0.902692, 1,
        6.2395, 5.88741, 5.59802, 5.2758, 4.95408, 4.53666, 4.24562, 3.88958, 3.63005, 3.33786, 3.07685, 2.81815, 2.5908, 2.30009, 1.97919, 1.65928, 1.18181, 0.704826, 0.322056, 0.0973426, 0;


    // Digitally extracted from Perry and Marusic 1995 Part 1, Figure 6b:
    ArrayXXd fig_6b_trace_1 = Eigen::ArrayXXd(2, 32);
    fig_6b_trace_1 << 0.0100224, 0.0122273, 0.0149175, 0.0184207, 0.0222046, 0.025661, 0.0316882, 0.0400868, 0.0501049, 0.0618777, 0.0755003, 0.0932444, 0.116558, 0.141373, 0.171473, 0.207977, 0.233234, 0.274473, 0.330876, 0.382356, 0.433898, 0.483524, 0.53235, 0.57557, 0.61855, 0.660734, 0.714332, 0.758475, 0.824999, 0.87607, 0.953028, 1,
        4.76633, 4.63086, 4.52908, 4.46094, 4.39289, 4.4092, 4.40843, 4.4581, 4.55834, 4.70915, 4.84315, 5.09502, 5.39736, 5.66613, 5.95174, 6.20367, 6.38851, 6.52265, 6.58934, 6.47091, 6.20098, 5.76268, 5.24022, 4.59994, 3.92598, 3.21837, 2.44335, 1.76944, 1.09545, 0.640487, 0.235966, 0;

    ArrayXXd fig_6b_trace_2 = Eigen::ArrayXXd(2, 28);
    fig_6b_trace_2 << 0.010021, 0.0119345, 0.0154645, 0.0193269, 0.0247441, 0.032258, 0.0430819, 0.0575372, 0.0773105, 0.0989925, 0.125994, 0.156538, 0.187581, 0.230263, 0.279257, 0.330591, 0.384334, 0.433528, 0.486049, 0.52875, 0.571723, 0.618184, 0.660416, 0.709782, 0.762838, 0.829823, 0.891929, 1,
        4.48001, 4.36148, 4.19211, 4.03971, 3.93776, 3.85258, 3.86836, 3.8673, 3.95043, 4.11795, 4.28549, 4.43628, 4.60404, 4.75487, 4.80469, 4.75354, 4.61825, 4.38202, 4.01108, 3.6234, 3.15153, 2.66282, 2.191, 1.66863, 1.14626, 0.674369, 0.337262, 0;

    ArrayXXd fig_6b_trace_3 = Eigen::ArrayXXd(2, 29);
    fig_6b_trace_3 << 0.00996, 0.0117196, 0.0143841, 0.0177616, 0.0220646, 0.0220646, 0.0279101, 0.0363852, 0.0477209, 0.0618401, 0.0777582, 0.100767, 0.125948, 0.163221, 0.207733, 0.251929, 0.303685, 0.361667, 0.410432, 0.457415, 0.506708, 0.561283, 0.610583, 0.664209, 0.718192, 0.7813, 0.844845, 0.913631, 1.00001,
        4.31161, 4.17628, 4.00711, 3.88844, 3.75291, 3.75291, 3.61732, 3.51529, 3.43009, 3.4123, 3.41146, 3.4442, 3.51075, 3.61086, 3.69418, 3.71032, 3.67595, 3.52374, 3.32117, 3.06814, 2.79829, 2.41054, 1.98918, 1.55098, 1.07911, 0.725116, 0.371145, 0.185594, 0.0168422;

    ArrayXXd fig_6b_trace_4 = Eigen::ArrayXXd(2, 23);
    fig_6b_trace_4 << 0.00995929, 0.0114395, 0.0132192, 0.0160306, 0.0206468, 0.0253411, 0.0318616, 0.0407915, 0.0528585, 0.0714453, 0.0948426, 0.123645, 0.162171, 0.211422, 0.270686, 0.346541, 0.417694, 0.503425, 0.585181, 0.672069, 0.758026, 0.855038, 1.00605,
        4.16003, 4.00795, 3.87268, 3.73724, 3.51736, 3.36504, 3.19578, 3.06014, 2.95813, 2.80545, 2.77073, 2.71924, 2.7014, 2.66675, 2.59847, 2.39547, 2.15899, 1.78778, 1.34933, 0.927773, 0.506279, 0.236363, 0;

    ArrayXXd fig_6b_trace_5 = Eigen::ArrayXXd(2, 21);
    fig_6b_trace_5 << 0.00995866, 0.0119317, 0.0146443, 0.0194381, 0.0262717, 0.0348725, 0.0449148, 0.0596205, 0.0835522, 0.112934, 0.149012, 0.199003, 0.262575, 0.340237, 0.422646, 0.515576, 0.581559, 0.659942, 0.744407, 0.829653, 1.00605,
        4.0253, 3.85621, 3.6702, 3.41653, 3.14596, 2.94281, 2.73978, 2.58717, 2.41751, 2.29851, 2.19645, 2.11118, 1.99227, 1.80606, 1.58631, 1.2319, 0.961987, 0.658365, 0.405292, 0.236474, 0;

    ArrayXXd fig_6b_trace_6 = Eigen::ArrayXXd(2, 21);
    fig_6b_trace_6 << 0.00995819, 0.0128256, 0.0167187, 0.0214031, 0.0291022, 0.0417765, 0.0568056, 0.0820406, 0.110887, 0.146307, 0.201365, 0.250153, 0.296128, 0.365661, 0.420007, 0.503206, 0.574513, 0.667872, 0.744342, 0.854977, 1,
        3.92424, 3.67069, 3.3834, 3.11302, 2.82558, 2.5211, 2.28418, 2.04705, 1.86068, 1.69125, 1.57218, 1.47034, 1.35183, 1.23316, 1.08107, 0.861463, 0.675714, 0.422531, 0.220028, 0.0847839, 0;


    // Digitally extracted from Perry and Marusic 1995 Part 1, Figure 6c:
    ArrayXXd fig_6c_trace_1 = Eigen::ArrayXXd(2, 33);
    fig_6c_trace_1 << 0, 0.0183246, 0.0366492, 0.0549738, 0.0746073, 0.0968586, 0.125654, 0.157068, 0.188482, 0.219895, 0.248691, 0.280105, 0.320681, 0.36911, 0.40445, 0.447644, 0.482984, 0.515707, 0.548429, 0.575916, 0.604712, 0.632199, 0.667539, 0.698953, 0.727749, 0.760471, 0.798429, 0.829843, 0.859948, 0.890052, 0.922775, 0.962042, 0.989529,
        0.995807, 1.58281, 2.02306, 2.45283, 2.80922, 3.12369, 3.48008, 3.7631, 4.00419, 4.15094, 4.27673, 4.36059, 4.38155, 4.32914, 4.20335, 4.01468, 3.80503, 3.57442, 3.30189, 3.02935, 2.73585, 2.44235, 2.09644, 1.78197, 1.5304, 1.2369, 0.974843, 0.754717, 0.597484, 0.461216, 0.314465, 0.199161, 0.146751;

    ArrayXXd fig_6c_trace_2 = Eigen::ArrayXXd(2, 26);
    fig_6c_trace_2 << 0, 0.0183246, 0.0366492, 0.0589005, 0.0837696, 0.112565, 0.14267, 0.179319, 0.227749, 0.280105, 0.331152, 0.382199, 0.424084, 0.479058, 0.524869, 0.566754, 0.606021, 0.650524, 0.700262, 0.748691, 0.786649, 0.82199, 0.861257, 0.908377, 0.950262, 1,
        0.974843, 1.41509, 1.76101, 2.07547, 2.40042, 2.64151, 2.85115, 3.02935, 3.18658, 3.27044, 3.23899, 3.12369, 2.96646, 2.74633, 2.44235, 2.18029, 1.88679, 1.56184, 1.22642, 0.943396, 0.72327, 0.545073, 0.398323, 0.24109, 0.146751, 0.0733753;

    ArrayXXd fig_6c_trace_3 = Eigen::ArrayXXd(2, 27);
    fig_6c_trace_3 << 0, 0.0196335, 0.0366492, 0.0641361, 0.0942408, 0.128272, 0.162304, 0.197644, 0.246073, 0.294503, 0.341623, 0.388743, 0.437173, 0.468586, 0.513089, 0.557592, 0.599476, 0.599476, 0.640052, 0.679319, 0.727749, 0.777487, 0.825916, 0.871728, 0.905759, 0.942408, 0.998691,
        0.995807, 1.3522, 1.58281, 1.87631, 2.10692, 2.27463, 2.42138, 2.51572, 2.55765, 2.56813, 2.49476, 2.36897, 2.21174, 2.07547, 1.87631, 1.65618, 1.42558, 1.42558, 1.20545, 0.995807, 0.754717, 0.555556, 0.377358, 0.262055, 0.178197, 0.104822, 0;

    ArrayXXd fig_6c_trace_4 = Eigen::ArrayXXd(2, 21);
    fig_6c_trace_4 << 0, 0.0235602, 0.0536649, 0.0876963, 0.13089, 0.174084, 0.212042, 0.264398, 0.323298, 0.374346, 0.446335, 0.498691, 0.557592, 0.608639, 0.660995, 0.71466, 0.763089, 0.810209, 0.865183, 0.91623, 1.00131,
        0.995807, 1.26834, 1.5304, 1.67715, 1.80294, 1.87631, 1.88679, 1.87631, 1.77149, 1.68763, 1.47799, 1.29979, 1.10063, 0.880503, 0.691824, 0.513627, 0.387841, 0.272537, 0.157233, 0.0733753, 0;

    ArrayXXd fig_6c_trace_5 = Eigen::ArrayXXd(2, 20);
    fig_6c_trace_5 << 0.0013089, 0.0379581, 0.0732984, 0.112565, 0.158377, 0.197644, 0.246073, 0.29712, 0.350785, 0.412304, 0.472513, 0.527487, 0.586387, 0.647906, 0.704188, 0.76178, 0.819372, 0.873037, 0.945026, 1,
        0.964361, 1.26834, 1.40461, 1.46751, 1.54088, 1.5304, 1.50943, 1.44654, 1.34172, 1.20545, 1.04822, 0.91195, 0.712788, 0.545073, 0.408805, 0.283019, 0.167715, 0.0943396, 0.0524109, 0.0104822;

    ArrayXXd fig_6c_trace_6 = Eigen::ArrayXXd(2, 22);
    fig_6c_trace_6 << 0, 0.0301047, 0.0693717, 0.103403, 0.143979, 0.188482, 0.227749, 0.270942, 0.318063, 0.352094, 0.399215, 0.447644, 0.497382, 0.541885, 0.589005, 0.638743, 0.685864, 0.744764, 0.798429, 0.853403, 0.912304, 1,
        1.00629, 1.12159, 1.19497, 1.21593, 1.22642, 1.20545, 1.16352, 1.12159, 1.04822, 1.00629, 0.922432, 0.828092, 0.712788, 0.618449, 0.513627, 0.429769, 0.324948, 0.220126, 0.146751, 0.104822, 0.0628931, 0;


    // Recreate P&M Figure 6a
    ScatterPlot trace1_6a = getTrace(fig_6a_trace_1, "#1f77b4", "dash");
    ScatterPlot trace2_6a = getTrace(fig_6a_trace_2, "#ff7f0e", "dash");
    ScatterPlot trace3_6a = getTrace(fig_6a_trace_3, "#2ca02c", "dash");
    ScatterPlot trace4_6a = getTrace(fig_6a_trace_4, "#d62728", "dash");
    ScatterPlot trace5_6a = getTrace(fig_6a_trace_5, "#9467bd", "dash");
    ScatterPlot trace6_6a = getTrace(fig_6a_trace_6, "#e377c2", "dash");
    Figure fig_6a = Figure();

    std::cout << "trace4 " << std::endl;

    // Traces from the paper
//    fig_6a.add(trace1_6a);
//    fig_6a.add(trace2_6a);
    fig_6a.add(trace3_6a);
    std::cout << "trace5 " << std::endl;
//    fig_6a.add(trace4_6a);
//    fig_6a.add(trace5_6a);
//    fig_6a.add(trace6_6a);

    // Traces recreated with es-flow
//    fig_6a.add(r11_trace_1);
//    fig_6a.add(r11_trace_2);
    fig_6a.add(r11_trace_3);
    std::cout << "trac6" << std::endl;
//    fig_6a.add(r11_trace_4);
//    fig_6a.add(r11_trace_5);
//    fig_6a.add(r11_trace_6);
    Layout lay_6a = Layout();
    lay_6a.xLog();
    lay_6a.xTitle("$z/\\delta_{c}$");
    lay_6a.yTitle("$\\overline{u_1^2} / U_{\\tau}^2$");
    fig_6a.setLayout(lay_6a);
    fig_6a.write("validation_perry_marusic_6a.json");


    // Recreate P&M Figure 6b
    ScatterPlot trace1_6b = getTrace(fig_6b_trace_1, "#1f77b4", "dash");
    ScatterPlot trace2_6b = getTrace(fig_6b_trace_2, "#ff7f0e", "dash");
    ScatterPlot trace3_6b = getTrace(fig_6b_trace_3, "#2ca02c", "dash");
    ScatterPlot trace4_6b = getTrace(fig_6b_trace_4, "#d62728", "dash");
    ScatterPlot trace5_6b = getTrace(fig_6b_trace_5, "#9467bd", "dash");
    ScatterPlot trace6_6b = getTrace(fig_6b_trace_6, "#e377c2", "dash");
    Figure fig_6b = Figure();

    // Traces from the paper
//    fig_6b.add(trace1_6b);
//    fig_6b.add(trace2_6b);
    fig_6b.add(trace3_6b);
//    fig_6b.add(trace4_6b);
//    fig_6b.add(trace5_6b);
//    fig_6b.add(trace6_6b);

    // Traces recreated with es-flow
//    fig_6b.add(r22_trace_1);
//    fig_6b.add(r22_trace_2);
    fig_6b.add(r22_trace_3);
//    fig_6b.add(r22_trace_4);
//    fig_6b.add(r22_trace_5);
//    fig_6b.add(r22_trace_6);

    Layout lay_6b = Layout();
    lay_6b.xLog();
    lay_6b.xTitle("$z/\\delta_{c}$");
    lay_6b.yTitle("$\\overline{u_2^2} / U_{\\tau}^2$");
    fig_6b.setLayout(lay_6b);
    fig_6b.write("validation_perry_marusic_6b.json");


    // Recreate P&M Figure 6c
    ScatterPlot trace1_6c = getTrace(fig_6c_trace_1, "#1f77b4", "dash");
    ScatterPlot trace2_6c = getTrace(fig_6c_trace_2, "#ff7f0e", "dash");
    ScatterPlot trace3_6c = getTrace(fig_6c_trace_3, "#2ca02c", "dash");
    ScatterPlot trace4_6c = getTrace(fig_6c_trace_4, "#d62728", "dash");
    ScatterPlot trace5_6c = getTrace(fig_6c_trace_5, "#9467bd", "dash");
    ScatterPlot trace6_6c = getTrace(fig_6c_trace_6, "#e377c2", "dash");
    Figure fig_6c = Figure();

    // Traces from the paper
//    fig_6c.add(trace1_6c);
//    fig_6c.add(trace2_6c);
    fig_6c.add(trace3_6c);
//    fig_6c.add(trace4_6c);
//    fig_6c.add(trace5_6c);
//    fig_6c.add(trace6_6c);

    // Traces recreated with es-flow
//    fig_6c.add(r33_trace_1);
//    fig_6c.add(r33_trace_2);
    fig_6c.add(r33_trace_3);
//    fig_6c.add(r33_trace_4);
//    fig_6c.add(r33_trace_5);
//    fig_6c.add(r33_trace_6);

    Layout lay_6c = Layout();
    lay_6c.xTitle("$z/\\delta_{c}$");
    lay_6c.yTitle("$\\overline{u_3^2} / U_{\\tau}^2$");
    fig_6c.setLayout(lay_6c);
    fig_6c.write("validation_perry_marusic_6c.json");

}



TEST_F(AdemTest, test_validate_spectra_perry_marusic) {

    // Run ADEM for the flow case 1 ("10APG") case in Perry and Marusic 1995
    double beta = 0.0;
    double delta_c = 1000.0;
    double kappa = 0.41;
    double pi_coles = 0.42;
    double shear_ratio = 23.6;
    double u_inf = 2.2;
    double zeta = 0.15;

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
    signature_b = signature_b / 4.0;

    AdemData data = adem(beta, delta_c, kappa, pi_coles, shear_ratio, u_inf, zeta, signature_a, signature_b);
    std::cout << data << std::endl;

    // Digitally extracted from Perry and Marusic 1995 Part 1, Figure 12a:
    ArrayXXd fig_12a_trace_1 = Eigen::ArrayXXd(2, 32);
    fig_12a_trace_1 << 4.12E-05, 8.02E-05, 0.000148271, 0.000234852, 0.000398128, 0.000640683, 0.000996147, 0.0013747, 0.0013747, 0.00196244, 0.00248839, 0.00315552, 0.00443334, 0.0060219, 0.00922047, 0.0146137, 0.0284457, 0.0436146, 0.0741215, 0.130451, 0.177631, 0.237797, 0.318387, 0.41189, 0.580655, 0.804595, 1.1533, 1.76971, 3.33369, 5.57179, 10.1384, 20.7821,
        0.0260674, 0.0575518, 0.100382, 0.18577, 0.302433, 0.484417, 0.686276, 0.868206, 0.868206, 1.07288, 1.19796, 1.31736, 1.43396, 1.51929, 1.58478, 1.62187, 1.61075, 1.56828, 1.48893, 1.34424, 1.23639, 1.12286, 0.997966, 0.890103, 0.739657, 0.600568, 0.484219, 0.379257, 0.260163, 0.149549, 0.0702154, 0.0193324;

    ArrayXXd fig_12a_trace_2 = Eigen::ArrayXXd(2, 27);
    fig_12a_trace_2 << 0.000435102, 0.000804204, 0.00146089, 0.0023544, 0.00354252, 0.00506441, 0.00642427, 0.008429, 0.0112502, 0.0147582, 0.0193601, 0.0258408, 0.0316697, 0.0510195, 0.0808822, 0.142211, 0.196914, 0.254651, 0.365146, 0.506135, 0.654988, 1.3455, 2.44975, 4.23649, 9.14984, 18.7578, 43.3119,
        0.0325798, 0.0640461, 0.115392, 0.180901, 0.277635, 0.36867, 0.462503, 0.587598, 0.709858, 0.849157, 0.988456, 1.10787, 1.22442, 1.32118, 1.33839, 1.27041, 1.18814, 1.10868, 0.963923, 0.799267, 0.665836, 0.447342, 0.319713, 0.211952, 0.0815427, 0.0221372, 0.0082277;

    ArrayXXd fig_12a_trace_3 = Eigen::ArrayXXd(2, 24);
    fig_12a_trace_3 << 0.00205775, 0.00414231, 0.00727267, 0.0125498, 0.0198724, 0.0309349, 0.0473061, 0.067595, 0.09494, 0.126767, 0.211421, 0.353091, 0.610885, 0.890521, 1.43805, 2.28183, 3.5594, 3.5594, 5.01454, 6.48601, 9.29529, 13.7856, 19.7415, 46.3733,
        0.0302866, 0.061783, 0.107436, 0.170127, 0.278243, 0.386353, 0.551274, 0.682081, 0.821404, 0.912414, 0.969412, 0.918457, 0.776606, 0.677308, 0.586568, 0.532754, 0.478934, 0.478934, 0.379624, 0.285965, 0.18382, 0.0731646, 0.0306778, 0.00825178;

    ArrayXXd fig_12a_trace_4 = Eigen::ArrayXXd(2, 32);
    fig_12a_trace_4 << 0.00407396, 0.00593025, 0.0104148, 0.0154179, 0.0220544, 0.0320872, 0.047469, 0.0655862, 0.0846449, 0.105586, 0.13859, 0.178843, 0.234856, 0.293107, 0.426738, 0.610951, 0.776234, 1.09254, 1.51144, 2.12703, 2.94331, 3.80507, 4.75485, 5.94215, 7.05486, 8.09591, 9.44955, 11.8078, 14.5038, 16.3558, 20.7806, 30.7791,
        0.0276865, 0.0420233, 0.0649489, 0.096337, 0.141917, 0.196026, 0.281391, 0.369573, 0.449207, 0.520307, 0.614153, 0.70231, 0.759224, 0.790552, 0.790684, 0.768084, 0.736918, 0.711471, 0.700221, 0.686137, 0.655002, 0.601116, 0.533013, 0.459229, 0.388267, 0.30593, 0.232121, 0.166859, 0.107273, 0.0561795, 0.0250141, 0.01663;

    ArrayXXd fig_12a_trace_5 = Eigen::ArrayXXd(2, 23);
    fig_12a_trace_5 << 0.00924652, 0.0192616, 0.0349961, 0.0563987, 0.0924389, 0.139087, 0.198797, 0.289064, 0.386092, 0.58118, 1.02023, 1.70129, 2.69709, 3.9968, 5.25675, 6.46026, 8.21688, 10.2791, 12.0038, 14.5023, 17.8154, 23.0248, 41.8533,
        0.0222936, 0.0481203, 0.0852624, 0.153612, 0.236172, 0.332906, 0.440986, 0.540549, 0.605991, 0.668635, 0.725651, 0.794012, 0.81122, 0.763064, 0.69782, 0.598462, 0.48207, 0.328741, 0.21516, 0.115796, 0.047687, 0.0165276, 0.0167383;

    ArrayXXd fig_12a_trace_6 = Eigen::ArrayXXd(2, 23);
    fig_12a_trace_6 << 0.0236457, 0.0444573, 0.0793951, 0.125789, 0.186143, 0.270724, 0.400488, 0.58238, 0.846729, 1.27406, 1.82174, 2.27293, 2.93554, 3.72944, 4.98921, 6.23683, 7.53905, 9.27112, 11.2049, 14.7493, 17.8141, 21.5025, 37.1384,
        0.0339882, 0.0654605, 0.11396, 0.179463, 0.2421, 0.324618, 0.412824, 0.506705, 0.614791, 0.708684, 0.785514, 0.839569, 0.862386, 0.836902, 0.777346, 0.680835, 0.538857, 0.388363, 0.26059, 0.130006, 0.0533687, 0.0250262, 0.0166962;

    ArrayXXd fig_12a_trace_7 = Eigen::ArrayXXd(2, 20);
    fig_12a_trace_7 << 0.044475, 0.0942317, 0.177111, 0.30558, 0.554582, 0.940007, 1.5657, 2.20044, 2.74582, 3.48663, 4.42924, 5.25694, 6.34859, 7.67054, 9.43214, 11.798, 13.084, 16.3564, 19.7444, 38.4299,
        0.0342109, 0.0685662, 0.125606, 0.199662, 0.324871, 0.452897, 0.612166, 0.703195, 0.745886, 0.760175, 0.740373, 0.694979, 0.626864, 0.521818, 0.377006, 0.232199, 0.152691, 0.0533386, 0.0193143, 0.0138673;

    ArrayXXd fig_12a_trace_8 = Eigen::ArrayXXd(2, 16);
    fig_12a_trace_8 << 0.0822453, 0.160028, 0.333283, 0.648106, 1.08055, 1.65402, 2.61997, 3.8783, 5.74598, 7.30284, 8.8203, 11.0239, 12.8644, 15.8044, 19.0788, 30.7802,
        0.0259049, 0.0460258, 0.0888977, 0.154473, 0.237039, 0.325256, 0.407804, 0.444874, 0.413763, 0.357029, 0.280392, 0.198085, 0.141321, 0.067531, 0.0306658, 0.0137891;

    ArrayXXd fig_12a_trace_9 = Eigen::ArrayXXd(2, 9);
    fig_12a_trace_9 << 0.357138, 0.649068, 1.12048, 1.96737, 3.51487, 5.96778, 9.15017, 14.5159, 25.0759,
        0.0235817, 0.0379967, 0.0694389, 0.10941, 0.12666, 0.121165, 0.0787019, 0.041933, 0.0193986;


    // Verification data
    Eigen::Tensor<double, 3> psi;
    psi = data.psi_a + data.psi_b;


    // Recreate P&M Figure 12, composite spectra (the bottom part)
    ScatterPlot trace1_12a = getTrace(fig_12a_trace_1, "#1f77b4", "dash", "$z/\\delta_{c}=0.01$");
    ScatterPlot trace2_12a = getTrace(fig_12a_trace_2, "#ff7f0e", "dash", "$z/\\delta_{c}=0.05$");
    ScatterPlot trace3_12a = getTrace(fig_12a_trace_3, "#2ca02c", "dash", "$z/\\delta_{c}=0.10$");
    ScatterPlot trace4_12a = getTrace(fig_12a_trace_4, "#d62728", "dash", "$z/\\delta_{c}=0.17$");
    ScatterPlot trace5_12a = getTrace(fig_12a_trace_5, "#9467bd", "dash", "$z/\\delta_{c}=0.27$");
    ScatterPlot trace6_12a = getTrace(fig_12a_trace_6, "#e377c2", "dash", "$z/\\delta_{c}=0.39$");
    ScatterPlot trace7_12a = getTrace(fig_12a_trace_7, "#7f7f7f", "dash", "$z/\\delta_{c}=0.54$");
    ScatterPlot trace8_12a = getTrace(fig_12a_trace_8, "#bcbd22", "dash", "$z/\\delta_{c}=0.72$");
    ScatterPlot trace9_12a = getTrace(fig_12a_trace_9, "#17becf", "dash", "$z/\\delta_{c}=0.93$");
    Figure fig_12a = Figure();
    fig_12a.add(trace1_12a);
    fig_12a.add(trace2_12a);
    fig_12a.add(trace3_12a);
    fig_12a.add(trace4_12a);
    fig_12a.add(trace5_12a);
    fig_12a.add(trace6_12a);
    fig_12a.add(trace7_12a);
    fig_12a.add(trace8_12a);
    fig_12a.add(trace9_12a);
    Layout lay_12a = Layout();
    lay_12a.xLog();
    lay_12a.xTitle("$k_{1}z$");
    lay_12a.yTitle("$k_{1}z\\Phi_{11}[k_{1}z] / U_\\tau^2$");
    fig_12a.setLayout(lay_12a);
    fig_12a.write("validation_perry_marusic_12a.json");

}
