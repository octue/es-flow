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
#include <math.h>
#include "gtest/gtest.h"

#include "cpplot.h"

#include "adem/adem.h"

#include "definitions.h"
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


/** @brief Helper to plot eddy intensities for different eddy types
 *
 * @param sig
 * @param type
 */
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


/** @brief Helper function to run ADEM for Reynolds Stress comparisons
 *
 * @param r11
 * @param r22
 * @param r33
 * @param r13
 * @param signature_a
 * @param signature_b
 * @param pi_coles
 * @param shear_ratio
 * @param zeta
 * @param beta
 * @param dash
 * @param name
 * @param color
 * @param type
 */
void getStressTraces(
    ScatterPlot &r11, ScatterPlot &r22, ScatterPlot &r33, ScatterPlot &r13,
    const EddySignature &signature_a, const EddySignature &signature_b,
    const double pi_coles, const double shear_ratio, const double zeta, const double beta,
    const std::string &dash, const std::string &name, const std::string &color,
    const std::string &type = "composite"
) {

    // TODO should use the K_tau values given in M&P with the shear ratio, to fix at least the ratio of delta_c to U1
    // Use these same parameters for all. The outputs are nondimensional so it shouldn't matter
    double delta_c = 1000.0;
    double kappa = KAPPA_VON_KARMAN;
    double u_inf = 2.2;

    AdemData data = adem(beta, delta_c, kappa, pi_coles, shear_ratio, u_inf, zeta, signature_a, signature_b, false);
    std::cout << "Computed stresses using ADEM for case " << name << std::endl;
    std::cout << data << std::endl;
    Eigen::ArrayXXd stress;
    if (type == "type_a") {
        stress = data.reynolds_stress_a;
        std::cout << "Using type a reynolds stresses" << std::endl;
    } else if (type == "type_b") {
        stress = data.reynolds_stress_b;
        std::cout << "Using type b reynolds stresses" << std::endl;
    } else {
        // type == composite by default
        stress = data.reynolds_stress;
        std::cout << "Using composite (type a + b) reynolds stresses" << std::endl;
    }
    r11.x = data.eta;
    r11.y = stress.col(0);
    r11.name = name;
    r11.setDash(dash);
    r11.setColor(color);

    std::cout << "trace1 " << name << std::endl;

    r22.x = data.eta;
    r22.y = stress.col(3);
    r22.name = name;
    r22.setDash(dash);
    r22.setColor(color);

    std::cout << "trace2 " << name << std::endl;

    r33.x = data.eta;
    r33.y = stress.col(5);
    r33.name = name;
    r33.setDash(dash);
    r33.setColor(color);
    std::cout << "trace3 " << name << std::endl;

    r13.x = data.eta;
    r13.y = stress.col(2);
    r13.name = name;
    r13.setDash(dash);
    r13.setColor(color);
    std::cout << "trace4 " << name << std::endl;

}


/** @brief Helper to make all 5 signature files, if no files already present
* TODO refactor out to signatures.h
 *
 * Signature files are named as "<data_path>/
 *
 * @param[in] data_path, std::string the path to the folder where signature files will be stored
 * @param[in] is_coarse, bool If true (the default), make coarse-grained signatures which are very quick for the purposes of unit testing, but which won't validate correctly.
 */
void make_signatures(const std::string &data_path, const bool is_coarse=true) {
    int n_lambda = 400;
    double dx = 0.005;
    std::string level = "fine";
    if (is_coarse) {
        n_lambda = 200;
        dx = 0.01;
        level = "coarse";
    }
    std::vector<std::string> types = {"A", "B1", "B2", "B3", "B4"};
    for (auto type : types) {
        EddySignature sig = EddySignature();
        sig.computeSignature(type, n_lambda, dx);
        sig.save(data_path + std::string("/signatures_"+type+"_"+level+".mat"));
    }
};


/** @brief Helper function to load eddy signatures for type A and ensembled Type B eddies
 *
 * @param sig_a, EddySignature to which the type A eddy data is loaded
 * @param sig_b
 * @param data_path
 * @param is_coarse
 */
void load_ensemble(EddySignature &sig_a, EddySignature &sig_b, const std::string &data_path, const bool is_coarse=true) {
    std::string level = "fine";
    if (is_coarse) {
        level = "coarse";
    }
    sig_a.load(data_path + "/signatures_A_" + level + ".mat");
    sig_b.load(data_path + "/signatures_B1_" + level + ".mat");
    EddySignature sig_bx = EddySignature();
    sig_bx.load(data_path + "/signatures_B2_" + level + ".mat");
    sig_b = sig_b + sig_bx;
    sig_bx.load(data_path + "/signatures_B3_" + level + ".mat");
    sig_b = sig_b + sig_bx;
    sig_bx.load(data_path + "/signatures_B4_" + level + ".mat");
    sig_b = sig_b + sig_bx;
    sig_b = sig_b / 4.0;
};


/** @brief Helper function to plot a trace digitised from a graph
 *
 * @param data
 * @param color
 * @param dash
 * @param name
 * @return
 */
ScatterPlot getTrace(const ArrayXXd &data, const std::string &color, const std::string &dash, std::string name="") {
    ScatterPlot p = ScatterPlot();
    p.x = data.row(0);
    p.y = data.row(1);
    p.name = name;
    p.setDash(dash);
    p.setColor(color);
    return p;
}
ScatterPlot getTrace(const ArrayXd &x, const ArrayXd &y, const std::string &color, const std::string &dash, std::string name="") {
    ScatterPlot p = ScatterPlot();
    p.x = x;
    p.y = y;
    p.name = name;
    p.setDash(dash);
    p.setColor(color);
    return p;
}


/** @brief Helper function to plot a trace of R13, calculated from parameters
 *
 * @param eta
 * @param pi_coles
 * @param shear_ratio
 * @param zeta
 * @param beta
 * @param color
 * @param dash
 * @param name
 * @return
 */
ScatterPlot getStress13Trace(
    const ArrayXd &eta,
    const double pi_coles, const double shear_ratio, const double zeta, const double beta,
    const std::string &color, const std::string &dash, const std::string &name
) {
    Eigen::ArrayXd r13;
    Eigen::ArrayXd r13_a;
    Eigen::ArrayXd r13_b;
    reynolds_stress_13(r13_a, r13_b, beta, eta, KAPPA_VON_KARMAN, pi_coles, shear_ratio, zeta);
    r13 = r13_a + r13_b;
    ScatterPlot p = ScatterPlot();
    p.x = eta;
    p.y = -1.0*r13;
    p.name = name;
    p.setDash(dash);
    p.setColor(color);
    return p;
}


TEST_F(AdemTest, test_save_load_signatures) {

    // Create a signature for test
    int n_lambda = 10;
    double dx = 0.05;
    EddySignature sig = EddySignature();
    sig.computeSignature("A", n_lambda, dx);

    // Save the signature
    std::string test_file_name = "test_signature_load_save.mat";
    sig.save(test_file_name);

    // Load from the file
    EddySignature sig2 = EddySignature();
    sig2.load(test_file_name);

    // Test for equality on all properties
    ASSERT_EQ(sig.eddy_type, sig2.eddy_type);
    ASSERT_TRUE(sig.lambda.isApprox(sig2.lambda));
    ASSERT_TRUE(sig.eta.isApprox(sig2.eta));
    ASSERT_TRUE(sig.domain_spacing.isApprox(sig2.domain_spacing));
    ASSERT_TRUE(sig.domain_extents.isApprox(sig2.domain_extents));
    // TODO Figure out how to assert equality of tensors
    // ASSERT_TRUE(sig.g.isApprox(sig2.g));
    ASSERT_TRUE(sig.j.isApprox(sig2.j));

};


TEST_F(AdemTest, test_get_type_a_eddy_intensity) {
    int n_lambda = 200;
    double dx = 0.01;
    EddySignature sig = EddySignature();
    sig.computeSignature("A", n_lambda, dx);
    plotEddyIntensities(sig, "A");
};


TEST_F(AdemTest, test_get_type_b1_eddy_intensity) {
    int n_lambda = 200;
    double dx = 0.01;
    EddySignature sig = EddySignature();
    sig.computeSignature("B1", n_lambda, dx);
    plotEddyIntensities(sig, "B1");
};


TEST_F(AdemTest, test_get_type_b2_eddy_intensity) {
    int n_lambda = 200;
    double dx = 0.01;
    EddySignature sig = EddySignature();
    sig.computeSignature("B2", n_lambda, dx);
    plotEddyIntensities(sig, "B2");
};


TEST_F(AdemTest, test_get_type_b3_eddy_intensity) {
    int n_lambda = 200;
    double dx = 0.01;
    EddySignature sig = EddySignature();
    sig.computeSignature("B3", n_lambda, dx);
    plotEddyIntensities(sig, "B3");
};


TEST_F(AdemTest, test_get_type_b4_eddy_intensity) {
    int n_lambda = 200;
    double dx = 0.01;
    EddySignature sig = EddySignature();
    sig.computeSignature("B4", n_lambda, dx);
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

        Eigen::Tensor<double, 2> slice_tens = psi.chip(j, 2);
        Eigen::ArrayXXd slice = tensor_to_array(slice_tens, dims[0], dims[1]);

        // Premultiply the spectrum with the wavenumber, which gives the spectral density
        slice = slice.transpose().eval();// * data.k1z;

        // Add to the vector of arrays (much easier to deal with than tensors for later plotting)
        spectra.push_back(slice);

        // Create a surface plot to show the spectrum term
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

        Layout lay = Layout();
        lay.xTitle("k1z");
        lay.yTitle("z");
        lay.zTitle("Psi");
        lay.xLog();

        Figure fig = Figure();
        fig.setLayout(lay);
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


TEST_F(AdemTest, test_reynolds_stress_13a_analytic) {

    // Digitally extracted from Perry and Marusic 1995, Figure 11 (eq. 51 for type a eddy shear stress):
    ArrayXXd fig_11_trace_1 = Eigen::ArrayXXd(2, 25);
    fig_11_trace_1 << 0.0, 0.012063814, 0.029733311, 0.046593536, 0.070664416, 0.093932194, 0.12360629, 0.156486629, 0.189363266, 0.224635643, 0.263913665, 0.301583013, 0.343258007, 0.384925599, 0.424186349, 0.467447811, 0.511511141, 0.557169574, 0.622848832, 0.677310555, 0.738973065, 0.797429331, 0.852680587, 0.921535246, 1.0,
                      1.0, 0.933830986, 0.864629672, 0.806207948, 0.741678448, 0.680224425, 0.617278925, 0.54973563, 0.486814803, 0.431616593, 0.373367579, 0.324351164, 0.272283946, 0.229461662, 0.192784164, 0.159219152, 0.124119486, 0.096736269, 0.063343967, 0.040651068, 0.024176974, 0.012300674, 0.003481346, 0.004011813, 0.001522318;

    Eigen::ArrayXd eta(25);
    eta = fig_11_trace_1.row(0);
    Eigen::ArrayXd r13_a(25);
    Eigen::ArrayXd r13_b(25);
    reynolds_stress_13(r13_a, r13_b, 0.0, eta, KAPPA_VON_KARMAN, 2.4, 20, 0.0); // params other than eta don't matter for r13a

    // Recreate P&M Figure 11
    Figure fig_11 = Figure();

    // Get the R13 profiles from our analytical function
    ScatterPlot r13a_validation_trace_1 = getTrace(fig_11_trace_1, "#1f77b4", "dash", "P&M 1995, Eqn. 51");

    ScatterPlot r13a_analytic_trace_1 = ScatterPlot();
    r13a_analytic_trace_1.x = eta;
    r13a_analytic_trace_1.y = -1.0 * r13_a;
    r13a_analytic_trace_1.name="es-flow, Eqn. 51 reproduction";
    r13a_analytic_trace_1.setDash("solid");

    // Traces from the paper
    fig_11.add(r13a_validation_trace_1);

    // Traces from es-flow
    fig_11.add(r13a_analytic_trace_1);

    Layout lay_11 = Layout();
    lay_11.xTitle("$z/\\delta_{c}$");
    lay_11.yTitle("$-\\overline{u_1u_3} / U_{\\tau}^2$");
    fig_11.setLayout(lay_11);
    fig_11.write("validation_perry_and_marusic_11.json");
}


TEST_F(AdemTest, test_sensitivity_to_eta_distribution) {
    /// The aim of this test is to determine how sensitive the Reynolds Stress calculation is to smoothness and
    /// monotonicity in the eta distribution

    // Digitally extracted from Perry and Marusic 1995, Figure 11 (eq. 51 for type a eddy shear stress):
    ArrayXXd fig_11_trace_1 = Eigen::ArrayXXd(2, 25);
    fig_11_trace_1 << 0.0, 0.012063814, 0.029733311, 0.046593536, 0.070664416, 0.093932194, 0.12360629, 0.156486629, 0.189363266, 0.224635643, 0.263913665, 0.301583013, 0.343258007, 0.384925599, 0.424186349, 0.467447811, 0.511511141, 0.557169574, 0.622848832, 0.677310555, 0.738973065, 0.797429331, 0.852680587, 0.921535246, 1.0,
        1.0, 0.933830986, 0.864629672, 0.806207948, 0.741678448, 0.680224425, 0.617278925, 0.54973563, 0.486814803, 0.431616593, 0.373367579, 0.324351164, 0.272283946, 0.229461662, 0.192784164, 0.159219152, 0.124119486, 0.096736269, 0.063343967, 0.040651068, 0.024176974, 0.012300674, 0.003481346, 0.004011813, 0.001522318;

    // Use the parameter set for a non-equilibrium case, to allow us to evaluate how sensitive we are to all the
    // different terms in R13. Params other than eta don't matter for r13a, but they do for r13b, so we choose the same
    // parameter set as in Figure 7, M&P199 Part 2.
    double pi_coles = 2.46;
    double shear_ratio = 34.5;
    double zeta = 8.01;
    double beta = 4.48;

    // Recreate P&M Figure 11 with different eta distributions...

    // From the paper
    Eigen::ArrayXd eta_paper(25);
    eta_paper = fig_11_trace_1.row(0);
    Eigen::ArrayXd r13_a_paper(25);
    r13_a_paper = fig_11_trace_1.row(1);

    // A coarse logarithmic distribution, with eta=0 prepended
    Eigen::ArrayXd eta_log_coarse(15);
    eta_log_coarse.setZero();
    eta_log_coarse.bottomRows(14) = Eigen::ArrayXd::LinSpaced(14, 5, 0).exp().inverse();
    std::cout << "ETA LOG COARSE: " << eta_log_coarse << std::endl;

    // A fine logarithmic distribution, with eta=0 prepended
    Eigen::ArrayXd eta_log_fine(200);
    eta_log_fine.setZero();
    eta_log_fine.bottomRows(199) = Eigen::ArrayXd::LinSpaced(199, 5, 0).exp().inverse();

    // A coarse linear distribution
    Eigen::ArrayXd eta_lin_coarse = Eigen::ArrayXd::LinSpaced(15, 0.0, 1.0);

    // A fine linear distribution
    Eigen::ArrayXd eta_lin_fine = Eigen::ArrayXd::LinSpaced(200, 0.0, 1.0);

    // Initialise Reynolds Stresses
    Eigen::ArrayXd r13_a_log_coarse(15);
    Eigen::ArrayXd r13_b_log_coarse(15);
    Eigen::ArrayXd r13_a_log_fine(200);
    Eigen::ArrayXd r13_b_log_fine(200);
    Eigen::ArrayXd r13_a_lin_coarse(15);
    Eigen::ArrayXd r13_b_lin_coarse(15);
    Eigen::ArrayXd r13_a_lin_fine(200);
    Eigen::ArrayXd r13_b_lin_fine(200);
    reynolds_stress_13(r13_a_log_coarse, r13_b_log_coarse, beta, eta_log_coarse, KAPPA_VON_KARMAN, pi_coles, shear_ratio, zeta);
    reynolds_stress_13(r13_a_log_fine,   r13_b_log_fine,   beta, eta_log_fine,   KAPPA_VON_KARMAN, pi_coles, shear_ratio, zeta);
    reynolds_stress_13(r13_a_lin_coarse, r13_b_lin_coarse, beta, eta_lin_coarse, KAPPA_VON_KARMAN, pi_coles, shear_ratio, zeta);
    reynolds_stress_13(r13_a_lin_fine,   r13_b_lin_fine,   beta, eta_lin_fine,   KAPPA_VON_KARMAN, pi_coles, shear_ratio, zeta);

    // Get scatter plot traces for each, plus the validation case
    ScatterPlot p_paper = getTrace(eta_paper, r13_a_paper, "#1f77b4", "solid", "P&M 1995, Eqn. 51");

    ScatterPlot p_log_coarse_a = getTrace(eta_log_coarse, -1.0*r13_a_log_coarse, "#2ca02c", "dash", "eta_log_coarse a");
    ScatterPlot p_log_coarse_b = getTrace(eta_log_coarse, -1.0*r13_b_log_coarse, "#2ca02c", "dot",  "eta_log_coarse b");

    ScatterPlot p_log_fine_a = getTrace(eta_log_fine, -1.0*r13_a_log_fine, "#d62728", "dash", "eta_log_fine a");
    ScatterPlot p_log_fine_b = getTrace(eta_log_fine, -1.0*r13_b_log_fine, "#d62728", "dot",  "eta_log_fine b");

    ScatterPlot p_lin_coarse_a = getTrace(eta_lin_coarse, -1.0*r13_a_lin_coarse, "#9467bd", "dash", "eta_lin_coarse a");
    ScatterPlot p_lin_coarse_b = getTrace(eta_lin_coarse, -1.0*r13_b_lin_coarse, "#9467bd", "dot",  "eta_lin_coarse b");

    ScatterPlot p_lin_fine_a = getTrace(eta_lin_fine, -1.0*r13_a_lin_fine, "#e377c2", "dash", "eta_lin_fine a");
    ScatterPlot p_lin_fine_b = getTrace(eta_lin_fine, -1.0*r13_b_lin_fine, "#e377c2", "dot",  "eta_lin_fine b");

    // Add Traces
    Figure fig = Figure();
    fig.add(p_paper);
    fig.add(p_log_coarse_a);
    fig.add(p_log_coarse_b);
    fig.add(p_log_fine_a);
    fig.add(p_log_fine_b);
    fig.add(p_lin_coarse_a);
    fig.add(p_lin_coarse_b);
    fig.add(p_lin_fine_a);
    fig.add(p_lin_fine_b);
    Layout lay = Layout();
    lay.xTitle("$z/\\delta_{c}$");
    lay.yTitle("$-\\overline{u_1u_3} / U_{\\tau}^2$");
    fig.setLayout(lay);
    fig.write("validation_sensitivity_to_eta_distribution.json");

}


TEST_F(AdemTest, test_reynolds_stress_13_analytic) {

    // Digitally extracted from Marusic and Perry 1995, Figure 4:
    ArrayXXd fig_4_trace_1 = Eigen::ArrayXXd(2, 28);
    fig_4_trace_1 << 0.000738437, 0.014580769, 0.02882141, 0.05082235, 0.070625881, 0.099317505, 0.128027031, 0.164496856, 0.212096937, 0.263049072, 0.325171742, 0.380706662, 0.435176442, 0.486343395, 0.526447224, 0.566555528, 0.603343105, 0.647935734, 0.681411533, 0.717107118, 0.751697286, 0.78738392, 0.818608606, 0.857629394, 0.89327575, 0.931110564, 0.965584372, 0.99778469,
                     1.004045738, 1.300138737, 1.553912596, 1.819642418, 2.055170176, 2.278439885, 2.47751572, 2.688547517, 2.857038645, 2.995227013, 3.036438497, 2.980995323, 2.865087605, 2.712949495, 2.512625031, 2.306252098, 2.087842646, 1.820904473, 1.578361565, 1.335778379, 1.087166864, 0.856680615, 0.656517264, 0.419922128, 0.243872094, 0.110121059, 0.01876972, 4.02784E-05;

    ArrayXXd fig_4_trace_2 = Eigen::ArrayXXd(2, 24);
    fig_4_trace_2 << 0.000733961, 0.016845309, 0.036684643, 0.06427085, 0.088549755, 0.136122983, 0.195954262, 0.235865649, 0.293557699, 0.342415359, 0.399078073, 0.443554342, 0.491378192, 0.537018058, 0.593797131, 0.641688111, 0.69961288, 0.743055338, 0.789827474, 0.836586185, 0.883326993, 0.928917631, 0.973358097, 0.998890107,
                     0.99799727, 1.239613775, 1.426753787, 1.643995167, 1.831054622, 2.035836559, 2.173863814, 2.233623487, 2.262818591, 2.231690125, 2.152032939, 2.042354942, 1.908422654, 1.726142899, 1.489225536, 1.264566224, 0.979240977, 0.766759158, 0.554216921, 0.35982009, 0.189617132, 0.073870527, 0.012580277, 0.006068607;

    ArrayXXd fig_4_trace_3 = Eigen::ArrayXXd(2, 20);
    fig_4_trace_3 << 0.000738437, 0.023518092, 0.065515004, 0.11755018, 0.166313858, 0.236192352, 0.2983553, 0.376123878, 0.460609994, 0.546259706, 0.601893083, 0.656421043, 0.723198102, 0.77883148, 0.77883148, 0.827778648, 0.866705453, 0.898946049, 0.92005639, 0.998890107,
                     1.004045738, 1.221347535, 1.462520978, 1.636980017, 1.732869386, 1.792085301, 1.77886057, 1.674626866, 1.491642239, 1.236055853, 1.047546376, 0.85300857, 0.603813019, 0.415303542, 0.415303542, 0.263205711, 0.15362841, 0.080462754, 0.049837768, 0.006068607;

    ArrayXXd fig_4_trace_4 = Eigen::ArrayXXd(2, 22);
    fig_4_trace_4 << 0.000733961, 0.015833874, 0.043514064, 0.06900132, 0.113352279, 0.153286044, 0.195466446, 0.236541431, 0.283188257, 0.354288527, 0.412079035, 0.469878494, 0.532144376, 0.624431068, 0.705605406, 0.757850925, 0.827868155, 0.862315111, 0.897867484, 0.998890107, 0.998890107, 0.996670322,
                     0.99799727, 1.106567612, 1.196791157, 1.25076417, 1.310443286, 1.339960617, 1.333146859, 1.320304773, 1.277119649, 1.185103717, 1.081232518, 0.965264383, 0.812924881, 0.587460001, 0.380341919, 0.270522947, 0.142236345, 0.087175815, 0.038143615, 0.006068607, 0.006068607, 0.006108886;

    ArrayXXd fig_4_trace_5 = Eigen::ArrayXXd(2, 21);
    fig_4_trace_5 << 0.000733961, 0.022520083, 0.056890958, 0.091284209, 0.134565553, 0.181198953, 0.238935756, 0.301138983, 0.353371076, 0.418930833, 0.47228177, 0.532296538, 0.584542057, 0.584542057, 0.64789098, 0.697903287, 0.745686858, 0.805688201, 0.854559288, 0.925623755, 0.998890107,
                     0.99799727, 1.070155967, 1.117919399, 1.135440489, 1.140703529, 1.11566381, 1.084374231, 1.016713285, 0.925039719, 0.821027546, 0.717236904, 0.607276959, 0.497457987, 0.497457987, 0.381389156, 0.289755868, 0.210259795, 0.118445255, 0.069171384, 0.025543199, 0.006068607;

    ArrayXXd fig_4_trace_6 = Eigen::ArrayXXd(2, 18);
    fig_4_trace_6 << 0.000733961, 0.030347513, 0.090308577, 0.090308577, 0.161381996, 0.206923404, 0.262458323, 0.322441764, 0.400214818, 0.463545839, 0.525766967, 0.556875294, 0.617982054, 0.681308599, 0.746850455, 0.811268992, 0.869014746, 0.997780214,
                     0.99799727, 0.991384905, 0.954006579, 0.954006579, 0.898281456, 0.849068003, 0.793624829, 0.726004162, 0.61572199, 0.523847032, 0.431992213, 0.389089038, 0.303302826, 0.217476336, 0.137658037, 0.076005281, 0.032618765, 0.006088747;

    // We want to ensure that the analytical function to get R13 works appropriately with nonuniformly spaced eta
    double lambda_max = log(1/0.01);
    double lambda_min = log(1/1);
    Eigen::ArrayXd lambda = Eigen::ArrayXd::LinSpaced(10000, lambda_min, lambda_max);
    Eigen::ArrayXd eta(10001);
    eta.setZero(10001);
    eta.topRows(10000) = lambda.exp().inverse();
    eta.reverseInPlace();

    // Get the R13 profiles from our analytical function
    ScatterPlot r13_analytic_trace_1 = getStress13Trace(eta, 3.23, 38.4, 15.32, 7.16, "#1f77b4", "solid", "es-flow, R_theta = 7257, zeta=15.32, beta=7.16");
    ScatterPlot r13_analytic_trace_2 = getStress13Trace(eta, 2.46, 34.5, 8.01,  4.48, "#ff7f0e", "solid", "es-flow, R_theta = 6395, zeta=8.01, beta=4.48");
    ScatterPlot r13_analytic_trace_3 = getStress13Trace(eta, 1.87, 31.5, 4.64,  2.90, "#2ca02c", "solid", "es-flow, R_theta = 5395, zeta=4.64, beta=2.90");
    ScatterPlot r13_analytic_trace_4 = getStress13Trace(eta, 1.19, 28.1, 2.18,  1.45, "#d62728", "solid", "es-flow, R_theta = 4155, zeta=2.18, beta=1.45");
    ScatterPlot r13_analytic_trace_5 = getStress13Trace(eta, 0.68, 25.4, 0.94,  0.65, "#9467bd", "solid", "es-flow, R_theta = 3153, zeta=0.94, beta=0.65");
    ScatterPlot r13_analytic_trace_6 = getStress13Trace(eta, 0.42, 23.6, 0.15,  0.0,  "#e377c2", "solid", "es-flow, R_theta = 2206, zeta=0.15, beta=0.00");

    // Get the same R13 profiles from the paper
    ScatterPlot r13_validation_trace_1 = getTrace(fig_4_trace_1, "#1f77b4", "dash", "M&P 1995, R_theta = 7257, zeta=15.32, beta=7.16");
    ScatterPlot r13_validation_trace_2 = getTrace(fig_4_trace_2, "#ff7f0e", "dash", "M&P 1995, R_theta = 6395, zeta=8.01, beta=4.48");
    ScatterPlot r13_validation_trace_3 = getTrace(fig_4_trace_3, "#2ca02c", "dash", "M&P 1995, R_theta = 5395, zeta=4.64, beta=2.90");
    ScatterPlot r13_validation_trace_4 = getTrace(fig_4_trace_4, "#d62728", "dash", "M&P 1995, R_theta = 4155, zeta=2.18, beta=1.45");
    ScatterPlot r13_validation_trace_5 = getTrace(fig_4_trace_5, "#9467bd", "dash", "M&P 1995, R_theta = 3153, zeta=0.94, beta=0.65");
    ScatterPlot r13_validation_trace_6 = getTrace(fig_4_trace_6, "#e377c2", "dash", "M&P 1995, R_theta = 2206, zeta=0.15, beta=0.00");

    // Recreate M&P Figure 4
    Figure fig_4 = Figure();

    // Traces from the paper
    fig_4.add(r13_validation_trace_1);
    fig_4.add(r13_validation_trace_2);
    fig_4.add(r13_validation_trace_3);
    fig_4.add(r13_validation_trace_4);
    fig_4.add(r13_validation_trace_5);
    fig_4.add(r13_validation_trace_6);

    // Traces from es-flow
    fig_4.add(r13_analytic_trace_1);
    fig_4.add(r13_analytic_trace_2);
    fig_4.add(r13_analytic_trace_3);
    fig_4.add(r13_analytic_trace_4);
    fig_4.add(r13_analytic_trace_5);
    fig_4.add(r13_analytic_trace_6);

    Layout lay_4 = Layout("Marusic & Perry Figure 4, 10APG flow case");
    lay_4.xTitle("$z/\\delta_{c}$");
    lay_4.yTitle("$-\\overline{u_1u_3} / U_{\\tau}^2$");
    fig_4.setLayout(lay_4);
    fig_4.write("validation_marusic_and_perry_4.json");

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


TEST_F(AdemTest, test_validate_reynolds_stress_13_skare_krogstad) {
    /* The Skare and Krogstad data is a very useful validation case, as it is taken for an equilibrium flow case.
     *
     * This means that the value of zeta is zero. Considering the formulation for f2 in the appendix of P&M 1995, that
     * can be very sensitive to errors in the coles wake function. This case with zeta = 0 still contains some error
     * (from the integration of the velocity deficit distribution, f) but is not highly sensitive to it.
     *
     * Also, the case includes a before-and-after distribution showing both the analytically derived and the
     * reconstructed (i.e. roundtripped through convolution-deconvolution cycle) R13 profiles, so we get to compare how
     * our cycle compares to that of P&M.
     *
     */

    // Digitally extracted from Marusic and Perry 1995, Figure 15:

    // This is the input profile (solid)
    ArrayXXd fig_15_trace_1 = Eigen::ArrayXXd(2, 21);
    fig_15_trace_1 << 0.000804653, 0.044885426, 0.07285776, 0.113558991, 0.17122183, 0.221257033, 0.278096286, 0.331564281, 0.39863919, 0.47340445, 0.538861771, 0.601803369, 0.646902974, 0.690305709, 0.735419514, 0.779684884, 0.82649556, 0.875003106, 0.924347255, 0.962603762, 1.000001183,
        1.086328594, 3.424933587, 4.956187839, 6.849204518, 9.465982712, 11.63729092, 13.64127869, 15.00469184, 16.06131928, 16.08709183, 15.22175876, 13.57655147, 11.93184118, 10.23146784, 8.252495311, 6.245691262, 4.322381772, 2.454735322, 0.893472254, 0.335303549, 3.55E-15;

    // This is the reconstructed profile (dotted)
    ArrayXXd fig_15_trace_2 = Eigen::ArrayXXd(2, 22);
    fig_15_trace_2 << 0.000804653, 0.049124053, 0.088968565, 0.133050522, 0.169507209, 0.216152221, 0.263657502, 0.308632859, 0.357869327, 0.419024122, 0.471715863, 0.524432454, 0.57292935, 0.618884491, 0.666553069, 0.698897744, 0.742307579, 0.798468793, 0.859714703, 0.925186224, 0.980437471, 1.000833052,
        1.086328594, 3.647656745, 5.707828207, 8.018578013, 9.828147464, 11.80456404, 13.53026027, 14.81034452, 15.78390339, 16.20002722, 15.83644249, 14.88789885, 13.27094908, 11.48693918, 9.368619725, 7.974961099, 6.107456646, 4.072466083, 2.34374057, 1.144145263, 0.529792861, 0.417804126;

    // Get the R13 profiles from the paper as traces for plotting
    ScatterPlot r13_validation_trace_1 = getTrace(fig_15_trace_1, "#1f77b4", "solid", "P&M 1995 Figure 15.   Skare & Krogstad");
    ScatterPlot r13_validation_trace_2 = getTrace(fig_15_trace_2, "#1f77b4", "dash", "P&M 1995 Figure 15.   Skare & Krogstad (reconstructed)");

    // Get the R13 profiles from our analytical function, using the skare and krogstad parameters
    Eigen::ArrayXd eta(10001);
    eta = Eigen::ArrayXd::LinSpaced(100000, 0, 1);
    ScatterPlot r13_analytic_trace_1 = getStress13Trace(eta, 6.85, 59.4, 0.0, 19.0, "#ff7f0e", "solid", "es-flow reproduction. Skare & Krogstad");

    // Run ADEM to do the roundtripping
    // Uncomment to make the signatures as part of the test
    // TODO Make only if not found, to avoid commenting/uncommenting
    bool is_coarse = false;
//    make_signatures(data_path, is_coarse);
    EddySignature signature_a = EddySignature();
    EddySignature signature_b = EddySignature();
    load_ensemble(signature_a, signature_b, data_path, is_coarse);
    ScatterPlot r11_trace_1 = ScatterPlot();
    ScatterPlot r22_trace_1 = ScatterPlot();
    ScatterPlot r33_trace_1 = ScatterPlot();
    ScatterPlot r13_trace_1 = ScatterPlot();
    getStressTraces(r11_trace_1, r22_trace_1, r33_trace_1, r13_trace_1, signature_a, signature_b, 6.85, 59.4, 0.0, 19.0, "dash", "es-flow reproduction. Skare & Krogstad (reconstructed)",  "#ff7f0e");

    // Negate to plot on the positive x axis as per the validation figure
    r13_trace_1.y = -1.0 * r13_trace_1.y;

    // Recreate Perry & Marusic 1995, Figure 15, with additional es-flow validation
    Figure fig_15 = Figure();
    fig_15.add(r13_validation_trace_1); // From the paper
    fig_15.add(r13_validation_trace_2); // From the paper
    fig_15.add(r13_analytic_trace_1);   // From es-flow
    fig_15.add(r13_trace_1);            // From es-flow
    Layout lay_15 = Layout("Validation against P&M 1995 Figure 15 (with R_theta = 49180)");
    lay_15.xTitle("$z/\\delta_{c}$");
    lay_15.yTitle("$-\\overline{u_1u_3} / U_{\\tau}^2$");
    fig_15.setLayout(lay_15);
    fig_15.write("validation_perry_and_marusic_15.json");

}


TEST_F(AdemTest, test_validate_signatures_perry_marusic) {

    // Load the signatures
    bool is_coarse = false;
    EddySignature signature_a = EddySignature();
    EddySignature signature_b = EddySignature();
    load_ensemble(signature_a, signature_b, data_path, is_coarse);

    ArrayXXd fig_20a_i11 = Eigen::ArrayXXd(2, 26);
    fig_20a_i11
        << 0, 0.012710542, 0.030316638, 0.05182384, 0.08898147, 0.1300581, 0.18483716, 0.2357049, 0.29636374, 0.35898846, 0.41573855, 0.4783607, 0.5585967, 0.6368718, 0.7073141, 0.76600707, 0.82468724, 0.87943816, 0.91658044, 0.9419913, 0.97133654, 1.0007021, 1.0516082, 1.1260028, 1.2239062, 1.50,
        1.6593906, 1.6418779, 1.611259, 1.5544281, 1.4713508, 1.3926289, 1.3051336, 1.2263831, 1.1476042, 1.08192, 1.0162529, 0.94620186, 0.8586324, 0.7667018, 0.6747941, 0.5829205, 0.469213, 0.33368298, 0.22440498, 0.1457287, 0.09760853, 0.084422655, 0.07117407, 0.040389933, 0.027004533, 0.0;

    ArrayXXd fig_20a_i13 = Eigen::ArrayXXd(2, 21);
    fig_20a_i13
        << 0, 0.023503713, 0.05093406, 0.08618198, 0.12728672, 0.17620902, 0.22121488, 0.28186864, 0.3620689, 0.45790172, 0.55568004, 0.64367235, 0.7238267, 0.80006206, 0.86259496, 0.90362066, 0.93487054, 0.9700394, 0.9993566, 1.0423712, 1.4980469,
        0, -0.061066702, -0.14832275, -0.22682244, -0.2878379, -0.340097, -0.38363394, -0.43149212, -0.4618262, -0.4702808, -0.46126255, -0.43917242, -0.39090434, -0.32954726, -0.24639612, -0.17204118, -0.102081485, -0.045210768, -0.014557855, -0.0030345472, -0.004372481;

    ArrayXXd fig_20a_i22 = Eigen::ArrayXXd(2, 23);
    fig_20a_i22
        << 0, 0.018947192, 0.03648182, 0.054011352, 0.0754471, 0.10665094, 0.17303112, 0.21991083, 0.28439146, 0.34693173, 0.41142765, 0.48179564, 0.57952243, 0.65573704, 0.7202125, 0.788579, 0.8452064, 0.9115866, 0.9701798, 1.0483195, 1.1304111, 1.2340457, 1.4980495,
        1.176464, 1.150218, 1.0934712, 1.0323671, 0.966883, 0.8926276, 0.79638124, 0.74817854, 0.6998736, 0.6646518, 0.6294187, 0.5985088, 0.5500107, 0.5016375, 0.4489753, 0.37886125, 0.3044581, 0.20821178, 0.14251183, 0.067983694, 0.028291034, 0.014617034, 0.0043686354;

    ArrayXXd fig_20a_i33 = Eigen::ArrayXXd(2, 32);
    fig_20a_i33
        << 0.0, 0.023483196, 0.041119818, 0.06269325, 0.09603633, 0.12936921, 0.17444114, 0.21754721, 0.25868744, 0.30175778, 0.35461667, 0.4133425, 0.46618098, 0.5346596, 0.608995, 0.6696255, 0.73806334, 0.7908406, 0.83576465, 0.8884807, 0.92553115, 0.95087314, 0.9898843, 1.0288852, 1.0737991, 1.1304673, 1.1793332, 1.2321156, 1.2849082, 1.3552966, 1.4374342, 1.5000051,
        0.0, 0.012935479, 0.043334138, 0.095496446, 0.1780915, 0.251972, 0.3301416, 0.39960802, 0.46037126, 0.49933675, 0.54696, 0.5945491, 0.6247433, 0.6504893, 0.6674866, 0.6714917, 0.66237944, 0.6402863, 0.59209496, 0.51771456, 0.42599595, 0.35613185, 0.26875913, 0.17267188, 0.115766004, 0.07622104, 0.05415063, 0.036414765, 0.027393447, 0.013912599, 0.013435401, 0.0043572737;

    ArrayXXd fig_20b_i11 = Eigen::ArrayXXd(2, 23);
    fig_20b_i11
        << 0, 0.07068844, 0.1666695, 0.27835685, 0.36659724, 0.4334307, 0.48672056, 0.5342177, 0.5837344, 0.62532634, 0.67665803, 0.74737716, 0.8137505, 0.8545449, 0.8894136, 0.93391985, 0.9706649, 1.00945, 1.0521617, 1.1184429, 1.2122818, 1.3199301, 1.5000205,
        0.016322415, 0.016322415, 0.018703382, 0.024518738, 0.03472676, 0.056264196, 0.091715895, 0.13412246, 0.18173045, 0.22154763, 0.25700384, 0.27592924, 0.25842202, 0.23056176, 0.19837682, 0.15315433, 0.11402357, 0.08182956, 0.050494146, 0.025177978, 0.011945607, 0.0073581133, 0.0017353186;

    ArrayXXd fig_20b_i13 = Eigen::ArrayXXd(2, 14);
    fig_20b_i13
        << 0, 0.2741542, 0.34076267, 0.42109883, 0.47993597, 0.56440324, 0.6233631, 0.7057494, 0.7918987, 0.8779049, 0.94425774, 1.0145372, 1.1240618, 1.5000103,
        0, -2.3333919E-4, -0.002682268, -0.0068347994, -0.014507807, -0.036872122, -0.054957043, -0.06691398, -0.06584696, -0.05263271, -0.033390157, -0.015006201, -0.0034729026, -8.676593E-4;

    ArrayXXd fig_20b_i22 = Eigen::ArrayXXd(2, 23);
    fig_20b_i22
        << 0, 0.14884579, 0.29180932, 0.4054613, 0.46828458, 0.515565, 0.5411826, 0.5746176, 0.6335705, 0.6963938, 0.74149364, 0.7689007, 0.8080032, 0.8490102, 0.89000964, 0.93290585, 0.9660264, 0.99715817, 1.0224689, 1.049784, 1.0810461, 1.1632366, 1.500023,
        0, 0.0023498295, 0.0038403252, 0.012338308, 0.030489888, 0.06258152, 0.08079781, 0.09726137, 0.12063707, 0.13878864, 0.14566672, 0.14474949, 0.13772498, 0.12461021, 0.110625885, 0.08968175, 0.07049375, 0.047830947, 0.031265218, 0.019913385, 0.012032943, 0.005803057, 0.0026086513;

    ArrayXXd fig_20b_i33 = Eigen::ArrayXXd(2, 20);
    fig_20b_i33
        << 0, 0.21936293, 0.31337926, 0.39961416, 0.44083738, 0.48011523, 0.5254081, 0.56480104, 0.64147526, 0.64147526, 0.70819265, 0.7767404, 0.8334498, 0.8939988, 0.9446548, 0.9952725, 1.0459669, 1.0967507, 1.1632637, 1.5000255,
        0, 0.002830162, 0.004290453, 0.009236455, 0.016043989, 0.023722952, 0.0409085, 0.05637945, 0.07693768, 0.07693768, 0.08626908, 0.08693706, 0.081578515, 0.07101401, 0.053551547, 0.033491757, 0.018626625, 0.009821928, 0.0053009023, 0.0017315528;

    // Plot the intensity functions as figures
    Figure fig_a = Figure();
    Figure fig_b = Figure();
    std::vector<std::string> labels = {"$I_{11}$", "$I_{12}$", "$I_{13}$", "$I_{22}$", "$I_{23}$", "$I_{33}$"};
    std::vector<std::string> dashes = {"solid", "na", "dot", "dash", "na", "longdashdot"};
    for (auto col = 0; col < 6; col++) {
        if (col == 0 || col == 2 || col == 3 || col == 5) {
            // Type A signatures into figure a
            ScatterPlot pa = ScatterPlot();
            pa.x = signature_a.eta;
            pa.y = signature_a.j.col(col) * 0.0440;
            pa.setDash(dashes[col]);
            pa.setColor("#ff7f0e");
            pa.name = labels[col];
            fig_a.add(pa);

            // Type b signatures into figure b
            ScatterPlot pb = ScatterPlot();
            pb.x = signature_a.eta;
            pb.y = signature_b.j.col(col) * 0.0342;
            pb.setDash(dashes[col]);
            pb.setColor("#ff7f0e");
            pb.name = labels[col];
            fig_b.add(pb);
        }
    }

    // Add digitised plots from the paper
    ScatterPlot pa_11 = getTrace(fig_20a_i11, "#1f77b4", "solid",       "$\\text{P&M 1995, }I_{11}$");
    ScatterPlot pa_13 = getTrace(fig_20a_i13, "#1f77b4", "dot",         "$\\text{P&M 1995, }I_{13}$");
    ScatterPlot pa_22 = getTrace(fig_20a_i22, "#1f77b4", "dash",        "$\\text{P&M 1995, }I_{22}$");
    ScatterPlot pa_33 = getTrace(fig_20a_i33, "#1f77b4", "longdashdot", "$\\text{P&M 1995, }I_{22}$");
    ScatterPlot pb_11 = getTrace(fig_20b_i11, "#1f77b4", "solid",       "$\\text{P&M 1995, }I_{11}$");
    ScatterPlot pb_13 = getTrace(fig_20b_i13, "#1f77b4", "dot",         "$\\text{P&M 1995, }I_{13}$");
    ScatterPlot pb_22 = getTrace(fig_20b_i22, "#1f77b4", "dash",        "$\\text{P&M 1995, }I_{22}$");
    ScatterPlot pb_33 = getTrace(fig_20b_i33, "#1f77b4", "longdashdot", "$\\text{P&M 1995, }I_{22}$");
    fig_a.add(pa_11);
    fig_a.add(pa_13);
    fig_a.add(pa_22);
    fig_a.add(pa_33);
    fig_b.add(pb_11);
    fig_b.add(pb_13);
    fig_b.add(pb_22);
    fig_b.add(pb_33);

    // Add axis labeling
    Layout lay_a = Layout("Type A Eddy Intensity");
    lay_a.xTitle("$\\eta$");
    lay_a.yTitle("$I_{ij}$");
    fig_a.setLayout(lay_a);
    Layout lay_b = Layout("Type B Eddy Intensity");
    lay_b.xTitle("$\\eta$");
    lay_b.yTitle("$I_{ij}$");
    fig_b.setLayout(lay_b);

    // Write figures
    fig_a.write("validation_perry_marusic_20_a.json");
    fig_b.write("validation_perry_marusic_20_b.json");

}


TEST_F(AdemTest, test_validate_stresses_perry_marusic) {

    // Uncomment to make the signatures as part of the test
    // TODO Make only if not found, to avoid commenting/uncommenting
    bool is_coarse = false;
//    make_signatures(data_path, is_coarse);

    // Load the eddy signatures
    EddySignature signature_a = EddySignature();
    EddySignature signature_b = EddySignature();
    load_ensemble(signature_a, signature_b, data_path, is_coarse);
    signature_a = signature_a * 0.044;
    signature_b = signature_b * 0.0342;

    // Run ADEM for all measurement locations in the "10APG" flow case in Marusic and Perry 1995, for reconstruction of their Figure 6
    ScatterPlot r11_trace_1 = ScatterPlot();
    ScatterPlot r22_trace_1 = ScatterPlot();
    ScatterPlot r33_trace_1 = ScatterPlot();
    ScatterPlot r13_trace_1 = ScatterPlot();
    getStressTraces(r11_trace_1, r22_trace_1, r33_trace_1, r13_trace_1, signature_a, signature_b, 3.23, 38.4, 15.32, 7.16, "solid", "$\\text{es-flow,} R_\\theta = 7257$", "#1f77b4");

    ScatterPlot r11_trace_2 = ScatterPlot();
    ScatterPlot r22_trace_2 = ScatterPlot();
    ScatterPlot r33_trace_2 = ScatterPlot();
    ScatterPlot r13_trace_2 = ScatterPlot();
    getStressTraces(r11_trace_2, r22_trace_2, r33_trace_2, r13_trace_2, signature_a, signature_b, 2.46, 34.5, 8.01, 4.48, "solid", "$\\text{es_flow,} R_\\theta = 6395$", "#ff7f0e");

    ScatterPlot r11_trace_3 = ScatterPlot();
    ScatterPlot r22_trace_3 = ScatterPlot();
    ScatterPlot r33_trace_3 = ScatterPlot();
    ScatterPlot r13_trace_3 = ScatterPlot();
    getStressTraces(r11_trace_3, r22_trace_3, r33_trace_3, r13_trace_3, signature_a, signature_b, 1.87, 31.5, 4.64, 2.90, "solid", "$\\text{es_flow,} R_\\theta = 5395$", "#2ca02c");

    ScatterPlot r11_trace_4 = ScatterPlot();
    ScatterPlot r22_trace_4 = ScatterPlot();
    ScatterPlot r33_trace_4 = ScatterPlot();
    ScatterPlot r13_trace_4 = ScatterPlot();
    getStressTraces(r11_trace_4, r22_trace_4, r33_trace_4, r13_trace_4, signature_a, signature_b, 1.19, 28.1, 2.18, 1.45, "solid", "$\\text{es_flow,} R_\\theta = 4155$", "#d62728");

    ScatterPlot r11_trace_5 = ScatterPlot();
    ScatterPlot r22_trace_5 = ScatterPlot();
    ScatterPlot r33_trace_5 = ScatterPlot();
    ScatterPlot r13_trace_5 = ScatterPlot();
    getStressTraces(r11_trace_5, r22_trace_5, r33_trace_5, r13_trace_5, signature_a, signature_b, 0.68, 25.4, 0.94, 0.65, "solid", "$\\text{es_flow,} R_\\theta = 3153$", "#9467bd");

    ScatterPlot r11_trace_6 = ScatterPlot();
    ScatterPlot r22_trace_6 = ScatterPlot();
    ScatterPlot r33_trace_6 = ScatterPlot();
    ScatterPlot r13_trace_6 = ScatterPlot();
    getStressTraces(r11_trace_6, r22_trace_6, r33_trace_6, r13_trace_6, signature_a, signature_b, 0.42, 23.6, 0.15, 0.0, "solid", "$\\text{es_flow,} R_\\theta = 2206$", "#e377c2");

    // Rerun ADEM for the Pi=2.46 location in the "10APG" flow case in Marusic and Perry 1995, for reconstruction of their Figure 7
    ScatterPlot r11_fig7_composite = ScatterPlot();
    ScatterPlot r22_fig7_composite = ScatterPlot();
    ScatterPlot r33_fig7_composite = ScatterPlot();
    ScatterPlot r13_fig7_composite = ScatterPlot();
    getStressTraces(r11_fig7_composite, r22_fig7_composite, r33_fig7_composite, r13_fig7_composite, signature_a, signature_b, 2.46, 34.5, 8.01, 4.48, "solid", "es-flow, composite", "#ff7f0e");
    r13_fig7_composite.y = -1.0*r13_fig7_composite.y;

    ScatterPlot r11_fig7_type_a = ScatterPlot();
    ScatterPlot r22_fig7_type_a = ScatterPlot();
    ScatterPlot r33_fig7_type_a = ScatterPlot();
    ScatterPlot r13_fig7_type_a = ScatterPlot();
    getStressTraces(r11_fig7_type_a, r22_fig7_type_a, r33_fig7_type_a, r13_fig7_type_a, signature_a, signature_b, 2.46, 34.5, 8.01, 4.48, "dash", "es-flow, Type A only", "#ff7f0e", "type_a");
    r13_fig7_type_a.y = -1.0*r13_fig7_type_a.y;

    ScatterPlot r11_fig7_type_b = ScatterPlot();
    ScatterPlot r22_fig7_type_b = ScatterPlot();
    ScatterPlot r33_fig7_type_b = ScatterPlot();
    ScatterPlot r13_fig7_type_b = ScatterPlot();
    getStressTraces(r11_fig7_type_b, r22_fig7_type_b, r33_fig7_type_b, r13_fig7_type_b, signature_a, signature_b, 2.46, 34.5, 8.01, 4.48, "dot", "es-flow, Type B only", "#ff7f0e", "type_b");
    r13_fig7_type_b.y = -1.0*r13_fig7_type_b.y;


    // Digitally extracted from Marusic and Perry 1995 Part 2, Figure 6a:
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


    // Digitally extracted from Marusic and Perry 1995 Part 2, Figure 6b:
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


    // Digitally extracted from Marusic and Perry 1995 Part 2, Figure 6c:
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


    // Digitally extracted from Marusic and Perry 1995 Part 2 Figure 7a:
    ArrayXXd fig_7a_trace_1 = Eigen::ArrayXXd(2, 20);
    fig_7a_trace_1 << 0.009985513, 0.015040939, 0.025460285, 0.038780015, 0.0525426, 0.087827416, 0.153750845, 0.242368482, 0.302693743, 0.369466454, 0.420463673, 0.507941934, 0.565361995, 0.621904302, 0.709029073, 0.789252447, 0.847588291, 0.899251072, 0.942719924, 0.999728212,
                      8.430250861, 8.076242007, 7.837038859, 7.831726513, 7.944171176, 8.286522381, 8.977127398, 9.436497787, 9.317412691, 8.907919331, 8.441170684, 7.159714707, 5.937432366, 4.773438269, 3.143876045, 1.863453025, 1.10673881, 0.640875553, 0.349581899, 0.174422036;


    ArrayXXd fig_7a_trace_2 = Eigen::ArrayXXd(2, 18);
    fig_7a_trace_2 << 0.010029037, 0.013761695, 0.020738304, 0.026524657, 0.034321259, 0.043897496, 0.069324284, 0.069324284, 0.094003308, 0.123064505, 0.163018943, 0.223631363, 0.306780216, 0.306780216, 0.387712311, 0.507481904, 0.656471592, 1,
                      5.639498278, 5.112247909, 4.467535662, 4.057452041, 3.705361535, 3.295277914, 2.591834727, 2.591834727, 2.181013281, 1.828775209, 1.418248893, 1.065420561, 0.712592228, 0.712592228, 0.47707821, 0.241121495, 0.063453025, 3.55E-15;

    ArrayXXd fig_7a_trace_3 = Eigen::ArrayXXd(2, 30);
    fig_7a_trace_3 << 0.009957505, 0.011863571, 0.015884154, 0.021018279, 0.021018279, 0.027486194, 0.039453415, 0.051580318, 0.074017763, 0.109966201, 0.152336418, 0.194541076, 0.251404954, 0.29952892, 0.356994144, 0.40143965, 0.441068108, 0.441068108, 0.496340411, 0.496340411, 0.53924427, 0.600094132, 0.668052548, 0.734931608, 0.798748823, 0.857708851, 0.909988507, 0.942719924, 0.999728212, 0.999728212,
                      2.732611904, 2.846679784, 3.075553369, 3.362715199, 3.362715199, 3.708165273, 4.284997541, 4.804869651, 5.556123955, 6.539498278, 7.407476636, 7.985784555, 8.447663551, 8.561731431, 8.443236596, 8.15105755, 7.742892277, 7.742892277, 6.985587801, 6.985587801, 6.286866699, 5.180865716, 3.842302017, 2.620167241, 1.688883424, 0.990309887, 0.524446631, 0.349581899, 0.174422036, 0.174422036;


    // Digitally extracted from Marusic and Perry 1995 Part 2 Figure 7b:
    ArrayXXd fig_7b_trace_1 = Eigen::ArrayXXd(2, 30);
    fig_7b_trace_1 << 0.010041006, 0.011780178, 0.013820786, 0.015616988, 0.018668711, 0.022212254, 0.027828982, 0.034700683, 0.045559599, 0.045559599, 0.066321695, 0.087884441, 0.115907277, 0.145855159, 0.191453291, 0.236445735, 0.292041284, 0.350746841, 0.40770645, 0.454336627, 0.506373463, 0.55390373, 0.603110226, 0.650551333, 0.69844786, 0.749870745, 0.824235936, 0.893085952, 0.953946793, 1.004603406,
                      4.513144026, 4.400145581, 4.27478222, 4.162187115, 4.073716833, 3.985296968, 3.896322509, 3.856858133, 3.866298828, 3.866298828, 3.911725085, 4.044714109, 4.214848298, 4.385486664, 4.59276602, 4.738876223, 4.798432013, 4.722276229, 4.498094372, 4.224906609, 3.828069683, 3.381974763, 2.861740762, 2.341607597, 1.809159933, 1.276712268, 0.706917767, 0.347478096, 0.173663422, 0.086554414;
    
    ArrayXXd fig_7b_trace_2 = Eigen::ArrayXXd(2, 21);
    fig_7b_trace_2 << 0.01000225, 0.011196887, 0.012831892, 0.015198018, 0.017665334, 0.021929294, 0.027869799, 0.035089575, 0.045226948, 0.060528242, 0.081006308, 0.109425699, 0.138407062, 0.180067793, 0.233164521, 0.299092186, 0.364325674, 0.467312253, 0.593800739, 0.743942051, 1.009433058,
                      3.808394215, 3.671170114, 3.533693924, 3.33404023, 3.146953123, 2.909700504, 2.647465964, 2.360602427, 2.110581551, 1.798332753, 1.486083955, 1.247923819, 1.047564279, 0.834537316, 0.646290604, 0.470509643, 0.357107857, 0.230786562, 0.116931018, 0.040321475, 5.04176E-05;

    ArrayXXd fig_7b_trace_3 = Eigen::ArrayXXd(2, 25);
    fig_7b_trace_3 << 0.010038748, 0.012634005, 0.015457921, 0.020292237, 0.025296678, 0.033204117, 0.042973209, 0.059666194, 0.075422785, 0.093565435, 0.12221108, 0.15887917, 0.202715614, 0.261101985, 0.330134976, 0.409779474, 0.483151488, 0.548783293, 0.614720622, 0.675751523, 0.728896072, 0.778827818, 0.847978715, 0.918718682, 0.981283589,
                      0.704800228, 0.776519263, 0.885635553, 1.018725412, 1.189464612, 1.421473801, 1.665999159, 2.021052499, 2.340019442, 2.671552972, 3.089086323, 3.506670091, 3.88736078, 4.21849097, 4.327254336, 4.176556131, 3.828573859, 3.295521184, 2.614240766, 1.957891851, 1.450123602, 0.991915853, 0.508776601, 0.235891344, 0.099171418;


    // Digitally extracted from Marusic and Perry 1995 Part 2 Figure 7c:
    ArrayXXd fig_7c_trace_1 = Eigen::ArrayXXd(2, 30);
    fig_7c_trace_1 << 0.00204232, 0.012115212, 0.028305497, 0.04551455, 0.069912376, 0.094329335, 0.12285963, 0.152442174, 0.18715204, 0.232116551, 0.281223096, 0.330382253, 0.380603226, 0.380603226, 0.421645717, 0.463730892, 0.504835562, 0.540817693, 0.576818956, 0.608716448, 0.645750827, 0.677657885, 0.710588494, 0.740405403, 0.780491305, 0.818482274, 0.856463678, 0.904671029, 0.949759896, 1.000985288,
                      1.000019132, 1.168111117, 1.350202797, 1.536957853, 1.700281237, 1.844912855, 1.970814441, 2.068668809, 2.157129465, 2.226802694, 2.24970824, 2.221211426, 2.15532151, 2.15532151, 2.056807094, 1.939591345, 1.780328684, 1.625786795, 1.452553139, 1.288703629, 1.106114523, 0.93291913, 0.759714171, 0.628594387, 0.46466835, 0.347490865, 0.239659263, 0.141077886, 0.089254625, 0.042046911;

    ArrayXXd fig_7c_trace_2 = Eigen::ArrayXXd(2, 25);
    fig_7c_trace_2 << 0.001018768, 0.009278922, 0.032892345, 0.059581205, 0.085232164, 0.122170885, 0.15498192, 0.189844841, 0.225740879, 0.26263177, 0.307734986, 0.34667298, 0.394827718, 0.436836366, 0.488061757, 0.541319903, 0.595596816, 0.652949167, 0.706192963, 0.761483862, 0.82393962, 0.872056094, 0.925280759, 0.96928868, 1.001013985,
                      1.000009566, 0.939165662, 0.868851518, 0.793835734, 0.732848342, 0.643718074, 0.587336662, 0.526263177, 0.455834242, 0.413433393, 0.347591307, 0.305171325, 0.257992309, 0.21554363, 0.168335916, 0.135127896, 0.106583252, 0.073336968, 0.054147774, 0.034939448, 0.01566416, 0.005868679, 0.005371253, 0.009632861, 0.01400926;

    ArrayXXd fig_7c_trace_3 = Eigen::ArrayXXd(2, 29);
    fig_7c_trace_3 << 4.78295E-06, 0.008006658, 0.018055635, 0.033241501, 0.045361495, 0.063608449, 0.079813082, 0.107243299, 0.134716562, 0.164236928, 0.192786355, 0.229509843, 0.274450439, 0.320471981, 0.37165911, 0.424950736, 0.477271423, 0.528630737, 0.571796859, 0.620114217, 0.658181714, 0.697263196, 0.742500335, 0.788732327, 0.844157149, 0.892369282, 0.942595037, 0.982556582, 1.003027607,
                      0.009345884, 0.182168207, 0.364317282, 0.52772676, 0.686491993, 0.859218657, 1.027291511, 1.227969733, 1.386591479, 1.545194092, 1.652403911, 1.773556027, 1.866593966, 1.903547036, 1.893722857, 1.827804243, 1.710492835, 1.532442748, 1.35914213, 1.153083089, 0.961138533, 0.778530295, 0.581845836, 0.413189462, 0.263138763, 0.159884444, 0.089321586, 0.046892039, 0.046700721;


    // Digitally extracted from Marusic and Perry 1995 Part 2 Figure 7d:
    ArrayXXd fig_7d_trace_1 = Eigen::ArrayXXd(2, 22);
    fig_7d_trace_1 << 1.90243E-05, 0.015894783, 0.04092756, 0.072140057, 0.111483846, 0.153988836, 0.194573959, 0.25154057, 0.31370715, 0.379050784, 0.440419929, 0.497838367, 0.553225962, 0.607625881, 0.658999355, 0.706216024, 0.759588632, 0.806702253, 0.860959489, 0.918227318, 0.966233494, 0.999889025,
                      1.018662815, 1.407217811, 1.850063335, 2.230604354, 2.634347063, 2.936951962, 3.122945973, 3.238700769, 3.253286047, 3.15118117, 2.948049457, 2.620562389, 2.285330856, 1.919010484, 1.52163298, 1.202080622, 0.843552276, 0.625090167, 0.39874091, 0.21900113, 0.124941936, 0.108866423;
    
    ArrayXXd fig_7d_trace_2 = Eigen::ArrayXXd(2, 17);
    fig_7d_trace_2 << 1.10975E-05, 0.044913146, 0.104132548, 0.159266487, 0.216423341, 0.287891207, 0.35630092, 0.421620773, 0.494084243, 0.554291322, 0.626746865, 0.702236779, 0.766481761, 0.82257167, 0.884793737, 0.93780964, 0.999992073,
                      1.010886642, 0.940203528, 0.845969945, 0.759575949, 0.688702592, 0.57872642, 0.468797808, 0.39002145, 0.303357943, 0.240213199, 0.161325865, 0.10571949, 0.081392197, 0.057191731, 0.017343799, 0.00874324, 0.007776173;
    
    ArrayXXd fig_7d_trace_3 = Eigen::ArrayXXd(2, 28);
    fig_7d_trace_3 << 0.001035238, 0.009940185, 0.019864515, 0.032870779, 0.05191408, 0.070989088, 0.092158352, 0.117357591, 0.148720697, 0.185172796, 0.22776498, 0.279618817, 0.334610074, 0.387673537, 0.444909658, 0.500194206, 0.54320651, 0.586306009, 0.625327971, 0.669423073, 0.710452096, 0.760703137, 0.810874911, 0.852812343, 0.895761233, 0.943775336, 0.982567421, 0.999904879,
                      0.0155682, 0.248679002, 0.51291035, 0.753765618, 1.072287495, 1.35970468, 1.592656946, 1.872202837, 2.104996568, 2.345487204, 2.562554199, 2.693940609, 2.747517728, 2.692260132, 2.543625044, 2.30948376, 2.114413585, 1.833805505, 1.55326084, 1.295965427, 1.046493747, 0.750222346, 0.531712676, 0.391091565, 0.258230774, 0.156395407, 0.10135976, 0.093314077;


    // Recreate Figure 6a
    ScatterPlot trace1_6a = getTrace(fig_6a_trace_1, "#1f77b4", "dash", "$\\text{M&P 1995,} R_\\theta=7257$");
    ScatterPlot trace2_6a = getTrace(fig_6a_trace_2, "#ff7f0e", "dash", "$\\text{M&P 1995,} R_\\theta=6395$");
    ScatterPlot trace3_6a = getTrace(fig_6a_trace_3, "#2ca02c", "dash", "$\\text{M&P 1995,} R_\\theta=5395$");
    ScatterPlot trace4_6a = getTrace(fig_6a_trace_4, "#d62728", "dash", "$\\text{M&P 1995,} R_\\theta=4155$");
    ScatterPlot trace5_6a = getTrace(fig_6a_trace_5, "#9467bd", "dash", "$\\text{M&P 1995,} R_\\theta=3153$");
    ScatterPlot trace6_6a = getTrace(fig_6a_trace_6, "#e377c2", "dash", "$\\text{M&P 1995,} R_\\theta=2206$");

    // Add traces from the paper and es-flow to the figure
    Figure fig_6a = Figure();
    fig_6a.add(trace1_6a);
    fig_6a.add(trace2_6a);
    fig_6a.add(trace3_6a);
    fig_6a.add(trace4_6a);
    fig_6a.add(trace5_6a);
    fig_6a.add(trace6_6a);
    fig_6a.add(r11_trace_1);
    fig_6a.add(r11_trace_2);
    fig_6a.add(r11_trace_3);
    fig_6a.add(r11_trace_4);
    fig_6a.add(r11_trace_5);
    fig_6a.add(r11_trace_6);
    Layout lay_6a = Layout();
    lay_6a.xLog();
    lay_6a.xTitle("$z/\\delta_{c}$");
    lay_6a.yTitle("$\\overline{u_1^2} / U_{\\tau}^2$");
    fig_6a.setLayout(lay_6a);
    fig_6a.write("validation_perry_marusic_6a.json");

    // Recreate Figure 6b
    ScatterPlot trace1_6b = getTrace(fig_6b_trace_1, "#1f77b4", "dash", "$\\text{M&P 1995,} R_\\theta=7257$");
    ScatterPlot trace2_6b = getTrace(fig_6b_trace_2, "#ff7f0e", "dash", "$\\text{M&P 1995,} R_\\theta=6395$");
    ScatterPlot trace3_6b = getTrace(fig_6b_trace_3, "#2ca02c", "dash", "$\\text{M&P 1995,} R_\\theta=5395$");
    ScatterPlot trace4_6b = getTrace(fig_6b_trace_4, "#d62728", "dash", "$\\text{M&P 1995,} R_\\theta=4155$");
    ScatterPlot trace5_6b = getTrace(fig_6b_trace_5, "#9467bd", "dash", "$\\text{M&P 1995,} R_\\theta=3153$");
    ScatterPlot trace6_6b = getTrace(fig_6b_trace_6, "#e377c2", "dash", "$\\text{M&P 1995,} R_\\theta=2206$");

    // Add traces from the paper and es-flow to the figure
    Figure fig_6b = Figure();
    fig_6b.add(trace1_6b);
    fig_6b.add(trace2_6b);
    fig_6b.add(trace3_6b);
    fig_6b.add(trace4_6b);
    fig_6b.add(trace5_6b);
    fig_6b.add(trace6_6b);
    fig_6b.add(r22_trace_1);
    fig_6b.add(r22_trace_2);
    fig_6b.add(r22_trace_3);
    fig_6b.add(r22_trace_4);
    fig_6b.add(r22_trace_5);
    fig_6b.add(r22_trace_6);
    Layout lay_6b = Layout();
    lay_6b.xLog();
    lay_6b.xTitle("$z/\\delta_{c}$");
    lay_6b.yTitle("$\\overline{u_2^2} / U_{\\tau}^2$");
    fig_6b.setLayout(lay_6b);
    fig_6b.write("validation_perry_marusic_6b.json");

    // Recreate Figure 6c
    ScatterPlot trace1_6c = getTrace(fig_6c_trace_1, "#1f77b4", "dash", "$\\text{M&P 1995,} R_\\theta=7257$");
    ScatterPlot trace2_6c = getTrace(fig_6c_trace_2, "#ff7f0e", "dash", "$\\text{M&P 1995,} R_\\theta=6395$");
    ScatterPlot trace3_6c = getTrace(fig_6c_trace_3, "#2ca02c", "dash", "$\\text{M&P 1995,} R_\\theta=5395$");
    ScatterPlot trace4_6c = getTrace(fig_6c_trace_4, "#d62728", "dash", "$\\text{M&P 1995,} R_\\theta=4155$");
    ScatterPlot trace5_6c = getTrace(fig_6c_trace_5, "#9467bd", "dash", "$\\text{M&P 1995,} R_\\theta=3153$");
    ScatterPlot trace6_6c = getTrace(fig_6c_trace_6, "#e377c2", "dash", "$\\text{M&P 1995,} R_\\theta=2206$");

    // Add traces from the paper and es-flow to the figure
    Figure fig_6c = Figure();
    fig_6c.add(trace1_6c);
    fig_6c.add(trace2_6c);
    fig_6c.add(trace3_6c);
    fig_6c.add(trace4_6c);
    fig_6c.add(trace5_6c);
    fig_6c.add(trace6_6c);
    fig_6c.add(r33_trace_1);
    fig_6c.add(r33_trace_2);
    fig_6c.add(r33_trace_3);
    fig_6c.add(r33_trace_4);
    fig_6c.add(r33_trace_5);
    fig_6c.add(r33_trace_6);
    Layout lay_6c = Layout();
    lay_6c.xTitle("$z/\\delta_{c}$");
    lay_6c.yTitle("$\\overline{u_3^2} / U_{\\tau}^2$");
    fig_6c.setLayout(lay_6c);
    fig_6c.write("validation_perry_marusic_6c.json");

    // Recreate P&M Figure 7a
    ScatterPlot trace1_7a = getTrace(fig_7a_trace_1, "#1f77b4", "solid", "M&P 1995, composite");
    ScatterPlot trace2_7a = getTrace(fig_7a_trace_2, "#1f77b4", "dash", "M&P 1995, Type A only");
    ScatterPlot trace3_7a = getTrace(fig_7a_trace_3, "#1f77b4", "dot", "M&P 1995, Type B only");

    // Add traces from the paper and es-flow to the figure
    Figure fig_7a = Figure();
    fig_7a.add(trace1_7a);
    fig_7a.add(trace2_7a);
    fig_7a.add(trace3_7a);
    fig_7a.add(r11_fig7_composite);
    fig_7a.add(r11_fig7_type_a);
    fig_7a.add(r11_fig7_type_b);
    Layout lay_7a = Layout();
    lay_7a.xTitle("$z/\\delta_{c}$");
    lay_7a.yTitle("$\\overline{u_1^2} / U_{\\tau}^2$");
    lay_7a.xLog();
    fig_7a.setLayout(lay_7a);
    fig_7a.write("validation_perry_marusic_7a.json");

    // Recreate P&M Figure 7b
    ScatterPlot trace1_7b = getTrace(fig_7b_trace_1, "#1f77b4", "solid", "M&P 1995, composite");
    ScatterPlot trace2_7b = getTrace(fig_7b_trace_2, "#1f77b4", "dash", "M&P 1995, Type A only");
    ScatterPlot trace3_7b = getTrace(fig_7b_trace_3, "#1f77b4", "dot", "M&P 1995, Type B only");

    // Add traces from the paper and es-flow to the figure
    Figure fig_7b = Figure();
    fig_7b.add(trace1_7b);
    fig_7b.add(trace2_7b);
    fig_7b.add(trace3_7b);
    fig_7b.add(r22_fig7_composite);
    fig_7b.add(r22_fig7_type_a);
    fig_7b.add(r22_fig7_type_b);
    Layout lay_7b = Layout();
    lay_7b.xTitle("$z/\\delta_{c}$");
    lay_7b.yTitle("$\\overline{u_2^2} / U_{\\tau}^2$");
    lay_7b.xLog();
    fig_7b.setLayout(lay_7b);
    fig_7b.write("validation_perry_marusic_7b.json");

    // Recreate P&M Figure 7c
    ScatterPlot trace1_7c = getTrace(fig_7c_trace_1, "#1f77b4", "solid", "M&P 1995, composite");
    ScatterPlot trace2_7c = getTrace(fig_7c_trace_2, "#1f77b4", "dash", "M&P 1995, Type A only");
    ScatterPlot trace3_7c = getTrace(fig_7c_trace_3, "#1f77b4", "dot", "M&P 1995, Type B only");

    // Add traces from the paper and es-flow to the figure
    Figure fig_7c = Figure();
    fig_7c.add(trace1_7c);
    fig_7c.add(trace2_7c);
    fig_7c.add(trace3_7c);
    fig_7c.add(r13_fig7_composite);
    fig_7c.add(r13_fig7_type_a);
    fig_7c.add(r13_fig7_type_b);
    Layout lay_7c = Layout();
    lay_7c.xTitle("$z/\\delta_{c}$");
    lay_7c.yTitle("$-\\overline{u_1u_3} / U_{\\tau}^2$");
    fig_7c.setLayout(lay_7c);
    fig_7c.write("validation_perry_marusic_7c.json");

    // Recreate P&M Figure 7d
    ScatterPlot trace1_7d = getTrace(fig_7d_trace_1, "#1f77b4", "solid", "M&P 1995, composite");
    ScatterPlot trace2_7d = getTrace(fig_7d_trace_2, "#1f77b4", "dash", "M&P 1995, Type A only");
    ScatterPlot trace3_7d = getTrace(fig_7d_trace_3, "#1f77b4", "dot", "M&P 1995, Type B only");

    // Add traces from the paper and es-flow to the figure
    Figure fig_7d = Figure();
    fig_7d.add(trace1_7d);
    fig_7d.add(trace2_7d);
    fig_7d.add(trace3_7d);
    fig_7d.add(r33_fig7_composite);
    fig_7d.add(r33_fig7_type_a);
    fig_7d.add(r33_fig7_type_b);
    Layout lay_7d = Layout();
    lay_7d.xTitle("$z/\\delta_{c}$");
    lay_7d.yTitle("$\\overline{u_3^2} / U_{\\tau}^2$");
    fig_7d.setLayout(lay_7d);
    fig_7d.write("validation_perry_marusic_7d.json");
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

    // Load the eddy signatures
    bool is_coarse = false;
    EddySignature signature_a = EddySignature();
    EddySignature signature_b = EddySignature();
    load_ensemble(signature_a, signature_b, data_path, is_coarse);

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
