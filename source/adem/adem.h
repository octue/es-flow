/*
 * adem.h Implementation of the Attached-Detached Eddy Method
 *
 * Author:              Tom Clark  (thclark @ github)
 *
 * Copyright (c) 2019 Octue Ltd. All Rights Reserved.
 *
 */

#ifndef SOURCE_ADEM_H_
#define SOURCE_ADEM_H_

#include <boost/algorithm/string/classification.hpp> // Include boost::for is_any_of
#include <boost/algorithm/string/split.hpp> // Include for boost::split
#include <Eigen/Dense>
#include <Eigen/Core>
#include <stdexcept>
#include <unsupported/Eigen/CXX11/Tensor>
#include <unsupported/Eigen/FFT>

#include "cpplot.h"

#include "adem/signature.h"
#include "profile.h"
#include "variable_readers.h"
#include "relations/stress.h"
#include "relations/velocity.h"
#include "utilities/conv.h"
#include "utilities/filter.h"
#include "utilities/tensors.h"

using namespace utilities;
using namespace cpplot;

namespace es {


/** @brief Data container for ADEM input parameters and results.
 *
 * # Upcoming refactor
 *
 * The AdemData class is treated as a struct, whose contents are updated by a number of functions... But, the
 * functions operating on this should be refactored into class methods. That way, appropriate validation of
 * the contents can be undertaken prior to application of each function.
 *
 * This refactor is captured in [issue #34](https://github.com/octue/es-flow/issues/34).
 *
 */
class AdemData {
public:

    /// Eddy types used to create the results.
    std::vector<std::string> eddy_types = {"A", "B1+B2+B3+B4"};

    //// Clauser parameter @f$ \Beta @f$, representing acceleration/decelaration of the boundary layer.
    double beta;

    /// Atmospheric boundary layer thickness @f$ \delta_c @f$ [m].
    double delta_c;

    /// von Karman constant @f$ \kappa @f$. Typically @f$ \kappa = 0.41 @f$.
    double kappa;

    /// Coles wake parameter @f$ \Pi @f$.
    double pi_coles;

    /// Ratio between free-stream and skin friction velocities @f$ S = U_{inf}/U_{\tau} @f$.
    double shear_ratio;

    /// Free-stream velocity @f$ U_{inf}|_{z = \delta_c} @f$) [m/s].
    double u_inf;

    /// Skin friction velocity [m/s].
    double u_tau;

    /// Scaled streamwise derivative @f$ \zeta @f$ of the Coles wake parameter @f$ \Pi @f$.
    double zeta;

    /// Vertical coordinates used in the analysis [m].
    Eigen::VectorXd z;

    /// Nondimensionalised vertical coordinates used in the analysis @f$ \eta = z/\delta_{c} @f$
    Eigen::VectorXd eta;

    /// Parameterised nondimensional vertical coordinates used in the analysis
    Eigen::VectorXd lambda_e;

    /// Horizontal mean velocity varying with vertical coordinate [m/s]
    Eigen::VectorXd u_horizontal;

    /// Reynolds stress profiles from all eddy types
    Eigen::ArrayXXd reynolds_stress;

    /// Reynolds Stress profiles, contributions from Type A eddies only
    Eigen::ArrayXXd reynolds_stress_a;

    /// Reynolds Stress profiles, contributions from Type B eddies only
    Eigen::ArrayXXd reynolds_stress_b;

    /// Wavenumbers @f$ k_{1}z @f$ for which spectra are defined at each vertical coordinate
    Eigen::ArrayXXd k1z;

    /// Turbulent Spectra @f$ \Psi @f$ (corresponding to wavenumber `k1z`) at each vertical coordinate
    Eigen::Tensor<double, 3> psi;

    /// Turbulent Spectra (corresponding to wavenumber `k1z`) at each vertical coordinate from Type A eddies only
    Eigen::Tensor<double, 3> psi_a;

    /// Turbulent Spectra (corresponding to wavenumber `k1z`) at each vertical coordinate from Type B eddies only
    Eigen::Tensor<double, 3> psi_b;

    /// Negated convolution function @f$ -(T^2)\omega @f$ encapsulating variation of eddy strength and scale for Type A eddies
    Eigen::VectorXd t2wa;

    /// Negated convolution function @f$ -(T^2)\omega @f$ encapsulating variation of eddy strength and scale for Type B eddies
    Eigen::VectorXd t2wb;

    /// Fit residuals from the deconvolution of `t2wa`
    Eigen::VectorXd residual_a;

    /// Fit residuals from the deconvolution of `t2wb`
    Eigen::VectorXd residual_b;

    /** @brief Load data from a *.mat file containing eddy signature data.
     *
     * @param[in] file_name File name (including relative or absolute path)
     * @param[in] print_var Boolean, default true. Print variables as they are read in (not advised except for debugging!)
     */
    void load(std::string file_name, bool print_var = true) {
        std::cout << "Reading adem data from file " << file_name << std::endl;

        // Open the MAT file for reading
        mat_t *matfp = Mat_Open(file_name.c_str(), MAT_ACC_RDONLY);
        if (matfp == NULL) {
            std::string msg = "Error reading MAT file: ";
            throw std::invalid_argument(msg + file_name);
        }

        // Use the variable readers to assist
        k1z                 = readArrayXXd(matfp, "k1z", print_var);
        beta                = readDouble(matfp, "beta", print_var);
        delta_c             = readDouble(matfp, "delta_c", print_var);
        kappa               = readDouble(matfp, "kappa", print_var);
        pi_coles            = readDouble(matfp, "pi_coles", print_var);
        shear_ratio         = readDouble(matfp, "shear_ratio", print_var);
        u_inf               = readDouble(matfp, "u_inf", print_var);
        u_tau               = readDouble(matfp, "u_tau", print_var);
        zeta                = readDouble(matfp, "zeta", print_var);
        z                   = readVectorXd(matfp, "z", print_var);
        eta                 = readVectorXd(matfp, "eta", print_var);
        lambda_e            = readVectorXd(matfp, "lambda_e", print_var);
        u_horizontal        = readVectorXd(matfp, "u_horizontal", print_var);
        reynolds_stress     = readArrayXXd(matfp, "reynolds_stress", print_var);
        reynolds_stress_a   = readArrayXXd(matfp, "reynolds_stress_a", print_var);
        reynolds_stress_b   = readArrayXXd(matfp, "reynolds_stress_b", print_var);
        k1z                 = readArrayXXd(matfp, "k1z", print_var);
        psi                 = readTensor3d(matfp, "psi", print_var);
        psi_a               = readTensor3d(matfp, "psi_a", print_var);
        psi_b               = readTensor3d(matfp, "psi_b", print_var);
        t2wa                = readVectorXd(matfp, "t2wa", print_var);
        t2wb                = readVectorXd(matfp, "t2wb", print_var);
        residual_a          = readVectorXd(matfp, "residual_a", print_var);
        residual_b          = readVectorXd(matfp, "residual_b", print_var);

        // Special handling to split the string of types into a vector of strings
        std::string typestr = readString(matfp, "eddy_types", print_var);
        boost::split(eddy_types, typestr, boost::is_any_of(", "), boost::token_compress_on);

        // Close the file
        Mat_Close(matfp);
        std::cout << "Finished reading adem data" << std::endl;
    }

    /** @brief Save eddy signature data to a *.mat file.
     *
     * @param[in] filename File name (including relative or absolute path)
     */
    void save(std::string filename) {
        std::cout << "Writing signature data..." << std::endl;
        throw std::invalid_argument("Error writing mat file - function not implemented");
    }

};


/** @brief Print information about the AdemData attributes to ostream using the << operator.
 *
 * @param os The ostream to print to
 * @param data The AdemDataclass instance
 * @return os The same ostream, with the data representation added to the stream
 */
std::ostream &operator<<(std::ostream &os, AdemData const &data) {
    os << "AdemData() with fields:" << std::endl;
    os << "    eddy_types:        ";
    for (std::vector<std::string>::const_iterator i = data.eddy_types.begin(); i != data.eddy_types.end(); ++i) {
        os << *i << ", ";
    }
    os << std::endl;
    os << "    beta:              " << data.beta << std::endl
       << "    delta_c:           " << data.delta_c << std::endl
       << "    kappa:             " << data.kappa << std::endl
       << "    pi_coles:          " << data.pi_coles << std::endl
       << "    shear_ratio:       " << data.shear_ratio << std::endl
       << "    u_inf:             " << data.u_inf << std::endl
       << "    u_tau:             " << data.u_tau << std::endl
       << "    zeta:              " << data.zeta << std::endl
       << "    z:                 [" << data.z.size() << " x 1]" << std::endl
       << "    eta:               [" << data.eta.size() << " x 1]" << std::endl
       << "    lambda_e:          [" << data.lambda_e.size() << " x 1]" << std::endl
       << "    u_horizontal:      [" << data.u_horizontal.size() << " x 1]" << std::endl
       << "    reynolds_stress:   [" << data.reynolds_stress.rows() << " x " << data.reynolds_stress.cols() << "]" << std::endl
       << "    reynolds_stress_a: [" << data.reynolds_stress_a.rows() << " x " << data.reynolds_stress_a.cols() << "]" << std::endl
       << "    reynolds_stress_b: [" << data.reynolds_stress_b.rows() << " x " << data.reynolds_stress_b.cols() << "]" << std::endl
       << "    k1z:               [" << data.k1z.rows() << " x " << data.k1z.cols() << "]" << std::endl
       << "    psi:               [" << tensor_dims(data.psi) << "]" << std::endl
       << "    psi_a:             [" << tensor_dims(data.psi_a) << "]" << std::endl
       << "    psi_b:             [" << tensor_dims(data.psi_b) << "]" << std::endl
       << "    t2wa:              [" << data.t2wa.size() << " x 1]" << std::endl
       << "    t2wb:              [" << data.t2wb.size() << " x 1]" << std::endl
       << "    residual_a:        [" << data.residual_a.size() << " x 1]" << std::endl
       << "    residual_b:        [" << data.residual_b.size() << " x 1]" << std::endl;
    return os;
}


/** @brief Get the T^2w distributions from the eddy signatures by deconvolution.
 *
 * These distributions are used for calculation of Spectra and Stress terms.
 *
 * Updates `t2wa`, `t2wb`, `lambda_e`, `residual_a` and `residual_b` in the input data structure.
 *
 * Legacy MATLAB equivalent is:
 * @code
 *    [t2wa, t2wb, lambda_e, residual_a, residual_b] = get_tw(pi_coles, s, beta, zeta, JA(:,3), JB(:,3), lambda);
 * @endcode
 *
 * @param data
 * @param[in] signature_a
 * @param[in] signature_b
 */
void get_t2w(AdemData& data, const EddySignature& signature_a, const EddySignature& signature_b) {

    // Define a range for lambda_e, choose from the smallest grid unit size to full b.l. scale (eta=1)
    Eigen::ArrayXd lambda_e = Eigen::ArrayXd::LinSpaced(500, 0, signature_a.lambda.maxCoeff());

    // Re-express as eta and flip so that eta ascends (required for the integration in reynolds_stress_13). Add a zero
    // first element, so that the integration goes from zero.
    Eigen::ArrayXd eta;
    eta = -1.0 * lambda_e; // *-1 inverts 1/eta in the subsequent exp() operator
    eta = eta.exp();

    Eigen::ArrayXd eta_with_zero = Eigen::ArrayXd(eta.rows()+1);
    eta_with_zero.setZero();
    eta_with_zero.bottomRows(eta.rows()) = eta.reverse();

    Figure figc = Figure();
    ScatterPlot pc = ScatterPlot();
    pc.x = Eigen::ArrayXd::LinSpaced(eta_with_zero.rows(), 1, eta_with_zero.rows());
    pc.y = eta_with_zero;
    figc.add(pc);
    figc.write("check_that_eta_ascends.json");

    Figure figc1 = Figure();
    ScatterPlot pc1 = ScatterPlot();
    pc1.x = Eigen::ArrayXd::LinSpaced(lambda_e.rows(), 1, lambda_e.rows());
    pc1.y = lambda_e;
    figc1.add(pc1);
    figc1.write("check_that_lambda_e_ascends.json");

    // Get the Reynolds Stresses, trim the zero point, and flip back
    Eigen::ArrayXd r13a;
    Eigen::ArrayXd r13b;
    reynolds_stress_13(r13a, r13b, data.beta, eta_with_zero, data.kappa, data.pi_coles, data.shear_ratio, data.zeta);
    r13a = r13a.bottomRows(eta.rows());
    r13b = r13b.bottomRows(eta.rows());
    r13a.reverseInPlace();
    r13b.reverseInPlace();

    Figure figr = Figure();
    ScatterPlot pr = ScatterPlot();
    pr.x =  Eigen::ArrayXd::LinSpaced(r13a.rows(), 1, r13a.rows());
    pr.y = r13a;
    figr.add(pr);
    figr.write("check_that_r13a_ascends_to_1_at_the_wall.json");
    Figure figrb = Figure();
    ScatterPlot prb = ScatterPlot();
    prb.x =  Eigen::ArrayXd::LinSpaced(r13b.rows(), 1, r13b.rows());
    prb.y = r13b;
    figr.add(prb);
    figr.write("check_that_r13b_tends_to_0_at_the_wall.json");

    // Show Reynolds Stresses behaving as part of validation
    cpplot::Figure fig = cpplot::Figure();
    cpplot::Layout lay = cpplot::Layout();
    cpplot::ScatterPlot pa = cpplot::ScatterPlot();
    pa.x = eta;
    pa.y = -1.0*r13a;
    pa.name = "Type A";
    cpplot::ScatterPlot pb = cpplot::ScatterPlot();
    pb.x = eta;
    pb.y = -1.0*r13b;
    pb.name = "Type B";
    cpplot::ScatterPlot pab = cpplot::ScatterPlot();
    pab.x = eta;
    pab.y = -1.0*(r13a + r13b);
    pab.name = "Total";
    fig.add(pa);
    fig.add(pb);
    fig.add(pab);
    lay.xTitle("$z/\\delta_{c}$");
    lay.yTitle("$-\\overline{u_1u_3}/U_\\tau^2$");
    fig.setLayout(lay);
    fig.write("check_r13_analytic.json");

    // Double check plot of J13
    Figure figc2 = Figure();
    ScatterPlot pc2 = ScatterPlot();
    pc2.x = signature_a.lambda;
    pc2.y = signature_a.j.col(2);
    pc2.name = "Type A";
    figc2.add(pc2);
    figc2.write("check_that_j13a_behaves.json");
    Figure figc3 = Figure();
    ScatterPlot pc3 = ScatterPlot();
    pc3.x = signature_a.lambda;
    pc3.y = signature_b.j.col(2);
    pc3.name = "Type B";
    figc3.add(pc3);
    Layout layc3 = Layout();
    layc3.xTitle("$\\lambda");
    layc3.xTitle("$J_{13}$");
    figc3.write("check_that_j13b_behaves.json");

    // Extract J13 signature terms and trim so that array sizes match after the deconv
//    auto len = signature_a.j.rows() - 2;
//    Eigen::ArrayXd j13a = signature_a.j.col(2).tail(len);
//    Eigen::ArrayXd j13b = signature_b.j.col(2).tail(len);
    auto len = signature_a.j.rows();
    Eigen::ArrayXd j13a = signature_a.j.col(2);
    Eigen::ArrayXd j13b = signature_b.j.col(2);
    Eigen::ArrayXd j13a_tmp(len);
    Eigen::ArrayXd j13b_tmp(len);



//    for (auto i = 0; i<len; i++) {
//        if (j13a(i) > -0.00001) {
//            j13a(i) = -0.00001;
//        }
//        if (j13b(i) > -0.00001) {
//            j13b(i) = -0.00001;
//        }
//    }

    // We multiply both Reynolds Stresses and Signatures by -1. The distributions are entirely negative, which
    // destabilises the deconvolution algorithm (it iteratively normalises on the first coefficient).
//    r13a = -1.0 * r13a;
//    r13b = -1.0 * r13b;
//    j13a = -1.0 * j13a * pow((0.25 / M_PI), 2.0);
//    j13b = -1.0 * j13b * pow((0.25 / M_PI), 2.0);

//    // Deconvolution is unstable where there are small values of j (a divide-by-0 error, effectively) so limit the
//    // kernel to sensible magnitude elements.
//    double cutoff_factor = 0.001;
//    double cutoff_j13a = j13a.minCoeff() * cutoff_factor;
//    double cutoff_j13b = j13b.minCoeff() * cutoff_factor;
//    Eigen::Index ctr = 0;
//    for (auto i = 0; i < j13a.rows(); i++) {
//        if ((j13a(i) <= cutoff_j13a) && (j13b(i) <= cutoff_j13b)) {
//            j13a_tmp(ctr) = j13a(i);
//            j13b_tmp(ctr) = j13b(i);
//            ctr++;
//        }
//    };
//    j13a = j13a_tmp.head(ctr);
//    j13b = j13b_tmp.head(ctr);

    // Deconvolve out the A and B structure contributions to the Reynolds Stresses
    // NOTE: it's actually -1*T^2w in this variable
    Eigen::ArrayXd minus_t2wa;
    Eigen::ArrayXd minus_t2wb;



    minus_t2wa = utilities::lowpass_fft_deconv(r13a, j13a, "Type_A");
    minus_t2wb = utilities::lowpass_fft_deconv(r13b, j13b, "Type_B");
//    deconv(minus_t2wa, r13a, j13a);
//    deconv(minus_t2wb, r13b, j13b);
//    std::cout << "r13a = [" << r13a << "];" << std::endl;
//    std::cout << "r13b = [" << r13b << "];" << std::endl;
//    std::cout << "j13a = [" << j13a << "];" << std::endl;
//    std::cout << "j13b = [" << j13b << "];" << std::endl;
//    std::cout << "minus_t2wa = [" << minus_t2wa << "];" << std::endl;
//    std::cout << "minus_t2wb = [" << minus_t2wb << "];" << std::endl;

    // Create a plot to show the t2w terms
    Figure figt = Figure();
    ScatterPlot pt = ScatterPlot();
    pt.x = lambda_e;
    pt.y = minus_t2wa;
    pt.name = "t2wa";
    figt.add(pt);
    ScatterPlot pt2 = ScatterPlot();
    pt2.x = lambda_e;
    pt2.y = minus_t2wb;
    pt2.name = "t2wb";
    figt.add(pt2);
    Layout layt = Layout();
    layt.xTitle("lambda_e");
    figt.setLayout(layt);
    figt.write("check_t2w_plot.json");



    // Extend by padding out to the same length as lambda_e.
    // The T^2w distributions converge to a constant at high lambda (close to the wall) so padding with the last value
    // in the vector is physically valid. This will result in the Reynolds Stresses and Spectra, which are obtained by
    // convolution, having the same number of points in the z direction as the axes variables (eta, lambda_e, etc)
//    double pad_value_a = minus_t2wa(minus_t2wa.rows()-1);
//    double pad_value_b = minus_t2wb(minus_t2wb.rows()-1);
//
//    auto last_n = lambda_e.rows() - minus_t2wa.rows();
//    minus_t2wa.conservativeResize(lambda_e.rows(),1);
//    minus_t2wb.conservativeResize(lambda_e.rows(),1);
//    minus_t2wa.tail(last_n) = pad_value_a;
//    minus_t2wb.tail(last_n) = pad_value_b;

    // Store in the data object
    data.eta = eta;
    data.lambda_e = lambda_e;
    data.t2wa = minus_t2wa.matrix();
    data.t2wb = minus_t2wb.matrix();

}


/** @brief Get the mean speed profile and update the data structure with it.
 *
 * @param data
 */
void get_mean_speed(AdemData& data) {
    data.u_horizontal = lewkowicz_speed(data.eta, data.pi_coles, data.kappa, data.u_inf, data.u_tau);
}


/** @brief Get the Reynolds Stress distributions from T2w and J distributions.
 *
 * The ouptut Reynolds Stress matrix is of size
 *   output_dim_size = input_dim_size - kernel_dim_size + 1 (requires: input_dim_size >= kernel_dim_size).
 * Legacy MATLAB equivalent is:
 * @code
 *      [R, RA, RB] = getReynoldsStresses(T2wA, T2wB, JA, JB);
 * @endcode
 *
 * @param data
 */
void get_reynolds_stresses(AdemData& data, const EddySignature& signature_a, const EddySignature& signature_b){

    // Map the input signature data to tensors (shared memory)
    auto input_rows = data.t2wa.rows();
    auto kernel_rows = signature_a.j.rows() - 2; // removes the first two rows of j for stability

    // Compute cumulative length of input and kernel;
    auto N = input_rows + kernel_rows - 1 ;
    auto M = fft_next_good_size(N);

    // Allocate the output and set zero
    data.reynolds_stress_a = Eigen::ArrayXXd(M, 6);
    data.reynolds_stress_b = Eigen::ArrayXXd(M, 6);
    data.reynolds_stress_a.setZero();
    data.reynolds_stress_b.setZero();

    Eigen::FFT<double> fft;

    // Column-wise convolution of T2w with J
    // TODO refactor to use the conv() function
    for (int k = 0; k < 6; k++) {
        Eigen::VectorXd in_a1(M);
        Eigen::VectorXd in_a2(M);
        Eigen::VectorXd in_b1(M);
        Eigen::VectorXd in_b2(M);
        in_a1.setZero();
        in_a2.setZero();
        in_b1.setZero();
        in_b2.setZero();
        in_a1.topRows(input_rows) = data.t2wa.matrix();
        in_b1.topRows(input_rows) = data.t2wb.matrix();
        in_a2.topRows(kernel_rows) = signature_a.j.col(k).bottomRows(kernel_rows).matrix();
        in_b2.topRows(kernel_rows) = signature_b.j.col(k).bottomRows(kernel_rows).matrix();

        // Take the forward ffts
        Eigen::VectorXcd out_a1(M);
        Eigen::VectorXcd out_a2(M);
        Eigen::VectorXcd out_b1(M);
        Eigen::VectorXcd out_b2(M);
        fft.fwd(out_a1, in_a1);
        fft.fwd(out_a2, in_a2);
        fft.fwd(out_b1, in_b1);
        fft.fwd(out_b2, in_b2);

        // Convolve by element-wise multiplication
        Eigen::VectorXcd inter_a = out_a1.array() * out_a2.array();
        Eigen::VectorXcd inter_b = out_b1.array() * out_b2.array();

        // Inverse FFT to deconvolve
        Eigen::VectorXd out_a(M);
        Eigen::VectorXd out_b(M);
        fft.inv(out_a, inter_a);
        fft.inv(out_b, inter_b);
        data.reynolds_stress_a.col(k) = out_a;
        data.reynolds_stress_b.col(k) = out_b;
    }

    // Trim the zero-padded ends
    data.reynolds_stress_a = data.reynolds_stress_a.topRows(data.lambda_e.rows());
    data.reynolds_stress_b = data.reynolds_stress_b.topRows(data.lambda_e.rows());
//    std::cout << "BODGED HERE _ REMOVE!!!!!" << std::endl;
    data.reynolds_stress = data.reynolds_stress_a + data.reynolds_stress_b;

    // FlipUD to match the reversal of z compared to the lambda_e basis on which these are calculated
    data.reynolds_stress_a = data.reynolds_stress_a.colwise().reverse().eval();
    data.reynolds_stress_b = data.reynolds_stress_b.colwise().reverse().eval();
    data.reynolds_stress = data.reynolds_stress.colwise().reverse().eval();

}


/** @brief Get the premultiplied power spectra.
 *
 * Computes tensors  @f$ \Psi_{a} @f$ (`psi_a`) and @f$ \Psi_{b} @f$ (`psi_b`) from the scale functions and
 * the eddy signatures, and adds them to the AdemData object.
 *
 * The output spectral tensors have dimension [nz x nk x 6], where `nz` should agree with the number of elements in
 * the `t2w` scale functions, and `nk` represents the number of wavenumbers for which the spectrum is computed.
 *
 * @param[inout] data AdemData object, which must have properties `t2wa`, `t2wb' and `u_tau` defined already.
 * @param[in] signature_a Signature object for Type A eddies
 * @param[in] signature_b Signature object for Type B eddies
 */
void get_spectra(AdemData& data, const EddySignature& signature_a, const EddySignature& signature_b){

    // Initialise the output spectrum tensors
    Eigen::array<Eigen::Index, 3> dims = signature_a.g.dimensions();
    Eigen::array<Eigen::Index, 3> psi_dims = {data.t2wa.rows(), dims[1], dims[2]};
    Eigen::Tensor<double, 3> psi_a(psi_dims);
    Eigen::Tensor<double, 3> psi_b(psi_dims);

    // Map the t2w arrays as vectors, for use with the conv function
    Eigen::Map<Eigen::VectorXd> t2wa_vec(data.t2wa.data(), data.t2wa.rows());
    Eigen::Map<Eigen::VectorXd> t2wb_vec(data.t2wb.data(), data.t2wb.rows());

    std::cout << "DIMS " << dims[0] << " " << dims[1] << " " << dims[2] << std::endl;
    std::cout << "PSI_DIMS " << psi_dims[0] << " " << psi_dims[1] << " " << psi_dims[2] << std::endl;

    // For each of the 6 auto / cross spectra terms
    for (Eigen::Index j = 0; j < dims[2]; j++) {

        auto page_offset_sig = j * dims[0] * dims[1];
        auto page_offset_psi = j * psi_dims[0] * psi_dims[1];

        // For each of the wavenumbers
        for (Eigen::Index i = 0; i < dims[1]; i++) {

            auto elements_offset_sig = (page_offset_sig + i * dims[0]);
            auto elements_offset_psi = (page_offset_psi + i * psi_dims[0]);

            Eigen::Map<Eigen::VectorXd> g_a_vec((double *)signature_a.g.data() + elements_offset_sig, dims[0]);
            Eigen::Map<Eigen::VectorXd> g_b_vec((double *)signature_b.g.data() + elements_offset_sig, dims[0]);

            Eigen::Map<Eigen::VectorXd> psi_a_vec((double *)psi_a.data() + elements_offset_psi, psi_dims[0]);
            Eigen::Map<Eigen::VectorXd> psi_b_vec((double *)psi_b.data() + elements_offset_psi, psi_dims[0]);

            psi_a_vec = conv(t2wa_vec, g_a_vec).reverse();
            psi_b_vec = conv(t2wb_vec, g_b_vec).reverse();

        }
    }

    // Don't store the sum of the spectra - they're memory hungry and we can always add the two together

    // Premultiply by the u_tau^2 term (see eq. 43)
    // auto u_tau_sqd = pow(data.u_tau, 2.0);
    // data.psi_a = u_tau_sqd * psi_a;
    // data.psi_b = u_tau_sqd * psi_b;
    data.psi_a = psi_a;
    data.psi_b = psi_b;

}


/** @brief Compute full turbulent properties from given Attached-Detached Eddy Model parameters.
 *
 * Uses the using the Lewkowicz (1982) formulation (Perry and Marusic eq.9) of the Coles wake function
 * to determine u_h(z), spectra and Reynolds Stresses from input parameters.
 *
 * @param[in] beta          The Clauser parameter, representing acceleration/decelaration of the boundary layer
 * @param[in] delta_c       The boundary layer thickness in m
 * @param[in] kappa         The von Karman constant, typically 0.41.
 * @param[in] pi_coles      Coles wake parameter Pi
 * @param[in] shear_ratio   Ratio between free stream and friction velocity S = u_inf/u_tau
 * @param[in] u_inf         The free stream speed in m/s
 * @param[in] zeta          Represents a scaled streamwise derivative of Pi
 * @param[in] signature_a   EddySignature class loaded with data for type A eddies
 * @param[in] signature_b   EddySignature class loaded with data for an ensemble of type B eddies
 */
AdemData adem(const double beta,
              const double delta_c,
              const double kappa,
              const double pi_coles,
              const double shear_ratio,
              const double u_inf,
              const double zeta,
              const EddySignature& signature_a,
              const EddySignature& signature_b,
              bool compute_spectra=true){

    // Data will contain computed outputs and useful small variables from the large signature files
    AdemData data = AdemData();

    // Add parameter inputs to data structure
    data.beta = beta;
    data.delta_c = delta_c;
    data.kappa = kappa;
    data.pi_coles = pi_coles;
    data.shear_ratio = shear_ratio;
    data.u_inf = u_inf;
    data.u_tau = data.u_inf / data.shear_ratio;
    data.zeta = zeta;

    // Deconvolve for T^2w and update the data structure with the outputs
    get_t2w(data, signature_a, signature_b);

    // Get vertical points through the boundary layer (for which the convolution functions were defined)
    data.z = data.eta.array() * data.delta_c;

    // Add k1z to the data structure
    Eigen::ArrayXd eta = data.eta.array();
    data.k1z = signature_a.k1z(eta);

    // Get the mean speed profile at the same vertical points and update the data structure
    get_mean_speed(data);

    // Determine Reynolds Stresses by convolution and add them to the data structure
    get_reynolds_stresses(data, signature_a, signature_b);

    // Determine Spectra by convolution and add them to the data structure
    if (compute_spectra) {
        get_spectra(data, signature_a, signature_b);
    }

    return data;

}

} /* namespace es */

#endif /* SOURCE_ADEM_H_ */
