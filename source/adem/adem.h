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
#include "io/variable_readers.h"
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

    /// Reynolds Stress profile @f$R_{13A}@f$ determined analytically from the parameter set
    Eigen::ArrayXXd r13a_analytic;

    /// Reynolds Stress profile @f$R_{13B}@f$ determined analytically from the parameter set
    Eigen::ArrayXXd r13b_analytic;

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

    Eigen::Index start_idx;




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

    /* Let's talk eddy scales.
     * lambda_e is the scale of a given eddy, which contributes to the stresses and spectra in the boundary layer.
     * lambda, as far as the signatures are concerned, is mapped for a domain spanning the actual eddy.
     * For a signature, values lambda < 0 are valid; since an eddy exerts an influence a short distance above itself.
     * lambda, as far as the deconvolution is concerned, is mapped over the whole boundary layer.
     * For a boundary layer profile, values < 0 are valid but worthless: R_ij and S_ij must definitively be 0 outside
     * the boundary layer.
     *
     * Here, we need a lambda distribution which is valid. So we have two strategies:
     *  - Take the part of signature.lambda which is >0 and use that for the deconvolution.
     *  - Take the whole of signature.lambda and ensure R_ij_ and S_ij are 0 outside the boundary layer.
     *
     * Let's do the first; as it'll lead to shorter t^2w vector and therefore be more computationally efficient at the
     * next step (which involves many convolutions).
     *
     * Why don't we define a totally separate lambda distribution? Because for the deconvolution we need to ensure the
     * distributions are mapped consistently (i.e. d_lambda is the same) to avoid biasing the result.
     */
    data.start_idx = 0;
    for (Eigen::Index i=0; i<signature_a.lambda.rows(); i++) {
        if (signature_a.lambda(i) >= 0) {
            data.start_idx = i;
            break;
        }
    }

    // TODO this is specific to the r_theta = 5395 P&M95b case
    double k_t = 1265;

//    double lambda_1 = 100 * nu / u_tau

//    Eigen::Index n_lambda = data.end_idx - data.start_idx + 1;
//    Eigen::ArrayXd lambda = signature_a.lambda.segment(data.start_idx, n_lambda);
//    std::cout << "here 2" << std::endl;
//    Eigen::ArrayXd eta = signature_a.eta.segment(data.start_idx, n_lambda);
//    std::cout << "here 3" << std::endl;

    /* So the smallest scale we want to go to is lambda_1, which is defined by the largest value in the lambda_e vector
     * (of the signature). When we deconvolve, we'd ideally do so at the same values of lambda as defined in the
     * signature lambda_e. But, that's unstable. Any residual from the deconvolution (e.g. where combinations of the
     * input signature kernel can't be used to exactly create the output solution) gets distributed across all scales,
     * creating a solution which is extremely unstable.
     *
     * To fix this, we upsample the solution function by some factor (10, in this case) compared to the signature.
     * The residual energy is then focused at the 10th harmonic, allowing it to be effectively removed with a low pass
     * filter.
     *
     * There still remains an issue of the smallest scales, where an instability remains (we can't resolve that scale
     * according to the Nyquist criterion). So, to achieve that, we also have to extend our solution function to half
     * the scale we are trying to resolve (i.e. double the value of lambda).
     */
    Eigen::Index n_lambda_signature = signature_a.lambda.rows() - data.start_idx;
    Eigen::ArrayXd lambda_signature = signature_a.lambda.bottomRows(n_lambda_signature);
    Eigen::ArrayXd eta_signature = signature_a.eta.bottomRows(n_lambda_signature);
    double d_lambda_signature = signature_a.domain_spacing(2);
    double d_lambda_fine = d_lambda_signature/10.0;
    double lambda_min = lambda_signature(0);
    double lambda_max = 2.0*lambda_signature(n_lambda_signature-1);
    Eigen::Index n_lambda_fine = round((lambda_max - lambda_min)/d_lambda_fine);

    Eigen::ArrayXd lambda_fine = Eigen::ArrayXd::LinSpaced(n_lambda_fine, lambda_min, lambda_max);
    Eigen::ArrayXd eta_fine = lambda_fine.exp().inverse();

    // Get an ascending version of eta, containing the point eta=0, for use in determining analytic R13 distributions
    Eigen::ArrayXd eta_with_zero = Eigen::ArrayXd(eta_fine.rows()+1);
    eta_with_zero.setZero();
    eta_with_zero.bottomRows(eta_fine.rows()) = eta_fine.reverse();

    // Produce visual check that eta is ascending exponentially
    Figure figc = Figure();
    ScatterPlot pc = ScatterPlot();
    pc.x = Eigen::ArrayXd::LinSpaced(eta_with_zero.rows(), 1, eta_with_zero.rows());
    pc.y = eta_with_zero;
    figc.add(pc);
    Layout layc = Layout("Check that eta ascends exponentially");
    layc.xTitle("Row number");
    layc.yTitle("$\\eta$");
    figc.setLayout(layc);
    figc.write("check_that_eta_ascends.json");

    // Produce a visual check that lambda is ascending linearly
    Figure figc1 = Figure();
    ScatterPlot pc1 = ScatterPlot();
    pc1.x = Eigen::ArrayXd::LinSpaced(lambda_fine.rows(), 1, lambda_fine.rows());
    pc1.y = lambda_fine;
    figc1.add(pc1);
    Layout layc1 = Layout("Check that lambda ascends linearly");
    layc1.xTitle("Row number");
    layc1.yTitle("$\\lambda$");
    figc1.setLayout(layc1);
    figc1.write("check_that_lambda_ascends.json");

    // Get Reynolds Stresses, trim the zero point, reverse back so the ordering is consistent with lambda coordinates
    Eigen::ArrayXd r13a_fine;
    Eigen::ArrayXd r13b_fine;
    reynolds_stress_13(r13a_fine, r13b_fine, data.beta, eta_with_zero, data.kappa, data.pi_coles, data.shear_ratio, data.zeta);
    r13a_fine = r13a_fine.bottomRows(eta_fine.rows());
    r13b_fine = r13b_fine.bottomRows(eta_fine.rows());
    r13a_fine.reverseInPlace();
    r13b_fine.reverseInPlace();

    /* On signatures:
     *
     * The signatures as produced by EddySignature().computeSignatures() extend out beyond lambda_e = 1. This is because
     * it's physically reasonable that eddies have some velocity influence above themselves. In effect, they accelerate
     * flow above them, and decelerate flow behind.
     *
     * In our implementation, signatures are computed up to z/delta = 1.5, beyond which the influence is small enough to
     * be considered zero. Having an eddy whose scale is of the same size as the boundary layer, would therefore create
     * an influence outside the boundary layer, breaking all assumptions about boundary layers ever.
     *
     * For us to maintain the fundamental boundary layer top condition that U = U_inf for all z/delta > 1, we must
     * alter our eddy model such that eddies can have no influence above their own scale. Otherwise, we can have no
     * eddies larger than z/delta=2/3 in our boundary layer. See the figures in the validation folder
     *      "... with_outer_influence.json"
     * for what happens in this case.
     *
     * So we make that assumption, which is equivalent to clipping the eddy signatures for values of lambda < 0.
     */
    Eigen::ArrayXXd ja = signature_a.j.bottomRows(n_lambda_signature);
    Eigen::ArrayXXd jb = signature_b.j.bottomRows(n_lambda_signature);
    Eigen::ArrayXd j13a = ja.col(2);
    Eigen::ArrayXd j13b = jb.col(2);

    // Double check J13 (although this should've already been checked at the eddy signature stage)
    Figure figj = Figure();
    ScatterPlot pja = ScatterPlot();
    pja.x = lambda_signature;
    pja.y = j13a;
    pja.name = "Type A";
    figj.add(pja);
    ScatterPlot pjb = ScatterPlot();
    pjb.x = lambda_signature;
    pjb.y = j13b;
    pjb.name = "Type B";
    figj.add(pjb);
    Layout layj = Layout();
    layj.xTitle("$\\lambda$");
    layj.yTitle("$J_{13}$");
    figj.setLayout(layj);
    figj.write("check_j13.json");

    // We multiply both Reynolds Stresses and Signatures by -1. The distributions are entirely negative, which
    // destabilises the deconvolution algorithm (it iteratively normalises on the first coefficient).
//    r13a = -1.0 * r13a;
//    r13b = -1.0 * r13b;
//    j13a = -1.0 * j13a * pow((0.25 / M_PI), 2.0);
//    j13b = -1.0 * j13b * pow((0.25 / M_PI), 2.0);

    // Deconvolve out the A and B structure contributions to the Reynolds Stresses
    // NOTE: it's actually -1*T^2w in this variable
    Eigen::ArrayXd minus_t2wa_fine;
    Eigen::ArrayXd minus_t2wb_fine;
    double stability = 0.005;
    minus_t2wa_fine = utilities::lowpass_fft_deconv(r13a_fine, j13a, "Type_A", stability);
    minus_t2wb_fine = utilities::lowpass_fft_deconv(r13b_fine, j13b, "Type_B", stability);

    // Now we want to map fine variables back onto the space where we're working
    Eigen::ArrayXd minus_t2wa(n_lambda_signature);
    Eigen::ArrayXd minus_t2wb(n_lambda_signature);
    Eigen::ArrayXd r13a(n_lambda_signature);
    Eigen::ArrayXd r13b(n_lambda_signature);
    for (auto ind=0; ind < n_lambda_signature; ind++) {
        minus_t2wa(ind) = minus_t2wa_fine(ind*10);
        minus_t2wb(ind) = minus_t2wb_fine(ind*10);
        r13a(ind) = r13a_fine(ind*10);
        r13b(ind) = r13b_fine(ind*10);
    };

    // Show Reynolds Stresses behaving as part of validation
    Figure fig = Figure();
    Layout lay = Layout();
    ScatterPlot pa = ScatterPlot();
    pa.x = eta_fine;
    pa.y = -1.0*r13a_fine;
    pa.name = "Type A";
    ScatterPlot pb = ScatterPlot();
    pb.x = eta_fine;
    pb.y = -1.0*r13b_fine;
    pb.name = "Type B";
    ScatterPlot pab = ScatterPlot();
    pab.x = eta_fine;
    pab.y = -1.0*(r13a_fine + r13b_fine);
    pab.name = "Total";
    fig.add(pa);
    fig.add(pb);
    fig.add(pab);
    lay.xTitle("$z/\\delta_{c}$");
    lay.yTitle("$-\\overline{u_1u_3}/U_\\tau^2$");
    fig.setLayout(lay);
    fig.write("check_r13_analytic.json");

    // Create a plot to show the t2w terms
    Figure figt = Figure();
    ScatterPlot pt = ScatterPlot();
    pt.x = lambda_signature;
    pt.y = minus_t2wa;
    pt.name = "$-T_{A}^2(\\lambda - \\lambda_E)\\omega_{A}(\\lambda-\\lambda_E)$";
    figt.add(pt);
    ScatterPlot pt2 = ScatterPlot();
    pt2.x = lambda_signature;
    pt2.y = minus_t2wb;
    pt2.name = "$-T_{B}^2(\\lambda - \\lambda_E)\\omega_{B}(\\lambda-\\lambda_E)$";
    figt.add(pt2);
    Layout layt = Layout();
    layt.xTitle("$\\lambda$");
    figt.setLayout(layt);
    figt.write("check_t2w.json");

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
    data.eta = eta_signature;
    data.r13a_analytic = r13a;
    data.r13b_analytic = r13b;
    // TODO it's not lambda e!!!
    data.lambda_e = lambda_signature;
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

    // Get eddy signatures and trim the initial terms. See get_t2w() comments for why we trim.
    Eigen::Index n_lambda = signature_a.j.rows() - data.start_idx;
    Eigen::ArrayXXd ja = signature_a.j.bottomRows(n_lambda);
    Eigen::ArrayXXd jb = signature_b.j.bottomRows(n_lambda);

    // Map the input signature data to tensors (shared memory)
    auto input_rows = data.t2wa.rows();
    auto kernel_rows = n_lambda;

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
        in_a2.topRows(kernel_rows) = ja.col(k).bottomRows(kernel_rows).matrix();
        in_b2.topRows(kernel_rows) = jb.col(k).bottomRows(kernel_rows).matrix();

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
    data.reynolds_stress_a = data.reynolds_stress_a.topRows(n_lambda);
    data.reynolds_stress_b = data.reynolds_stress_b.topRows(n_lambda);
    data.reynolds_stress = data.reynolds_stress_a + data.reynolds_stress_b;

    // Check that R13 stacks up against the analytical solution
    Figure fig = Figure();
    Layout lay = Layout();
    ScatterPlot paa = ScatterPlot();
    paa.x = data.eta;
    paa.y = -1.0*data.r13a_analytic;
    paa.name = "Type A - analytic";
    ScatterPlot pba = ScatterPlot();
    pba.x = data.eta;
    pba.y = -1.0*data.r13b_analytic;
    pba.name = "Type B - analytic";
    ScatterPlot par = ScatterPlot();
    par.x = data.eta;
    par.y = -1.0*(data.reynolds_stress_a.col(2));
    par.name = "Type A - reconstructed";
    ScatterPlot pbr = ScatterPlot();
    pbr.x = data.eta;
    pbr.y = -1.0*(data.reynolds_stress_b.col(2));
    pbr.name = "Type B - reconstructed";
    fig.add(paa);
    fig.add(pba);
    fig.add(par);
    fig.add(pbr);
    lay.xTitle("$z/\\delta_{c}$");
    lay.yTitle("$-\\overline{u_1u_3}/U_\\tau^2$");
    fig.setLayout(lay);
    fig.write("check_r13_analytic_and_reconstructed.json");
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
