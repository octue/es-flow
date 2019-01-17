/*
 * adem.h
 *
 *  Created on: 29 Nov 2016
 *      Author: thc29
 */

#ifndef SOURCE_ADEM_H_
#define SOURCE_ADEM_H_

#include <boost/algorithm/string/classification.hpp> // Include boost::for is_any_of
#include <boost/algorithm/string/split.hpp> // Include for boost::split
#include <Eigen/Dense>
#include <Eigen/Core>
#include <stdexcept>
#include <unsupported/Eigen/CXX11/Tensor>

#include "variable_readers.h"
#include "profile.h"
#include "relations/velocity.h"

namespace es {


/** Data container for ADEM input parameters and results
 *
 */
class AdemData {
public:

    // Eddy types used to create the results
    std::vector<std::string> eddy_types = {"A", "B"};

    // Parameters from which the results were created
    double beta;
    double delta_c;
    double kappa = 0.41;
    double pi_coles;
    double shear_ratio;
    double u_inf;
    double u_tau;
    double zeta;

    // Mapped vertical coordinates used in the analysis
    Eigen::VectorXd z;
    Eigen::VectorXd eta;
    Eigen::VectorXd lambda_e;

    // Horizontal mean velocity
    Eigen::VectorXd u_horizontal;

    // Reynolds stresses
    Eigen::ArrayXXd reynolds_stress;
    Eigen::ArrayXXd reynolds_stress_a;
    Eigen::ArrayXXd reynolds_stress_b;

    // Wavenumber
    Eigen::ArrayXXd k1z;

    // Spectra
    Eigen::Tensor<double, 3> psi;
    Eigen::Tensor<double, 3> psi_a;
    Eigen::Tensor<double, 3> psi_b;

    // Convolution functions
    Eigen::VectorXd t2wa;
    Eigen::VectorXd t2wb;

    // Fit residuals
    Eigen::VectorXd residual_a;
    Eigen::VectorXd residual_b;

    /** Load data from a *.mat file containing eddy signature data
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

    /** Save eddy signature data to a *.mat file
     *
     * @param[in] filename File name (including relative or absolute path)
     */
    void save(std::string filename) {
        std::cout << "Writing signature data..." << std::endl;
        throw std::invalid_argument("Error writing mat file - function not implemented");
    }

};


/** Data container for Eddy signature tensors
 *
 */
class EddySignature {
public:

    // Eddy types used to create the results
    std::string eddy_type;

    // Mapped vertical coordinates used in the analysis (e.g. 1 x 50)
    Eigen::VectorXd lambda;

    // Wavenumber (wavenumber space for each vertical coord, e.g. 50 x 801)
    Eigen::ArrayXXd k1z;

    // g (6 coefficients at each vertical coord and wavenumber, e.g 50 x 801 x 6)
    Eigen::Tensor<double, 3> g;

    // J (6 coefficients at each vertical coord, e.g 50 x 6)
    Eigen::ArrayXXd j;

    /** Load data from a *.mat file containing eddy signature data
     *
     * @param[in] file_name File name (including relative or absolute path)
     * @param[in] print_var Boolean, default true. Print variables as they are read in (not advised except for debugging!)
     */
    void load(std::string file_name, bool print_var = true) {
        std::cout << "Reading eddy signature data from file " << file_name << std::endl;

        // Open the MAT file for reading
        mat_t *matfp = Mat_Open(file_name.c_str(), MAT_ACC_RDONLY);
        if (matfp == NULL) {
            std::string msg = "Error reading MAT file: ";
            throw std::invalid_argument(msg + file_name);
        }

        // Use the variable readers to assist
        eddy_type = readString(matfp, "type", print_var);
        lambda = readVectorXd(matfp, "lambda", print_var);
        k1z = readArrayXXd(matfp, "k1z", print_var);
        g = readTensor3d(matfp, "g", print_var);
        j = readArrayXXd(matfp, "J", print_var);

        // Close the file
        Mat_Close(matfp);
        std::cout << "Finished reading eddy signature (Type " + eddy_type + ")" << std::endl;
    }

    /** Save eddy signature data to a *.mat file
     *
     * @param[in] filename File name (including relative or absolute path)
     */
    void save(std::string filename) {
        std::cout << "Writing signature data..." << std::endl;
        throw std::invalid_argument("Error writing mat file - function not implemented");
    }

    EddySignature operator+(const EddySignature&) const;       // operator+()

};

/** Define overloaded + (plus) operator for eddy signatures
 *
 */
EddySignature EddySignature::operator+ (const EddySignature& c) const
{
    EddySignature result;
    result.eddy_type = "Merged";
    // TODO assert equality of lambda and k1z
    result.lambda = this->lambda;
    result.k1z = this->k1z;
    result.g = (this->g + c.g);
    result.j = (this->j + c.j);
    return result;
}


/** Get the T^2omega distributions from the eddy signatures by deconvolution
 *
 * Updates t2wa, t2wb, lambda_e, residual_a and residual_b in the input data structure.
 *
 * Matlab equivalent is:
 *    [t2wa, t2wb, lambda_e, residual_a, residual_b] = get_tw(pi_coles, s, beta, zeta, JA(:,3), JB(:,3), lambda);
 *
 * @param data
 * @param[in] signature_a
 * @param[in] signature_b
 */
void get_t2w(AdemData& data, const EddySignature& signature_a, const EddySignature& signature_b) {}


/** Get the mean speed profile and update the data structure with it
 *
 * @param data
 */
void get_mean_speed(AdemData& data) {
    data.u_horizontal = lewkowicz_speed(data.eta, data.pi_coles, data.kappa, data.u_inf, data.u_tau);
}


/** Get the Reynolds Stress distributions from T2w and J distributions
 * Matlab equivalent is:
 *      [R, RA, RB] = getReynoldsStresses(T2wA, T2wB, JA, JB);
 *
 * @param data
 */
void get_reynolds_stresses(AdemData& data){}


/** Get the premultiplied power spectra
 * Matlab equivalent is:
 *      [Psi, PsiA, PsiB] = getSpectra(T2wA, T2wB, gA, gB, U1, S);
 * @param data
 */
void get_spectra(AdemData& data){}


/// Compute full turbulent properties from given Attached-Detached Eddy Model parameters
/**
 * Uses the using the Lewkowicz (1982) formulation (Perry and Marusic eq.9) of the Coles wake function
 * to determine u_h(z), spectra and Reynolds Stresses from input parameters.
 *
 * @param[in]  z            Height(s) in m at which output profiles are required
 * @param[in]  delta_c      The boundary layer thickness in m
 * @param[in]  u_inf        The free stream speed in m/s
 * @param[in]  pi_coles     Coles wake parameter Pi
 * @param[in]  shear_ratio  Ratio between free stream and friction velocity S = u_inf/u_tau
 * @param[in]  zeta         Represents a scaled streamwise derivative of Pi
 * @param[in]  beta         The Clauser parameter, representing acceleration/decelaration of the boundary layer
 * @param[in]  signature_a  EddySignature class loaded with data for type A eddies
 * @param[in]  signature_b  EddySignature class loaded with data for an ensemble of type B eddies
 */
AdemData adem(const double beta,
              const double delta_c,
              const double kappa,
              const double pi_coles,
              const double shear_ratio,
              const double u_inf,
              const double zeta,
              const EddySignature& signature_a,
              const EddySignature& signature_b){

    // Data will contain computed outputs and useful small variables from the large signature files
    AdemData data = AdemData();

    // Add k1z to the data structure
    data.k1z = signature_a.k1z;

    // Add parameter inputs to data structure
    data.beta = beta;
    data.delta_c = delta_c;
    data.kappa = kappa;
    data.pi_coles = pi_coles;
    data.shear_ratio = shear_ratio;
    data.u_inf = u_inf;
    data.zeta = zeta;

    // Deconvolve for T^2w and update the data structure with the outputs
    get_t2w(data, signature_a, signature_b);

    // Get vertical points through the boundary layer (for which the convolution functions were defined)
    Eigen::VectorXd eta = data.lambda_e;
    eta.array().exp();
    eta.array().inverse();
    data.eta = eta;
    data.z = eta.array() * data.delta_c;

    // Get the mean speed profile at the same vertical points and update the data structure
    get_mean_speed(data);

    // Determine Reynolds Stresses by convolution and add them to the data structure
    get_reynolds_stresses(data);

    // Determine Spectra by convolution and add them to the data structure
    get_spectra(data);

    return data;
}


} /* namespace es */

#endif /* SOURCE_ADEM_H_ */
