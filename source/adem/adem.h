/*
 * adem.h
 *
 *  Created on: 29 Nov 2016
 *      Author: thc29
 */

#ifndef SOURCE_ADEM_H_
#define SOURCE_ADEM_H_

#include <Eigen/Dense>
#include <Eigen/Core>
#include <stdexcept>
#include <unsupported/Eigen/CXX11/Tensor>

#include "variable_readers.h"

using Eigen::Array;
using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::Dynamic;

#include "profile.h"

using namespace Eigen;

namespace es {


/// Data container for ADEM input parameters and results
/**
 *
 */
class AdemData {
public:

    // Eddy types used to create the results
    std::vector<std::string> eddy_types = {"A", "B1", "B2", "B3", "B4"};

    // Parameters from which the results were created
    Eigen::VectorXd z = Eigen::VectorXd::LinSpaced(3, 1.0, 3.0);
    double u_inf;
    double pi_coles;
    double s;
    double u_tau;
    double beta;
    double zeta;

    // Mapped vertical coordinates used in the analysis
    Eigen::VectorXd eta;
    Eigen::VectorXd lambda_e;

    // Horizontal mean velocity
    Eigen::VectorXd u_h;

    // Reynolds stresses
    Eigen::MatrixXd reynolds_stress;
    Eigen::MatrixXd reynolds_stress_a;
    Eigen::MatrixXd reynolds_stress_b;

    // Wavenumber
    Eigen::VectorXd k1z;

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

};

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
    Eigen::Matrix3Xd j;

    /** Load data from a *.mat file containing eddy signature data
     *
     * @param[in] filename File name (including relative or absolute path)
     * @param[in] print_var Boolean, default false. Print variables as they are read in (not advised except for debugging!)
     */
    void load(std::string filename, bool print_var = false) {
        std::cout << "Reading eddy signature data..." << std::endl;

        // Open the MAT file for reading
        mat_t *matfp = Mat_Open(filename.c_str(), MAT_ACC_RDONLY);
        if (matfp == NULL) {
            std::string msg = "Error reading MAT file: ";
            throw std::invalid_argument(msg + filename);
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

};


/// Compute full turbulent properties from given Attached-Detached Eddy Model parameters
/**
 * Uses the using the Lewkowicz (1982) formulation (Perry and Marusic eq.9) of the Coles wake function
 * to determine u_h(z), spectra and Reynolds Stresses from input parameters.
 *
 * @param[in]  z        Height(s) in m at which output profiles are required
 * @param[in]  delta_c  The boundary layer thickness in m
 * @param[in]  u_inf    The free stream speed in m/s
 * @param[in]  pi_coles Coles wake parameter Pi
 * @param[in]  s        Ratio between free stream and friction velocity S = u_inf/u_tau
 * @param[in]  zeta     Represents a scaled streamwise derivative of Pi
 * @param[in]  beta     The Clauser parameter, representing acceleration/decelaration of the boundary layer
 * @param[in]  eddyFile Name of the eddy signatures file. If not given, use default eddySignatures.mat in the repository. If not present, compute eddySignatures file.
 */
AdemData adem(const double delta_c,
              const double u_inf,
              const double pi_coles,
              const double s,
              const double beta,
              const double zeta,
              EddySignature signature){

//    // Load signature from file or memory
//    k1z = signature.k1z;
//    lambda = signature.lambda;
//
//    // Get the mean profile at the same vertical points
//    eta = 1./exp(lambdaE);
//    z = eta*deltac;
//    [Ux] = getMeanProfile(Pi, S, deltac, U1, z);
//
//    // Signatures of Type A and Type B eddies (the latter are ensemble averaged)
//    g_b = signature.g_b;
//    j_b = signature.j_b;
//
//    // Deconvolve for T^2w
//    [t2wa, t2wb, lambda_e, residual_a, residual_b] = get_tw(pi_coles, s, beta, zeta, JA(:,3), JB(:,3), lambda);


    // Store output results (both computed and useful small variables from the large signatures file)
    AdemData results = AdemData();
//    results.k1z = signature.k1z;
//    results.lambda_e = lambda_e


    return results;
}




} /* namespace es */

#endif /* SOURCE_ADEM_H_ */
