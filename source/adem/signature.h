/*
 * signature.h Eddy Signatures for use with the Attached-Detached Eddy Method
 *
 * Author:              Tom Clark  (thclark @ github)
 *
 * Copyright (c) 2019 Octue Ltd. All Rights Reserved.
 *
 */

#ifndef SOURCE_ADEM_SIGNATURE_H_
#define SOURCE_ADEM_SIGNATURE_H_

#include <boost/algorithm/string/classification.hpp> // Include boost::for is_any_of
#include <boost/algorithm/string/split.hpp> // Include for boost::split
#include <Eigen/Dense>
#include <Eigen/Core>
#include <math.h>
#include <stdexcept>
#include <unsupported/Eigen/CXX11/Tensor>
#include <unsupported/Eigen/FFT>

#include "adem/biot_savart.h"
#include "profile.h"
#include "relations/stress.h"
#include "relations/velocity.h"
#include "utilities/filter.h"
#include "utilities/conv.h"
#include "utilities/interp.h"
#include "utilities/tensors.h"
#include "utilities/trapz.h"
#include "io/variable_readers.h"
#include "io/variable_writers.h"

#include "cpplot.h"


typedef Eigen::Array<double, 5, 3> Array53d;
typedef Eigen::Array<double, 3, 3> Array33d;


using namespace utilities;
namespace es {


/** @brief Data container for Eddy signature tensors.
 *
 */
class EddySignature {
public:

    /// Eddy type (or types, if ensembled) used to create the results
    std::string eddy_type;

    /// Mapped vertical coordinates used in the analysis
    Eigen::ArrayXd lambda;

    /// Unmapped (cartesian space) coordinates corresponding to the points in ``lambda``.
    Eigen::ArrayXd eta;

    /// g (6 coefficients at each vertical coord and wavenumber, e.g 50 x 801 x 6)
    Eigen::Tensor<double, 3> g;

    /** Eddy intensity functions $Jij(lambda)$, which is [n_lambda x 6] in size, with columns ordered as
     * `[J11 J12 J13 J22 J23 J33]`.
     * Note that the lower diagonal of the 3x3 Reynolds Stress Tensor is not included due to symmetry.
     */
    Eigen::Array<double, Eigen::Dynamic, 6> j;

    /** Spacing of the regular grid used to create the eddy intensity signatures [dx, dy, dlambda] (note for the third
     * coorinate, points are linearly spaced in a logarithmic domain, so this is d_lambda, not d_z).
     */
    Eigen::Array3d domain_spacing;

    /// Extents of the regular grid placed over the unit eddy to create the eddy intensity signatures, in [x_min, x_max; y_min, y_max; z_min, z_max] form
    Eigen::Array<double, 3, 2> domain_extents;

    /** @brief Load data from a *.mat file containing eddy signature data.
     *
     * TODO overload with load(std::vector<std::string> file_names, bool print_var = false){} to load and average
     * multiple signature files
     *
     * @param[in] file_name File name (including relative or absolute path)
     * @param[in] print_var Boolean, default true. Print variables as they are read in (not advised except for debugging!)
     */
    void load(std::string file_name, bool print_var = false) {
        std::cout << "Reading eddy signature data from file " << file_name << std::endl;

        // Open the MAT file for reading
        mat_t *matfp = Mat_Open(file_name.c_str(), MAT_ACC_RDONLY);
        if (matfp == NULL) {
            std::string msg = "Error reading MAT file: ";
            throw std::invalid_argument(msg + file_name);
        }

        // Use the variable readers to assist
        eddy_type = readString(matfp, "eddy_type", print_var);
        lambda = readArrayXd(matfp, "lambda", print_var);
        eta = readArrayXd(matfp, "eta", print_var);
        domain_spacing = readArray3d(matfp, "domain_spacing", print_var);
        domain_extents = readArray32d(matfp, "domain_extents", print_var);
        g = readTensor3d(matfp, "g", print_var);
        j = readArrayXXd(matfp, "j", print_var);

        // Close the file
        Mat_Close(matfp);
        std::cout << "Finished reading eddy signature (Type " + eddy_type + ")" << std::endl;
    }

    /** @brief Save eddy signature data to a *.mat file.
     *
     * @param[in] filename File name (including relative or absolute path)
     */
    void save(std::string file_name) {
        std::cout << "Writing signature data to file " << file_name << std::endl;

        mat_t *matfp = Mat_CreateVer(file_name.c_str(), NULL, MAT_FT_MAT73);
        if ( NULL == matfp ) {
            throw std::runtime_error("Unable to create/overwrite MAT file '" + file_name + "'");
        }

        // Use the variable writers to assist
        writeString(matfp, "eddy_type", eddy_type);
        writeArrayXd(matfp, "lambda", lambda);
        writeArrayXd(matfp, "eta", eta);
        writeArray3d(matfp, "domain_spacing", domain_spacing);
        writeArray32d(matfp, "domain_extents", domain_extents);
        writeTensor3d(matfp, "g", g);
        writeArrayXXd(matfp, "j", j);

        // Close the file
        Mat_Close(matfp);
        std::cout << "Finished writing eddy signature (Type " + eddy_type + ")" << std::endl;

    }

    /** @brief Define overloaded + (plus) operator for eddy signatures.
     *
     * @param[in] c The EddySignature to add.
     * @return A new EddySignature() with combined signatures of the two eddy types.
     */
    EddySignature operator+(const EddySignature& c) const
    {
        EddySignature result;
        result.eddy_type = this->eddy_type + "+" + c.eddy_type;
        result.eta = this->eta;
        result.lambda = this->lambda;
        result.domain_spacing = this->domain_spacing;
        result.domain_extents = this->domain_extents;
        result.g = (this->g + c.g);
        result.j = (this->j + c.j);
        return result;
    }

    /** @brief Define overloaded * (multiply) operator for eddy signatures.
     *
     * @param[in] a The number to multiply by
     * @return new EddySignature() whose signature (g, j) is element-wise multiplied by input a.
     */
    EddySignature operator*(double a) const
    {
        EddySignature result;
        result.eddy_type = "(" + this->eddy_type + ")*" + std::to_string(a);
        result.eta = this->eta;
        result.lambda = this->lambda;
        result.domain_spacing = this->domain_spacing;
        result.domain_extents = this->domain_extents;
        result.g = this->g;
        result.j = this->j;
        result.g = result.g * a;
        result.j = result.j * a;
        return result;
    }

    /** @brief Define overloaded / (divide) operator for eddy signatures.
     *
     * @param[in] denom A number to divide by
     * @return A new EddySignature() whose signature (g, j) is element-wise divided by input denom.
     */
    EddySignature operator/(double denom) const
    {
        EddySignature result;
        result.eddy_type = "(" + this->eddy_type + ")/" + std::to_string(denom);
        result.eta = this->eta;
        result.lambda = this->lambda;
        result.domain_spacing = this->domain_spacing;
        result.domain_extents = this->domain_extents;
        result.g = this->g;
        result.j = this->j;
        result.g = result.g / denom;
        result.j = result.j / denom;
        return result;
    }

    /** @brief Get k1z wavenumber array (wavenumber space for each vertical coord, e.g. 50 x 801).
     *
     * @param[in] eta vertical heights at which to get the k1z value, normalised (i.e. z/delta)
     * @return k1z the wavenumber-length matrix
     */
    Eigen::ArrayXXd k1z(Eigen::ArrayXd &eta) const {

        double dx = domain_spacing[0];
        auto nx = Eigen::Index((domain_extents(0,1) - domain_extents(0,0)) / dx) + 1;
        Eigen::ArrayXd k1_delta = Eigen::ArrayXd::LinSpaced(nx, 0, nx-1) * 2.0 * M_PI / dx;
        Eigen::ArrayXXd k1z = k1_delta.replicate(1, eta.rows()) * eta.transpose().replicate(k1_delta.rows(), 1);

        std::cout << "k1_delta_start " << k1_delta(0) << std::endl;
        std::cout << "k1_delta_end " << k1_delta(nx-1) << std::endl;
        std::cout << "eta_start " << eta(0) << std::endl;
        std::cout << "eta_end " << eta(eta.rows()-1) << std::endl;
        std::cout << "k1z_n_rows " << k1z.rows() << std::endl;
        std::cout << "k1z_n_cols " << k1z.cols() << std::endl;

        return k1z;
    }

    /** @brief Interpolate the signature to new locations in lambda or eta
     *
     * To avoid warping the reconstruction, we need to do the convolution and the deconvolution with
     * both the input and the signature on the same, equally spaced, basis. What if we want to do it on a basis
     * different to the positions where the signature was calculated?
     *
     * This function allows you to get the signature array, j, at different vertical locations
     *
     * @param[in] locations [N x 1] The lambda (or optionally eta) values at which to retrieve the signature j
     * @param[in] linear boolean If true, locations are on a linear basis (i.e. they are values of eta). Otherwise
     * (default) they are on a logarithmic basis (i.e. they are values of lambda).
     * @return [N x 6] signature array Jij interpolated to input locations
     */
    Eigen::ArrayXXd getJ(const Eigen::ArrayXd & locations, const bool linear=false) const {
        Eigen::ArrayXXd j_fine(locations.rows(), 6);
        if (linear) {
            // The new locations are values of eta

            // Assert locations are in the valid range...
            if ((locations < 0.0).any()) {
                throw std::invalid_argument("Input locations (eta values) must be defined for eta >= 0.0");
            }

            // Interpolate on an eta basis
            for (int k = 0; k < 6; k++) {
                utilities::CubicSplineInterpolant s(this->eta.matrix(), this->j.col(k).matrix());
                j_fine.col(k) = s(locations);
            }

        } else {
            // The new locations are values of lambda

            // Interpolate on a lambda basis
            for (int k = 0; k < 6; k++) {
                utilities::CubicSplineInterpolant s(this->lambda.matrix(), this->j.col(k).matrix());
                j_fine.col(k) = s(locations);
            }
        }
        return j_fine;
     };

    /** @brief Calculate eddy intensity functions @f$J_{i,j}@f$ and @f$g_{i,j}@f$.
     *
     * Eddy intensity functions @f$J_{i,j}(\lambda)@f$ are computed for type ``A``, ``B1``, ``B2``, ``B3`` or ``B4``,
     * which are the eddies described in Ref 1. These functions are the signatures (in terms of turbulent fluctuations)
     * that an individual structure will contribute to Reynolds Stress in a boundary layer. They are used in the
     * convolution integral (36) of Perry and Marusic (1995a) to determine Reynolds Stress profiles.
     *
     * Eddy spectral functions @f$g_{i,j}(k1z,\lambda)@f$ are then computed from @f$J@f$. These functions are the
     * signatures (in terms of turbulent fluctuations) that an individual structure will contribute to turbulent spectra
     * in a boundary layer. They are used in the convolution integral (43) of Perry and Marusic (1995a) to determine
     * Spectral Tensor profiles.
     *
     * This method updates class properties ``j``, ``lambda``, ``eta``, ``g``, ``eddy_type``, ``domain_spacing``,
     * ``domain_extents``.
     *
     * It works by computing induced velocity [u, v, w] of a single eddy structure on a regular grid [x,y,lambda],
     * normalised to the eddy scale and surrounding the unit eddy, then integrating for @f$J@f$ and taking ffts for
     * @f$g@f$.
     *
     * Perry AE and Marusic I (1995a) A wall-wake model for turbulent boundary layers. Part 1. Extension of the
     * attached eddy hypothesis J Fluid Mech vol 298 pp 361-388
     *
     * @param[in] type, string one of 'A', 'B1', 'B2', 'B3', 'B4'
     * @param[in] n_lambda, int number of points logarithmically spaced in the z direction, between the outer part of
     * the domain and a location very close to the wall. Default 200.
     *
     */
    void computeSignature(const std::string &type, const int n_lambda=200, const double dx=0.002) {

        // Set the type string
        eddy_type = type;

        // Express eddy structure as line vortices.  See Figure 13 Perry and Marusic 1995.
        Array33d eddy;
        if (type == "A") {
            eddy << 0, -0.8, 0,
                    1,  0,   1,
                    0,  0.8, 0;

        } else if (type == "B1") {
            eddy << 0.2, -0.15, 1,
                    0,    0,    0.5,
                    0.09, 0.15, 1;

        } else if (type == "B2") {
            eddy << 0,     0.15, 0.5,
                    0.2,   0,    1,
                    0.11, -0.15, 0.5;

        } else if (type == "B3") {
            eddy << 0.09, -0.15, 1,
                    0,     0,    0.5,
                    0.2,   0.15, 1;

        } else if (type == "B4") {
            eddy << 0.11, 0.15, 0.5,
                    0.20, 0,    1,
                    0,   -0.15, 0.5;
        } else {
            // TODO sort out cpplot import ordering so that I can import exceptions from our own damn library!
//            es::InvalidEddyTypeException e;
//            throw(e);
        }

        // Determine start and end nodes of the eddy...
        Eigen::Array3Xd startNodes(3,4);
        Eigen::Array3Xd endNodes(3,4);
        startNodes.leftCols(2) = eddy.topRows(2).transpose();
        endNodes.leftCols(2) = eddy.bottomRows(2).transpose();

        // ... and its reflection in the wall
        eddy.col(2) *= -1;
        startNodes.rightCols(2) = eddy.bottomRows(2).transpose();
        endNodes.rightCols(2) = eddy.topRows(2).transpose();

        //    % PLOT EDDY STRUCTURE SHAPE
        //    raiseFigure(['Type ' type ' Eddy Structure'])
        //    clf
        //    subplot(1,3,1)
        //    plot3([startNodes(:,1)'; endNodes(:,1)'],[startNodes(:,2)'; endNodes(:,2)'],[startNodes(:,3)'; endNodes(:,3)'],'k-')
        //    axis equal
        //    xlabel('x/\delta')
        //    ylabel('y/\delta')
        //    zlabel('z/\delta')

        // Set up influence domain

        /* Let's just think a second about frequency. For an ABL, the boundary layer is (say) delta = 1km thick,
         * and the wind speed is (say) 20m/s, and we want to resolve spectra to a frequency of (say) 10Hz,
         * then what eddy scale do we need to go down to?
         *      20/10 = 2 m ...
         * nondimensionalising, that's an eddy scale z/delta of 0.002. So our largest value of lambda_e (the eddy scale,
         * smaller as lambda_e increases) should correspond to this z/delta (or even finer scale).
         *
         * We thus choose z/delta = 0.002 as our finest eddy scale - that's why input param dx is set at 0.002.
         *
         */
        //double lambda_max = log(1/dx);
        std::cout << "UNBODGE THIS" << std::endl;
        double lambda_max = log(500);
        double lambda_min = log(1/1.5);

        // The domain extents are set accordingly
        domain_extents = Eigen::Array<double, 3, 2>::Zero();
        domain_extents << -4, 4,
                          -2, 2,
                           (1.0 / exp(lambda_max)), (1.0 / exp(lambda_min));

        // Logarithmically space the lambda coordinates
        lambda = Eigen::ArrayXd::LinSpaced(n_lambda, lambda_min, lambda_max);
        double d_lambda = lambda(1)-lambda(0);

        // Store real space vertical coordinates too
        eta = lambda.exp().inverse();

        /* So now we need to set the spacing in the physical domain
         * Wavenumber k1 2*pi/L = 2*pi*f/U
         *
         * So for the above estimated parameters of f = 10Hz, U = 20m/s, L = 2. Normalised to x/delta, this is 0.002.
         * So we need the grid spacing to be of a similar order in order to capture the appropriate wavenumbers...
         */
        domain_spacing << dx, dx, d_lambda;

        // A constant expression used for scaling the Reynolds Stress contributions
        constexpr double u0_sqd = 1.0 / (4.0 * M_PI * M_PI);

        // Streamwise and transverse grids (eta is z)
        Eigen::Index nx = ceil((domain_extents(0,1) - domain_extents(0,0)) / domain_spacing(0));
        Eigen::Index ny = ceil((domain_extents(1,1) - domain_extents(1,0)) / domain_spacing(1));
        Eigen::ArrayXXd x_vec = Eigen::ArrayXd::LinSpaced(nx, domain_extents(0,0), domain_extents(0,1));
        Eigen::ArrayXXd y_vec = Eigen::ArrayXd::LinSpaced(ny, domain_extents(1,0), domain_extents(1,1));

        // Produce x, y matrices, with y being the quickest changing (down columns)
        // Eigen syntax to replicate and interleave like this is v. awkward, so we have to write a loop. Sigh.
        Eigen::Index n_locations = x_vec.rows() * y_vec.rows();
        Eigen::Array3Xd locations(3, n_locations);
        Eigen::Index ctr = 0;
        for (auto i=0; i<x_vec.rows(); i++) {
            for (auto j=0; j<y_vec.rows(); j++) {
                locations(0, ctr) = x_vec(i);
                locations(1, ctr) = y_vec(j);
                ctr++;
            };
        };

        // Set unit circulation and the core radius of 0.05 as recommended by Perry and Marusic.
        Eigen::VectorXd gamma = Eigen::VectorXd::Ones(4);
        Eigen::VectorXd effective_core_radius_squared = Eigen::VectorXd::Ones(4) * pow(0.05, 2.0);

        j = Eigen::ArrayXXd(n_lambda,6);
        j.setZero();

        // For each vertical coordinate, determine the contribution to the eddy signature and spectra
        for (Eigen::Index lam_ctr=0; lam_ctr<n_lambda; lam_ctr++) {
            locations.row(2) = Eigen::ArrayXXd::Constant(1, n_locations, eta(lam_ctr));
            std::cout << "Computing signature at eta=" << eta(lam_ctr) << " (" << lam_ctr+1 << " of " << n_lambda << ")" << std::endl;
            Eigen::Array3Xd induction = NaiveBiotSavart(
                startNodes.matrix(),
                endNodes.matrix(),
                locations.matrix(),
                gamma,
                effective_core_radius_squared
            );

            // We don't want to copy or reallocate the memory... map into the induction
            Eigen::Map<Eigen::ArrayXXd, 0, Eigen::InnerStride<3>> u (induction.data(), y_vec.rows(), x_vec.rows());
            Eigen::Map<Eigen::ArrayXXd, 0, Eigen::InnerStride<3>> v (induction.data()+1, y_vec.rows(), x_vec.rows());
            Eigen::Map<Eigen::ArrayXXd, 0, Eigen::InnerStride<3>> w (induction.data()+2, y_vec.rows(), x_vec.rows());

            // Compute and plot Reynolds Stresses produced by a single eddy
            Eigen::ArrayXXd uu = u * u / u0_sqd;
            Eigen::ArrayXXd uv = u * v / u0_sqd;
            Eigen::ArrayXXd uw = u * w / u0_sqd;
            Eigen::ArrayXXd vv = v * v / u0_sqd;
            Eigen::ArrayXXd vw = v * w / u0_sqd;
            Eigen::ArrayXXd ww = w * w / u0_sqd;

            // We have constant spacing, so use quick version of trapz to integrate first across X, then down Y
            j.block<1,1>(lam_ctr, 0) = trapz(trapz(uu,2) * domain_spacing(0),1) * domain_spacing(1);
            j.block<1,1>(lam_ctr, 1) = trapz(trapz(uv,2) * domain_spacing(0),1) * domain_spacing(1);
            j.block<1,1>(lam_ctr, 2) = trapz(trapz(uw,2) * domain_spacing(0),1) * domain_spacing(1);
            j.block<1,1>(lam_ctr, 3) = trapz(trapz(vv,2) * domain_spacing(0),1) * domain_spacing(1);
            j.block<1,1>(lam_ctr, 4) = trapz(trapz(vw,2) * domain_spacing(0),1) * domain_spacing(1);
            j.block<1,1>(lam_ctr, 5) = trapz(trapz(ww,2) * domain_spacing(0),1) * domain_spacing(1);
        };


    }
};

} /* namespace es */

#endif /* SOURCE_ADEM_SIGNATURE_H_ */
