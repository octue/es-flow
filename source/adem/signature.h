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
     * coordinate, points are linearly spaced in a logarithmic domain, so this is d_lambda, not d_z).
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

    /** @brief Interpolate the signature @f$J_{i,j}(\lambda)@f$ to new locations in lambda or eta.
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

    /** @brief Apply signature defaults.
     *
     * Eddy intensity functions @f$J_{i,j}(\lambda)@f$ is applied for type ``A`` or ``B`` eddies, by digitisation of
     * Figure 20 in Perry and Marusic (1995).
     *
     * Note that these signatures are incomplete:
     *  1. Only four terms are given (I11, I13, I22, I33, the remaining two are assumed to be zero)
     *
     *  2. Eddy spectral functions @f$g_{i,j}(k1z,\lambda)@f$ are not documented in P&M 1995, so can't be reproduced
     *     here. Thus, signatures defined by this method cannot be used for computation of spectra (and sets arbitrary
     *     ``domain_spacing`` and ``domain_extents`` values for the ``x`` and ``y`` directions).
     *
     * This method updates class properties ``j``, ``lambda``, ``eta``, ``eddy_type``, ``domain_spacing``,
     * ``domain_extents``.
     *
     * See Perry AE and Marusic I (1995) A wall-wake model for turbulent boundary layers. Part 1. Extension of the
     * attached eddy hypothesis J Fluid Mech vol 298 pp 361-388
     *
     * @param[in] type, string one of 'A', 'B1', 'B2', 'B3', 'B4'
     * @param[in] n_lambda, int number of points logarithmically spaced in the z direction, between the outer part of
     * the domain and a location very close to the wall. Default 200.
     *
     */
    void applySignature(const std::string &type, const int n_lambda=200) {

        // Set the type string
        eddy_type = type;

        // Set arbitrary domain extents
        double lambda_max = log(500);
        double lambda_min = log(1/1.5);
        domain_extents = Eigen::Array<double, 3, 2>::Zero();
        domain_extents << -4, 4,
            -2, 2,
            (1.0 / exp(lambda_max)), (1.0 / exp(lambda_min));

        // Logarithmically space lambda coordinates
        lambda = Eigen::ArrayXd::LinSpaced(n_lambda, lambda_min, lambda_max);
        domain_spacing << 1.0, 1.0, lambda(1)-lambda(0);

        // Store real space vertical coordinates too
        eta = lambda.exp().inverse();

        // Interpolate the signatures that are given.
        j = Eigen::ArrayXXd(n_lambda,6);
        j.setZero();
        Eigen::ArrayXXd i11;
        Eigen::ArrayXXd i13;
        Eigen::ArrayXXd i22;
        Eigen::ArrayXXd i33;
        if (type == "A") {
            i11 = Eigen::ArrayXXd(2, 26);
            i11 << 0, 0.012710542, 0.030316638, 0.05182384, 0.08898147, 0.1300581, 0.18483716, 0.2357049, 0.29636374, 0.35898846, 0.41573855, 0.4783607, 0.5585967, 0.6368718, 0.7073141, 0.76600707, 0.82468724, 0.87943816, 0.91658044, 0.9419913, 0.97133654, 1.0007021, 1.0516082, 1.1260028, 1.2239062, 1.50,
                1.6593906, 1.6418779, 1.611259, 1.5544281, 1.4713508, 1.3926289, 1.3051336, 1.2263831, 1.1476042, 1.08192, 1.0162529, 0.94620186, 0.8586324, 0.7667018, 0.6747941, 0.5829205, 0.469213, 0.33368298, 0.22440498, 0.1457287, 0.09760853, 0.084422655, 0.07117407, 0.040389933, 0.027004533, 0.0;

            i13 = Eigen::ArrayXXd(2, 21);
            i13 << 0, 0.023503713, 0.05093406, 0.08618198, 0.12728672, 0.17620902, 0.22121488, 0.28186864, 0.3620689, 0.45790172, 0.55568004, 0.64367235, 0.7238267, 0.80006206, 0.86259496, 0.90362066, 0.93487054, 0.9700394, 0.9993566, 1.0423712, 1.50,
                0, -0.061066702, -0.14832275, -0.22682244, -0.2878379, -0.340097, -0.38363394, -0.43149212, -0.4618262, -0.4702808, -0.46126255, -0.43917242, -0.39090434, -0.32954726, -0.24639612, -0.17204118, -0.102081485, -0.045210768, -0.014557855, -0.0030345472, -0.004372481;

            i22 = Eigen::ArrayXXd(2, 23);
            i22 << 0, 0.018947192, 0.03648182, 0.054011352, 0.0754471, 0.10665094, 0.17303112, 0.21991083, 0.28439146, 0.34693173, 0.41142765, 0.48179564, 0.57952243, 0.65573704, 0.7202125, 0.788579, 0.8452064, 0.9115866, 0.9701798, 1.0483195, 1.1304111, 1.2340457, 1.50,
                1.176464, 1.150218, 1.0934712, 1.0323671, 0.966883, 0.8926276, 0.79638124, 0.74817854, 0.6998736, 0.6646518, 0.6294187, 0.5985088, 0.5500107, 0.5016375, 0.4489753, 0.37886125, 0.3044581, 0.20821178, 0.14251183, 0.067983694, 0.028291034, 0.014617034, 0.0043686354;

            i33 = Eigen::ArrayXXd(2, 32);
            i33 << 0.0, 0.023483196, 0.041119818, 0.06269325, 0.09603633, 0.12936921, 0.17444114, 0.21754721, 0.25868744, 0.30175778, 0.35461667, 0.4133425, 0.46618098, 0.5346596, 0.608995, 0.6696255, 0.73806334, 0.7908406, 0.83576465, 0.8884807, 0.92553115, 0.95087314, 0.9898843, 1.0288852, 1.0737991, 1.1304673, 1.1793332, 1.2321156, 1.2849082, 1.3552966, 1.4374342, 1.50,
                0.0, 0.012935479, 0.043334138, 0.095496446, 0.1780915, 0.251972, 0.3301416, 0.39960802, 0.46037126, 0.49933675, 0.54696, 0.5945491, 0.6247433, 0.6504893, 0.6674866, 0.6714917, 0.66237944, 0.6402863, 0.59209496, 0.51771456, 0.42599595, 0.35613185, 0.26875913, 0.17267188, 0.115766004, 0.07622104, 0.05415063, 0.036414765, 0.027393447, 0.013912599, 0.013435401, 0.0043572737;

        } else if (type == "B") {
            i11 = Eigen::ArrayXXd(2, 23);
            i11 << 0, 0.07068844, 0.1666695, 0.27835685, 0.36659724, 0.4334307, 0.48672056, 0.5342177, 0.5837344, 0.62532634, 0.67665803, 0.74737716, 0.8137505, 0.8545449, 0.8894136, 0.93391985, 0.9706649, 1.00945, 1.0521617, 1.1184429, 1.2122818, 1.3199301, 1.50,
                0.016322415, 0.016322415, 0.018703382, 0.024518738, 0.03472676, 0.056264196, 0.091715895, 0.13412246, 0.18173045, 0.22154763, 0.25700384, 0.27592924, 0.25842202, 0.23056176, 0.19837682, 0.15315433, 0.11402357, 0.08182956, 0.050494146, 0.025177978, 0.011945607, 0.0073581133, 0.0017353186;

            i13 = Eigen::ArrayXXd(2, 14);
            i13 << 0, 0.2741542, 0.34076267, 0.42109883, 0.47993597, 0.56440324, 0.6233631, 0.7057494, 0.7918987, 0.8779049, 0.94425774, 1.0145372, 1.1240618, 1.50,
                0, -2.3333919E-4, -0.002682268, -0.0068347994, -0.014507807, -0.036872122, -0.054957043, -0.06691398, -0.06584696, -0.05263271, -0.033390157, -0.015006201, -0.0034729026, -8.676593E-4;

            i22 = Eigen::ArrayXXd(2, 23);
            i22 << 0, 0.14884579, 0.29180932, 0.4054613, 0.46828458, 0.515565, 0.5411826, 0.5746176, 0.6335705, 0.6963938, 0.74149364, 0.7689007, 0.8080032, 0.8490102, 0.89000964, 0.93290585, 0.9660264, 0.99715817, 1.0224689, 1.049784, 1.0810461, 1.1632366, 1.50,
                0, 0.0023498295, 0.0038403252, 0.012338308, 0.030489888, 0.06258152, 0.08079781, 0.09726137, 0.12063707, 0.13878864, 0.14566672, 0.14474949, 0.13772498, 0.12461021, 0.110625885, 0.08968175, 0.07049375, 0.047830947, 0.031265218, 0.019913385, 0.012032943, 0.005803057, 0.0026086513;

            i33 = Eigen::ArrayXXd(2, 20);
            i33 << 0, 0.21936293, 0.31337926, 0.39961416, 0.44083738, 0.48011523, 0.5254081, 0.56480104, 0.64147526, 0.64147526, 0.70819265, 0.7767404, 0.8334498, 0.8939988, 0.9446548, 0.9952725, 1.0459669, 1.0967507, 1.1632637, 1.5000255,
                0, 0.002830162, 0.004290453, 0.009236455, 0.016043989, 0.023722952, 0.0409085, 0.05637945, 0.07693768, 0.07693768, 0.08626908, 0.08693706, 0.081578515, 0.07101401, 0.053551547, 0.033491757, 0.018626625, 0.009821928, 0.0053009023, 0.0017315528;

        } else {
            throw std::invalid_argument("You can only apply signatures for type 'A' or type 'B' eddies");
        }

        utilities::LinearInterpolant s11(i11.row(0).transpose(), i11.row(1).transpose());
        utilities::LinearInterpolant s13(i13.row(0).transpose(), i13.row(1).transpose());
        utilities::LinearInterpolant s22(i22.row(0).transpose(), i22.row(1).transpose());
        utilities::LinearInterpolant s33(i33.row(0).transpose(), i33.row(1).transpose());

        j.col(0) = s11(eta);
        j.col(2) = s13(eta);
        j.col(3) = s22(eta);
        j.col(5) = s33(eta);
    }
};

} /* namespace es */

#endif /* SOURCE_ADEM_SIGNATURE_H_ */
