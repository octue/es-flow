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
#include <stdexcept>
#include <unsupported/Eigen/CXX11/Tensor>
#include <math.h>

#include "variable_readers.h"
#include "profile.h"
#include "relations/stress.h"
#include "relations/velocity.h"
#include "utilities/filter.h"

#include "cpplot.h"
#include <unsupported/Eigen/FFT>
#include "utilities/conv.h"
#include "utilities/tensors.h"

using namespace utilities;


namespace es {


/** @brief Data container for Eddy signature tensors.
 *
 */
class EddySignature {
public:

    /// Eddy types used to create the results
    std::string eddy_type;

    /// Mapped vertical coordinates used in the analysis (e.g. 1 x 50)
    Eigen::VectorXd lambda;

    /// g (6 coefficients at each vertical coord and wavenumber, e.g 50 x 801 x 6)
    Eigen::Tensor<double, 3> g;

    /// J (6 coefficients at each vertical coord, e.g 50 x 6)
    Eigen::ArrayXXd j;

    /// Spacing of the regular grid used to create the eddy intensity signatures, in [x, y, z] directions
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
        eddy_type = readString(matfp, "type", print_var);
        lambda = readVectorXd(matfp, "lambda", print_var);
        domain_spacing = readArray3d(matfp, "domain_spacing", print_var);
        domain_extents = readArray32d(matfp, "domain_extents", print_var);
        g = readTensor3d(matfp, "g", print_var);
        j = readArrayXXd(matfp, "J", print_var);

        // Close the file
        Mat_Close(matfp);
        std::cout << "Finished reading eddy signature (Type " + eddy_type + ")" << std::endl;
    }

    /** @brief Save eddy signature data to a *.mat file.
     *
     * @note NOT IMPLEMENTED YET
     *
     * @param[in] filename File name (including relative or absolute path)
     */
    void save(std::string filename) {
        std::cout << "Writing signature data..." << std::endl;
        throw std::invalid_argument("Error writing mat file - function not implemented");
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
        result.lambda = this->lambda;
        result.domain_spacing = this->domain_spacing;
        result.domain_extents = this->domain_extents;
        result.g = (this->g + c.g);
        result.j = (this->j + c.j);
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
        result.lambda = this->lambda;
        result.domain_spacing = this->domain_spacing;
        result.domain_extents = this->domain_extents;
        result.g = this->g;
        result.j = this->j;
        result.g = result.g / denom;
        result.j = result.j / denom;
        return result;
    }

    /** @brief Get k1z wavenumber array (wavenumber space for each vertical coord, e.g. 50 x 801)
     *
     * @param[in] eta vertical heights at which to get the k1z value, normalised (i.e. z/delta)
     * @return k1z the wavenumber-length matrix
     */
    Eigen::ArrayXXd k1z(Eigen::ArrayXd &eta) const {

        double dx = domain_spacing[0];
        auto nx = Eigen::Index((domain_extents(0,1) - domain_extents(0,0)) / dx) + 1;
        Eigen::ArrayXd k1_delta = Eigen::ArrayXd::LinSpaced(nx, 0, nx-1) * 2.0 * M_PI / dx;
        std::cout << "k1_delta_start " << k1_delta(0) << std::endl;
        std::cout << "k1_delta_end " << k1_delta(nx-1) << std::endl;
        std::cout << "eta_start " << eta(0) << std::endl;
        std::cout << "eta_end " << eta(eta.rows()-1) << std::endl;


        Eigen::ArrayXXd k1z = k1_delta.replicate(1, eta.rows()) * eta.transpose().replicate(k1_delta.rows(), 1);
        std::cout << "k1z_n_rows " << k1z.rows() << std::endl;
        std::cout << "k1z_n_cols " << k1z.cols() << std::endl;

//        Eigen::ArrayXXd k1z = k1_delta.transpose().replicate(eta.rows(), 1) * eta.replicate(1, k1_delta.rows());

        return k1z;
    }
};

} /* namespace es */

#endif /* SOURCE_ADEM_SIGNATURE_H_ */
