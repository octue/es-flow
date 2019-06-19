/*
 * interp.h Brief overview sentence
 * 
 * Author:              Tom Clark  (thclark @ github)
 *
 * Copyright (c) 2019 Octue Ltd. All Rights Reserved.
 *
 */

#ifndef ES_FLOW_INTERP_H
#define ES_FLOW_INTERP_H
#include <Eigen/Core>
#include <unsupported/Eigen/Splines>

#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <unsupported/Eigen/FFT>
#include <math.h>


namespace utilities {

class CubicSplineInterpolant {
public:
    /** @brief Performs cubic interpolation (or maximum degree possible, where inputs are shorter than 4 elements)
     *
     * Initialise an interpolant as follows:
     *
     *    Eigen::VectorXd xvals(3);
     *    Eigen::VectorXd yvals(xvals.rows());
     *
     *    xvals << 0, 15, 30;
     *    yvals << 0, 12, 17;
     *
     *    CubicSplineInterpolant s(xvals, yvals);
     *
     */
    CubicSplineInterpolant(Eigen::VectorXd const &x_vec, Eigen::VectorXd const &y_vec)
        : x_min(x_vec.minCoeff()),
          x_max(x_vec.maxCoeff()),
          spline_(Eigen::SplineFitting<Eigen::Spline<double, 1>>::Interpolate(
              y_vec.transpose(),
              // No more than cubic spline, but accept short vectors.
              std::min<long>(x_vec.rows() - 1, 3),
              scaled_values(x_vec))) {}

    /** @brief Evaluate interpolant at values xi whose yi values are unknown
     *
     * Performs cubic interpolation (or maximum degree possible, where inputs are shorter than 4 elements)
     *
     * Initialise and evaluate an interpolant:
     *
     *    // ... Initialise interpolant (see constructor) ...
     *    CubicSplineInterpolant s(xvals, yvals);
     *    std::cout << s(12.34) << std::endl;
     *
     *@param xi double (or eigen VectorXd) of target x values for the interpolation
     *@return double interpolated value yi corresponding to location xi
     */
    double operator()(double xi) const {
        // x values need to be scaled down in extraction as well.
        return spline_(scaled_value(xi))(0);
    }

    Eigen::VectorXd operator()(Eigen::VectorXd const &xi_vec) {
        Eigen::VectorXd yi_vec(xi_vec.rows());
        yi_vec = xi_vec.unaryExpr([this](double xi) { return spline_(scaled_value(xi))(0); }).transpose();
        return yi_vec;
    }

private:
    /** @brief Scales x values to the domain [0, 1] on which the interpolant was defined
     *
     *@param x double (or eigen VectorXd)
     *@return double (or eigen VectorXd) the scaled value(s)
     */
    double scaled_value(double x) const {
        return (x - x_min) / (x_max - x_min);
    }

    Eigen::RowVectorXd scaled_values(Eigen::VectorXd const &x_vec) const {
        return x_vec.unaryExpr([this](double x) { return scaled_value(x); }).transpose();
    }

    double x_min;
    double x_max;

    // Spline of one-dimensional "points."
    Eigen::Spline<double, 1> spline_;
};


}  /* namespace utilities */


#endif //ES_FLOW_INTERP_H
