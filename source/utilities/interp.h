/*
 * interp.h Interpolation tools
 * 
 * Author:              Tom Clark  (thclark @ github)
 *
 * Copyright (c) 2019 Octue Ltd. All Rights Reserved.
 *
 */

#ifndef ES_FLOW_INTERP_H
#define ES_FLOW_INTERP_H

#include <iostream>
#include <math.h>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/Splines>


namespace utilities {


class LinearInterpolant {
public:
    /** @brief Create a linear interpolant
     *
     * Initialise an interpolant as follows:
     *
     *    Eigen::ArrayXd xvals(3);
     *    Eigen::ArrayXd yvals(xvals.rows());
     *
     *    xvals << 0, 15, 30;
     *    yvals << 0, 12, 17;
     *
     *    LinearInterpolant s(xvals, yvals);
     *
     */
    LinearInterpolant(const Eigen::ArrayXd x_vec, const Eigen::ArrayXd y_vec)
      : x_min(x_vec.minCoeff()), x_max(x_vec.maxCoeff()), x(x_vec), y(y_vec) {}

    /** @brief Evaluate interpolant at values xi whose yi values are unknown
     *
     * Out of range xi values are constrained to the endpoints (i.e. nearest neighbour interpolation)
     *
     * Performs linear interpolation to evaluate values yi at
     *    // ... Initialise interpolant (see constructor) ...
     *    LinearInterpolant s(xvals, yvals);
     *    std::cout << s(12.34) << std::endl;
     *
     *@param xi double (or eigen ArrayXd) of target x values for the interpolation
     *@return double interpolated value yi corresponding to location xi
     */
    double operator()(double xi) const {
        double bounded_xi = std::max(std::min(xi, x_max), x_min);
        Eigen::Index node = findNode(bounded_xi);
        if (node == x.size()-1) {
            return y[node];
        }
        double proportion = (bounded_xi - x[node]) / (x[node + 1] - x[node]);
        return proportion * (y[node + 1] - y[node]) + y[node];
    }

    Eigen::ArrayXd operator()(Eigen::ArrayXd const &xi_vec) {
        Eigen::ArrayXd yi_vec(xi_vec.rows());
        yi_vec = xi_vec.unaryExpr([this](double xi) { return this->operator()(xi); }).transpose();
        return yi_vec;
    }

private:

    // Find the index of the data value in x_vec immediately before the required xi value
    Eigen::Index findNode(const double bounded_xi) const {

        // Run a binary search to find the adjacent node in the grid. This is O(log(N)) to evaluate one location.
        // https://stackoverflow.com/questions/6553970/find-the-first-element-in-a-sorted-array-that-is-greater-than-the-target
        Eigen::Index low = 0;
        Eigen::Index high = x.size();
        while (low != high) {
            Eigen::Index mid = (low + high) / 2; // Or a fancy way to avoid int overflow
            if (x[mid] <= bounded_xi) {
                low = mid + 1;
            } else {
                high = mid;
            }
        }
        return low - 1;
    }

    // Boundaries of the vector
    double x_min;
    double x_max;

    // Hold a copy of the data in case it changes outside during the life of the interpolant
    Eigen::ArrayXd x;
    Eigen::ArrayXd y;
};


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
