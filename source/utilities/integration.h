/*
 * integration.h Convenient wrappers around tbs1980's Numerical Integration module for Eigen
 *
 * Author:              Tom Clark  (thclark @ github)
 *
 * Copyright (c) 2019 Octue Ltd. All Rights Reserved.
 *
 */

#ifndef ES_FLOW_INTEGRATION_H
#define ES_FLOW_INTEGRATION_H

#include <Eigen/Dense>
#include <Eigen/Core>
#include "NumericalIntegration.h"


namespace utilities {


/** @brief Cumulative numerical integration of a functor in 1D, using gauss-Kronrod61 integration
 *
 * Operates in the first dimension (colwise). Applied by reference or returns by value.
 *
 * @tparam T_functor Type of the functor class instance supplied to the integrator
 * @param[in] x The domain of x over which to do a cumulative integration. Output values of the cumulative integral are given at these x locations.
 * @param[in] functor The integrand functor. This functor must, given a scalar value of x, return a scalar of the same type containing the value of the integrand function at the input point.
 * @param[in] max_subintervals The maximum number of subintervals allowed in the subdivision process of quadrature functions, for each segment in x. This corresponds to the amount of memory allocated for said functions.
 * @return The cumulative integration of f with respect to x, over the domain in input x
 */
template<typename T_functor>
Eigen::ArrayXd cumulative_integrate(const Eigen::ArrayXd &x, const T_functor &functor, const int max_subintervals=200)
{
    typedef double ScalarType;
    Eigen::ArrayXd integral(x.rows());
    integral.setZero();

    // Define the integrator and quadrature rule
    Eigen::Integrator<ScalarType> eigIntgtor(max_subintervals);
    Eigen::Integrator<ScalarType>::QuadratureRule quadratureRule = Eigen::Integrator<ScalarType>::GaussKronrod61;

    // Define the desired absolute and relative errors
    auto desAbsErr = ScalarType(0.0);
    ScalarType desRelErr = Eigen::NumTraits<ScalarType>::epsilon() * 50.;

    // Integrate to each value in eta
    for (Eigen::Index i = 0; i < x.rows()-1; i++) {
        integral[i+1] = eigIntgtor.quadratureAdaptive(functor, ScalarType(x[i]), ScalarType(x[i+1]), desAbsErr, desRelErr, quadratureRule) + integral[i];
    }
    return integral;
}


}  /* namespace utilities */

#endif //ES_FLOW_INTEGRATION_H


