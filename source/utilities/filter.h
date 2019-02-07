/*
 * filter.h One-dimensional polynomial based digital filter for eigen arrayXd
 *
 * Author:              Tom Clark  (thclark @ github)
 *
 * Copyright (c) 2019 Octue Ltd. All Rights Reserved.
 *
 */

#ifndef ES_FLOW_FILTER_H
#define ES_FLOW_FILTER_H

#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <stdexcept>


namespace utilities {


/** @brief One-dimensional digital filter
 *
 * Filters the input 1d array (or std::vector<>) `x` with the filter described by coefficient arrays `b` and `a`.
 *
 * See https://en.wikipedia.org/wiki/Digital_filter
 *
 * If a(0) is not equal to 1, the filter coefficients are normalised by a(0).
 *
 *   a(0)*y(n-1) = b(0)*x(n-1) + b(1)*x(n-2) + ... + b(nb)*x(n-nb-1)  - a(1)*y(n-2) - ... - a(na)*y(n-na-1)
 *
 */
template<typename DerivedOut, typename DerivedIn>
void filter(Eigen::ArrayBase<DerivedOut> &y,
            const Eigen::ArrayBase<DerivedIn> &b,
            const Eigen::ArrayBase<DerivedIn> &a,
            const Eigen::ArrayBase<DerivedIn> &x) {

    // Check for 1D inputs
    eigen_assert((a.cols() == 1) || (a.rows() == 1));
    eigen_assert((b.cols() == 1) || (a.rows() == 1));
    eigen_assert((x.cols() == 1) || (a.rows() == 1));

    // Initialise output
    Eigen::ArrayBase<DerivedOut> &y_ = const_cast< Eigen::ArrayBase<DerivedOut> & >(y);
    y_.derived().setZero(x.rows(), x.cols());

    // Direct Form II transposed standard difference equation
    typename DerivedOut::Scalar tmp;
    Eigen::Index i;
    Eigen::Index j;
    for (i = 0; i < x.size(); i++) {
        tmp = 0.0;
        j = 0;
        y_(i) = 0.0;

        for (j = 0; j < b.size(); j++) {
            if (i - j < 0) continue;
            tmp += b(j) * x(i - j);
        }

        for (j = 1; j < a.size(); j++) {
            if (i - j < 0) continue;
            tmp -= a(j) * y_(i - j);
        }

        tmp /= a(0);
        y_(i) = tmp;
    }
}

template<typename T>
void filter(std::vector<T> &y, const std::vector<T> &b, const std::vector<T> &a, const std::vector<T> &x) {

    y.resize(0);
    y.resize(x.size());

    for (int i = 0; i < x.size(); i++) {
        T tmp = 0.0;
        int j = 0;
        y[i] = 0.0;
        for (j = 0; j < b.size(); j++) {
            if (i - j < 0) continue;
            tmp += b[j] * x[i - j];
        }

        for (j = 1; j < a.size(); j++) {
            if (i - j < 0) continue;
            tmp -= a[j] * y[i - j];
        }

        tmp /= a[0];
        y[i] = tmp;
    }
}

/** @brief One-dimensional deconvolution
 *
 *    `z = deconv(b, a)` deconvolves vector `a` out of column `b`, returning a column vector
 *
 */
template<typename DerivedIn, typename DerivedOut>
void deconv(Eigen::ArrayBase<DerivedOut> &z,
            const Eigen::ArrayBase<DerivedIn> &b,
            const Eigen::ArrayBase<DerivedIn> &a) {

    eigen_assert(a(0) != 0);
    eigen_assert((a.cols() == 1) || (a.rows() == 1));
    eigen_assert((b.cols() == 1) || (b.rows() == 1));

    auto na = a.size();
    auto nb = b.size();

    // Initialise output
    Eigen::ArrayBase<DerivedOut> &z_ = const_cast< Eigen::ArrayBase<DerivedOut> & >(z);

    // Cannot deconvolve a longer signal out of a smaller one
    if (na > nb) {
        z_.derived().setZero(1, 1);
        return;
    }
    z_.derived().setZero(nb - na + 1, 1);

    // Deconvolution and polynomial division are the same operations as a digital filter's impulse response B(z)/A(z):
    Eigen::Array<typename DerivedOut::Scalar, Eigen::Dynamic, 1> impulse;
    impulse.setZero(nb - na + 1, 1);
    impulse(0) = 1;
    filter(z_, b, a, impulse);
}

}  /* namespace utilities */

#endif //ES_FLOW_FILTER_H
