/*
 * trapz.h Trapezoidal Numerical Integration for Eigen
 *
 * Author:              Tom Clark  (thclark @ github)
 *
 * Copyright (c) 2019 Octue Ltd. All Rights Reserved.
 *
 */

#ifndef ES_FLOW_TRAPZ_H
#define ES_FLOW_TRAPZ_H

#include <Eigen/Dense>
#include <stdexcept>

namespace es {

using namespace Eigen;

/** @brief Trapezoidal numerical integration with unit spacing.
 *
 * Operates in the first dimension (colwise). Applied by reference or returns by value.
 *
 * @tparam Derived Type of the input array
 * @tparam OtherDerived Type of the output array
 * @param out The output (integral) array
 * @param in The input (integrand) array
 */
template<typename Derived, typename OtherDerived>
EIGEN_STRONG_INLINE void trapz(ArrayBase<OtherDerived> const & out, const ArrayBase<Derived>& in)
{
    ArrayBase<OtherDerived>& out_ = const_cast< ArrayBase<OtherDerived>& >(out);
    Array<typename ArrayBase<Derived>::Scalar, ArrayBase<Derived>::RowsAtCompileTime, ArrayBase<Derived>::ColsAtCompileTime> inter;

    if (in.rows() == 1) {
        out_.derived().setZero(1, in.cols());
    } else {
        inter = (in.topRows(in.rows()-1) + in.bottomRows(in.rows()-1)) * 0.5;
        out_.derived() = inter.colwise().sum();
    }
};

/*
 * Overload method to return result by value.
 */
template<typename Derived>
EIGEN_STRONG_INLINE Array<typename ArrayBase<Derived>::Scalar, ArrayBase<Derived>::RowsAtCompileTime, ArrayBase<Derived>::ColsAtCompileTime> trapz(const ArrayBase<Derived>& y)
{
    Array<typename ArrayBase<Derived>::Scalar, ArrayBase<Derived>::RowsAtCompileTime, ArrayBase<Derived>::ColsAtCompileTime> z;
    trapz(z,y);
    return z;
};

/** @brief Trapezoidal numerical integration with non-unit spacing
 *
 * Operates in the first dimension (colwise). Applied by reference or returns by value.
 *
 * TODO make the function also accept an array the same size as in_y so we can integrate a different spacing for
 * each column.
 *
 * @tparam DerivedX Type of the in_x array
 * @tparam DerivedY Type of the in_y array
 * @tparam DerivedOut Type of the output array
 * @param out The output (integral) array
 * @param in_x The input spacing array, must be same size as in_y
 * @param in_y The input (integrand) array
 */
template<typename DerivedX, typename DerivedY, typename DerivedOut>
EIGEN_STRONG_INLINE void trapz(ArrayBase<DerivedOut> const & out, const ArrayBase<DerivedX>& in_x, const ArrayBase<DerivedY>& in_y)
{
    // Input size check
    eigen_assert(in_x.rows() == in_y.rows());
    eigen_assert(in_x.cols() == 1);

    // Get dx for each piece of the integration
    Array<typename ArrayBase<DerivedX>::Scalar, ArrayBase<DerivedX>::RowsAtCompileTime, ArrayBase<DerivedX>::ColsAtCompileTime> dx;
    dx = (in_x.bottomRows(in_x.rows()-1) - in_x.topRows(in_x.rows()-1));


    // Get the average heights of the trapezoids
    Array<typename ArrayBase<DerivedY>::Scalar, ArrayBase<DerivedY>::RowsAtCompileTime, ArrayBase<DerivedY>::ColsAtCompileTime> inter;
    inter = (in_y.topRows(in_y.rows()-1) + in_y.bottomRows(in_y.rows()-1)) * 0.5;

    // Multiply by trapezoid widths. NB Broadcasting with *= only works for arrayX types, not arrayXX types like dx
    //  inter *= dx;
    for (int i = 0 ; i < inter.cols() ; ++i)
    {
        for (int j = 0 ; j < inter.rows() ; ++j)
        {
            inter(j,i) *= dx(j,0);
        }
    }

    // Initialise output
    ArrayBase<DerivedOut>& out_ = const_cast< ArrayBase<DerivedOut>& >(out);
    out_.derived().setZero(1, in_y.cols());

    // Output the column-wise sum
    out_.derived() = inter.colwise().sum();
};

/*
 * Overload method to return result by value.
 */
template<typename DerivedX, typename DerivedY, typename DerivedOut>
EIGEN_STRONG_INLINE Array<typename ArrayBase<DerivedOut>::Scalar, ArrayBase<DerivedOut>::RowsAtCompileTime, ArrayBase<DerivedOut>::ColsAtCompileTime> trapz(const ArrayBase<DerivedX>& x, const ArrayBase<DerivedY>& y)
{
    Array<typename ArrayBase<DerivedOut>::Scalar, ArrayBase<DerivedOut>::RowsAtCompileTime, ArrayBase<DerivedOut>::ColsAtCompileTime> z;
    trapz(z, x, y);
    return z;
};

};

#endif //ES_FLOW_TRAPZ_H


