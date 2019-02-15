/*
 * cumtrapz.h Cumulative Trapezoidal Numerical Integration for Eigen
 *
 * Author:              Tom Clark  (thclark @ github)
 *
 * Copyright (c) 2019 Octue Ltd. All Rights Reserved.
 *
 */

#ifndef ES_FLOW_CUMTRAPZ_H
#define ES_FLOW_CUMTRAPZ_H

#include <Eigen/Dense>


namespace es {


/** @brief Columnwise cumulative sum of elements.
 *
 * @tparam Derived Type of the input array
 * @tparam OtherDerived Type of the output array
 * @param out The output (cumulative sum) array
 * @param in The input (integrand) array
 */
template<typename Derived, typename OtherDerived>
EIGEN_STRONG_INLINE void cumsum(Eigen::ArrayBase<OtherDerived> const & out, const Eigen::ArrayBase<Derived>& in)
{
    Eigen::ArrayBase<OtherDerived>& out_ = const_cast< Eigen::ArrayBase<OtherDerived>& >(out);
    out_ = in;
    for (int i = 0 ; i < out_.cols() ; ++i)
    {
        for (int j = 1 ; j < out_.rows() ; ++j)
        {
            out_.coeffRef(j,i) += out_.coeff(j-1,i);
        }
    }
}
// Remove template specialisation from doc (causes duplicate) @cond
template<typename Derived>
EIGEN_STRONG_INLINE Eigen::Array<typename Eigen::ArrayBase<Derived>::Scalar, Eigen::ArrayBase<Derived>::RowsAtCompileTime, Eigen::ArrayBase<Derived>::ColsAtCompileTime> cumsum(const Eigen::ArrayBase<Derived>& y)
{
    Eigen::Array<typename Eigen::ArrayBase<Derived>::Scalar, Eigen::ArrayBase<Derived>::RowsAtCompileTime, Eigen::ArrayBase<Derived>::ColsAtCompileTime> z;
    cumsum(z,y);
    return z;
}
// @endcond


/** @brief Trapezoidal numerical integration with unit spacing.Cumulative trapezoidal numerical integration with unit spacing.
 *
 * Operates in the first dimension (colwise). Applied by reference or returns by value.
 *
 * @tparam Derived Type of the input array
 * @tparam OtherDerived Type of the output array
 * @param out The output (integral) array
 * @param in The input (integrand) array
 */
template<typename Derived, typename OtherDerived>
EIGEN_STRONG_INLINE void cumtrapz(Eigen::ArrayBase<OtherDerived> const & out, const Eigen::ArrayBase<Derived>& in)
{
    Eigen::ArrayBase<OtherDerived>& out_ = const_cast< Eigen::ArrayBase<OtherDerived>& >(out);
    out_.derived().setZero(in.rows(), in.cols());
    cumsum(out_.bottomRows(out_.rows()-1), (in.topRows(in.rows()-1) + in.bottomRows(in.rows()-1)) * 0.5);
}
// Remove template specialisation from doc (causes duplicate) @cond
template<typename Derived>
EIGEN_STRONG_INLINE Eigen::Array<typename Eigen::ArrayBase<Derived>::Scalar, Eigen::ArrayBase<Derived>::RowsAtCompileTime, Eigen::ArrayBase<Derived>::ColsAtCompileTime> cumtrapz(const Eigen::ArrayBase<Derived>& y)
{
    Eigen::Array<typename Eigen::ArrayBase<Derived>::Scalar, Eigen::ArrayBase<Derived>::RowsAtCompileTime, Eigen::ArrayBase<Derived>::ColsAtCompileTime> z;
    cumtrapz(z,y);
    return z;
}
// @endcond


/** @brief Cumulative trapezoidal numerical integration with non-unit spacing.
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
EIGEN_STRONG_INLINE void cumtrapz(Eigen::ArrayBase<DerivedOut> const & out, const Eigen::ArrayBase<DerivedX>& in_x, const Eigen::ArrayBase<DerivedY>& in_y)
{
    // Input size check
    eigen_assert(in_x.rows() == in_y.rows());
    eigen_assert(in_x.cols() == 1);

    // Get dx for each piece of the integration
    Eigen::Array<typename Eigen::ArrayBase<DerivedX>::Scalar, Eigen::ArrayBase<DerivedX>::RowsAtCompileTime, Eigen::ArrayBase<DerivedX>::ColsAtCompileTime> dx;
    dx = (in_x.bottomRows(in_x.rows()-1) - in_x.topRows(in_x.rows()-1));

    // Get the average heights of the trapezoids
    Eigen::Array<typename Eigen::ArrayBase<DerivedY>::Scalar, Eigen::ArrayBase<DerivedY>::RowsAtCompileTime, Eigen::ArrayBase<DerivedY>::ColsAtCompileTime> inter;
    inter = (in_y.topRows(in_y.rows()-1) + in_y.bottomRows(in_y.rows()-1)) * 0.5;

    // Multiply by trapezoid widths to give areas of the trapezoids. 
    // NB Broadcasting with *= only works for arrayX types, not arrayXX types like dx
    //  inter *= dx;
    for (int i = 0 ; i < inter.cols() ; ++i)
    {
        for (int j = 0 ; j < inter.rows() ; ++j)
        {
            inter(j,i) *= dx(j,0);
        }
    }

    // Initialise output
    Eigen::ArrayBase<DerivedOut>& out_ = const_cast< Eigen::ArrayBase<DerivedOut>& >(out);
    out_.derived().setZero(in_y.rows(), in_y.cols());

    // Perform the cumulative sum down the columns
    cumsum(out_.bottomRows(out_.rows()-1), inter);
}
// Remove template specialisation from doc (causes duplicate) @cond
template<typename DerivedX, typename DerivedY, typename DerivedOut>
EIGEN_STRONG_INLINE Eigen::Array<typename Eigen::ArrayBase<DerivedOut>::Scalar, Eigen::ArrayBase<DerivedOut>::RowsAtCompileTime, Eigen::ArrayBase<DerivedOut>::ColsAtCompileTime> cumtrapz(const Eigen::ArrayBase<DerivedX>& x, const Eigen::ArrayBase<DerivedY>& y)
{
    Eigen::Array<typename Eigen::ArrayBase<DerivedOut>::Scalar, Eigen::ArrayBase<DerivedOut>::RowsAtCompileTime, Eigen::ArrayBase<DerivedOut>::ColsAtCompileTime> z;
    cumtrapz(z, x, y);
    return z;
}

}  /* namespace utilities */

#endif //ES_FLOW_CUMTRAPZ_H


