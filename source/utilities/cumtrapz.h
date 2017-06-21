/*
 * cumtrapz.h Brief overview sentence
 *
 * References:
 *
 *   [1]
 *
 * Future Improvements:
 *
 *   [1] 
 *
 * Author:                   T. Clark
 * Work address:             Ocean Array Systems Ltd
 *                           Hauser Forum
 *                           3 Charles Babbage Road
 *                           Cambridge
 *                           CB3 0GT
 * Email:                    tom.clark@oceanarraysystems.com
 * Website:                  www.oceanarraysystems.com
 *
 * Copyright (c) 2017 Ocean Array Systems. All Rights Reserved.
 *
 */

#ifndef ES_FLOW_CUMTRAPZ_H
#define ES_FLOW_CUMTRAPZ_H

#include <Eigen/Dense>

namespace es
{

    using namespace Eigen;

    /**
    * Cumulative sum of elements.
    */
    template<typename Derived, typename OtherDerived>
    EIGEN_STRONG_INLINE void cumsum(MatrixBase<OtherDerived> const & out, const MatrixBase<Derived>& in)
    {
        MatrixBase<OtherDerived>& out_ = const_cast< MatrixBase<OtherDerived>& >(out);
        out_ = in;
        for (int i = 0 ; i < out_.cols() ; ++i)
        {
            for (int j = 1 ; j < out_.rows() ; ++j)
            {
                out_.coeffRef(j,i) += out_.coeff(j-1,i);
            }
        }
    };

    /**
    * Overload method to return result by value.
    */
    template<typename Derived>
    EIGEN_STRONG_INLINE Matrix<typename MatrixBase<Derived>::Scalar, MatrixBase<Derived>::RowsAtCompileTime, MatrixBase<Derived>::ColsAtCompileTime>
    cumsum(const MatrixBase<Derived>& y)
    {
        Matrix<typename MatrixBase<Derived>::Scalar, MatrixBase<Derived>::RowsAtCompileTime, MatrixBase<Derived>::ColsAtCompileTime> z;
        cumsum(z,y);
        return z;
    };

    /*
    * Cumulative trapezoidal numerical integration with unit spacing.
    */
    template<typename Derived, typename OtherDerived>
    EIGEN_STRONG_INLINE void cumtrapz(MatrixBase<OtherDerived> const & out, const MatrixBase<Derived>& in)
    {
        MatrixBase<OtherDerived>& out_ = const_cast< MatrixBase<OtherDerived>& >(out);
        out_.derived().setZero(in.rows(), in.cols());
        cumsum(out_.middleRows(1,out_.rows()-1), (in.topRows(in.rows()-1) + in.bottomRows(in.rows()-1)) * 0.5);
    };

    /**
    * Overload method to return result by value.
    */
    template<typename Derived>
    EIGEN_STRONG_INLINE Matrix<typename MatrixBase<Derived>::Scalar, MatrixBase<Derived>::RowsAtCompileTime, MatrixBase<Derived>::ColsAtCompileTime>
    cumtrapz(const MatrixBase<Derived>& y)
    {
        Matrix<typename MatrixBase<Derived>::Scalar, MatrixBase<Derived>::RowsAtCompileTime, MatrixBase<Derived>::ColsAtCompileTime> z;
        cumtrapz(z,y);
        return z;
    };

    /*
     * Cumulative trapezoidal numerical integration with non-unit spacing. Operates in the first dimension
     */
    template<typename Derived, typename OtherDerived>
    EIGEN_STRONG_INLINE void cumtrapz(MatrixBase<OtherDerived> const & out, const MatrixBase<Derived>& in_x, const MatrixBase<Derived>& in_y)
    {
        // Input size check
        eigen_assert(in_x.rows() == in_y.rows());
        eigen_assert(in_x.cols() == in_y.cols());

        // Instantiate the dx matrix and set zero
        MatrixBase<Derived>& dx;
        dx.derived().setZero(in_y.rows(), in_y.cols());


        MatrixBase<OtherDerived>& out_ = const_cast< MatrixBase<OtherDerived>& >(out);
        out_.derived().setZero(in_y.rows(), in_y.cols());
        cumsum(out_.middleRows(1,out_.rows()-1), (in_y.topRows(in_y.rows()-1) + in_y.bottomRows(in_y.rows()-1)) * 0.5);

    };

    /**
    * Overload method to return result by value.
    */
    template<typename Derived>
    EIGEN_STRONG_INLINE Matrix<typename MatrixBase<Derived>::Scalar, MatrixBase<Derived>::RowsAtCompileTime, MatrixBase<Derived>::ColsAtCompileTime>
    cumtrapz(const MatrixBase<Derived>& y)
    {
        Matrix<typename MatrixBase<Derived>::Scalar, MatrixBase<Derived>::RowsAtCompileTime, MatrixBase<Derived>::ColsAtCompileTime> z;
        cumtrapz(z,y);
        return z;
    };

};


#endif // __btkEigenCumtrapz_h
#endif //ES_FLOW_CUMTRAPZ_H


