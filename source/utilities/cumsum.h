/*
 * cumsum.h Cumulative sum of a vector
 *
 * Author:              Tom Clark  (thclark @ github)
 *
 * Copyright (c) 2019 Octue Ltd. All Rights Reserved.
 *
 * ---
 *
 * This file is adapted from part of libigl, a simple c++ geometry processing library.
 *
 * Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public License
 * v. 2.0. If a copy of the MPL was not distributed with this file, You can
 * obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#ifndef ES_FLOW_CUMSUM_H
#define ES_FLOW_CUMSUM_H

#include <Eigen/Core>
#include <Eigen/Dense>


namespace utilities {

/* Computes a cumulative sum down the columns (or across the rows) of X
 *
 * @tparam DerivedX Type of matrix/array X
 * @tparam DerivedY Type of matrix/array Y
 * @param[in] X    m by n Matrix to be cumulatively summed.
 * @param[in] dim  dimension to take cumulative sum (1 goes down columns, 2 across rows). Default 1.
 * @param[out] Y   m by n Matrix containing cumulative sum.
 */
template <typename DerivedX, typename DerivedY>
EIGEN_STRONG_INLINE void cumsum(Eigen::PlainObjectBase<DerivedY > & y, const Eigen::MatrixBase<DerivedX > & x, const int dim=1) {

    y.resizeLike(x);

    Eigen::Index num_outer = (dim == 1 ? x.cols() : x.rows() );
    Eigen::Index num_inner = (dim == 1 ? x.rows() : x.cols() );

    // This has been optimized so that dim = 1 or 2 is roughly the same cost.
    // (Optimizations assume ColMajor order)
    if (dim == 1) {
        // TODO multithread the outer loop
        for(int o = 0; o<num_outer; o++)
        {
            typename DerivedX::Scalar sum = 0;
            for (int i = 0;i<num_inner;i++) {
                sum += x(i,o);
                y(i,o) = sum;
            }
        }
    } else {
        for (int i = 0; i<num_inner; i++) {
            // TODO multithread the inner loop
            for (int o = 0; o<num_outer; o++) {
                if (i == 0) {
                    y(o,i) = x(o,i);
                } else {
                    y(o,i) = y(o,i-1) + x(o,i);
                }
            }
        }
    }
}

}

#endif //ES_FLOW_CUMSUM_H
