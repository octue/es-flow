/*
 * tensors.h Utilities to help with Eigen::Tensors
 *
 * Author:              Tom Clark  (thclark @ github)
 *
 * Copyright (c) 2019 Octue Ltd. All Rights Reserved.
 *
 */

#ifndef ES_FLOW_TENSORS_H
#define ES_FLOW_TENSORS_H

#include <string>


/** @brief Return a string representation of tensor dimensions
 *
 * @param tensor An Eigen::Tensor
 * @return string The output string like "[2 x 3 x 4]" for a rank 3 tensor
 */
template<typename T>
std::string tensor_dims(T &tensor) {
    std::stringstream dims;
    for (auto i = tensor.dimensions().begin(); i != tensor.dimensions().end(); ++i) {
        dims << *i << " x ";
    }
    std::string dim_str = dims.str();
    dim_str.pop_back();
    dim_str.pop_back();
    dim_str.pop_back();
    return dim_str;
}


#endif //ES_FLOW_TENSORS_H
