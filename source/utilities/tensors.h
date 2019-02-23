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

#include <unsupported/Eigen/CXX11/Tensor>
#include <iostream>
#include <string>


namespace utilities {


template<typename T>
using  MatrixType = Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic>;


template<typename T>
using  ArrayType = Eigen::Array<T,Eigen::Dynamic, Eigen::Dynamic>;


/** @brief convert Eigen Tensor to a Matrix
 *
 * @code
 *    int main () {
 *        Eigen::Tensor<double,4> my_rank4 (2,2,2,2);
 *        my_rank4.setRandom();
 *
 *        Eigen::MatrixXd         mymatrix =  Tensor_to_Matrix(my_rank4, 4,4);
 *        Eigen::Tensor<double,3> my_rank3 =  Matrix_to_Tensor(mymatrix, 2,2,4);
 *
 *        std::cout << my_rank3 << std::endl;
 *
 *        return 0;
 *    }
 * @endcode
 *
 * @return
 */
template<typename Scalar,int rank, typename sizeType>
auto tensor_to_matrix(const Eigen::Tensor<Scalar,rank> &tensor,const sizeType rows,const sizeType cols)
{
    return Eigen::Map<const MatrixType<Scalar>> (tensor.data(), rows,cols);
}


/** @brief convert Eigen Matrix<> to a Tensor
 *
 * @code
 *    int main () {
 *        Eigen::Tensor<double,4> my_rank4 (2,2,2,2);
 *        my_rank4.setRandom();
 *
 *        Eigen::MatrixXd         mymatrix =  Tensor_to_Matrix(my_rank4, 4,4);
 *        Eigen::Tensor<double,3> my_rank3 =  Matrix_to_Tensor(mymatrix, 2,2,4);
 *
 *        std::cout << my_rank3 << std::endl;
 *
 *        return 0;
 *    }
 * @endcode
 *
 * @return
 */
template<typename Scalar, typename... Dims>
auto matrix_to_tensor(const MatrixType<Scalar> &matrix, Dims... dims)
{
    constexpr int rank = sizeof... (Dims);
    return Eigen::TensorMap<Eigen::Tensor<const Scalar, rank>>(matrix.data(), {dims...});
}


/** @brief Convert Eigen Tensor<> to an Array<>
 *
 * @code
 *    int main () {
 *        Eigen::Tensor<double,4> my_rank4 (2,2,2,2);
 *        my_rank4.setRandom();
 *
 *        Eigen::ArrayXd         myarray =  Tensor_to_Array(my_rank4, 4,4);
 *        Eigen::Tensor<double,3> my_rank3 =  Array_to_Tensor(mymatrix, 2,2,4);
 *
 *        std::cout << my_rank3 << std::endl;
 *
 *        return 0;
 *    }
 * @endcode
 *
 * @return
 */
template<typename Scalar,int rank, typename sizeType>
auto tensor_to_array(const Eigen::Tensor<Scalar,rank> &tensor,const sizeType rows,const sizeType cols)
{
    return Eigen::Map<const ArrayType<Scalar>> (tensor.data(), rows,cols);
}


/** @brief convert Eigen Tensor<> to an Array<>
 *
 * @code
 *    int main () {
 *        Eigen::Tensor<double,4> my_rank4 (2,2,2,2);
 *        my_rank4.setRandom();
 *
 *        Eigen::ArrayXd         myarray =  Tensor_to_Array(my_rank4, 4,4);
 *        Eigen::Tensor<double,3> my_rank3 =  Array_to_Tensor(mymatrix, 2,2,4);
 *
 *        std::cout << my_rank3 << std::endl;
 *
 *        return 0;
 *    }
 * @endcode
 *
 * @return
 */
template<typename Scalar, typename... Dims>
auto array_to_tensor(const ArrayType<Scalar> &matrix, Dims... dims)
{
    constexpr int rank = sizeof... (Dims);
    return Eigen::TensorMap<Eigen::Tensor<const Scalar, rank>>(matrix.data(), {dims...});
}


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

} /* namespace utilities */

#endif //ES_FLOW_TENSORS_H
