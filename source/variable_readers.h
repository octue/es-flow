/*
 * variable_readers.h Helper functions for reading data from mat files to Eigen arrays and vectors.
 *
 * References:
 *
 *   [1] Eigen: http://eigen.tuxfamily.org/dox/index.html
 *
 *   [2] matio: https://sourceforge.net/projects/matio/
 *
 *   [3] eigen-matio: https://github.com/tesch1/eigen-matio
 *
 * Future Improvements:
 *
 *   [1] Extension to structs and a range of different types
 *
 *   [2] More elegant implementation based on type, or possibly use of eigen-matio (ref [3])
 *
 * Author:                   Tom Clark (thclark @ github)
 *
 * Copyright (c) 2019 Octue Ltd. All Rights Reserved.
 *
 */

#ifndef ES_FLOW_VARIABLE_READERS_H
#define ES_FLOW_VARIABLE_READERS_H

#include <iostream>
#include <string>
#include "matio.h"
#include <eigen3/Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>

using Eigen::Array;
using Eigen::Array3d;
using Eigen::ArrayXXd;
using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::Tensor;
using Eigen::Dynamic;
typedef Eigen::Array<double, 3, 2> Array32d;


// TODO the Eigen variable readers are getting very unwieldy. Template them!

namespace es {


/** Return the size (i.e. number of elements) in an array variable
 * Comes from https://coderwall.com/p/nb9ngq/better-getting-array-size-in-c
 *
 * Example:
 *  int arr[] = {1, 2, 3, 4, 5};
 *  std::cout << array_size(arr) << std::endl;
 *
 */
template<size_t SIZE, class T> inline size_t array_size(T (&arr)[SIZE]) {
    return SIZE;
}

void checkVariableType(matvar_t *mat_var, int matvar_type) {

    // Throw error if the requested data type does not match the saved data type
    if (mat_var->class_type != matvar_type) {
        std::string msg = "Error reading mat file (" + std::string(mat_var->name) + " incorrect variable type)";
        throw std::invalid_argument(msg);
    }

}

/** Get pointer to a variable and check validity. Optionally print variable and check rank.
 *
 * @param matfp
 * @param var_name
 * @param print_var
 * @param max_rank The maximum rank of the variable. Default 2 as everything is an array in MATLAB.
 * @return
 */
matvar_t * getVariable(mat_t *matfp, const std::string var_name, bool print_var = true, int max_rank = 2) {

    // Get the variable's structure pointer and check validity
    matvar_t *mat_var = Mat_VarRead(matfp, var_name.c_str());
    if (mat_var == NULL) {
        std::string msg = "Error reading mat file (most likely incorrect file type: " + var_name + " variable is missing)";
        throw std::invalid_argument(msg);
    }

    // Optionally print the variable information to terminal
    if (print_var) {
        std::cout << "Reading variable: " << std::endl;
        Mat_VarPrint(mat_var, true);
    }

    // Check the rank
    if (mat_var->rank > max_rank) {
        std::string msg = "Rank of variable " + var_name + " exceeds required rank of " + std::to_string(max_rank) ;
        throw std::invalid_argument(msg);
    }

    return mat_var;
}

std::string readString(mat_t *matfp, const std::string var_name, bool print_var) {

    // Get the variable's structure pointer and check validity
    matvar_t *mat_var = getVariable(matfp, var_name, print_var);

    // Read the const char from the file, instantiate as std::string
    std::string var = std::string((const char *) mat_var->data, mat_var->dims[1]);

    // Free the data pointer and return the new variable
    Mat_VarFree(mat_var);
    return var;
}

Vector3d readVector3d(mat_t *matfp, const std::string var_name, bool print_var) {

    // Get the variable's structure pointer and check validity
    matvar_t *mat_var = getVariable(matfp, var_name, print_var);

    // Check for three elements always
    if (mat_var->dims[1]*mat_var->dims[0] != 3) {
        std::string msg = "Number of elements in variable '" + var_name + "' not equal to 3";
        throw std::invalid_argument(msg);
    }

    // Read a vector, with automatic transpose to column vector always.
    Vector3d var;
    if (mat_var->dims[0] == 1) {
        if (print_var) std::cout << "Converting row vector '" << mat_var->name << "' to column vector" << std::endl;
    }
    var = Vector3d();
    double *var_d = (double *) mat_var->data;
    long int i = 0;
    for (i = 0; i < 3; i++) {
        var[i] = var_d[i];
    }

    // Free the data pointer and return the new variable
    Mat_VarFree(mat_var);
    return var;

}

Array3d readArray3d(mat_t *matfp, const std::string var_name, bool print_var) {

    // Get the variable's structure pointer and check validity
    matvar_t *mat_var = getVariable(matfp, var_name, print_var);

    // Check for three elements always
    if (mat_var->dims[1]*mat_var->dims[0] != 3) {
        std::string msg = "Number of elements in variable '" + var_name + "' not equal to 3";
        throw std::invalid_argument(msg);
    }

    // Read a vector, with automatic transpose to column vector always.
    Eigen::Array3d var;
    if (mat_var->dims[0] == 1) {
        if (print_var) std::cout << "Converting row vector '" << mat_var->name << "' to column vector" << std::endl;
    }
    var =  Eigen::Array3d();
    double *var_d = (double *) mat_var->data;
    long int i = 0;
    for (i = 0; i < 3; i++) {
        var[i] = var_d[i];
    }

    // Free the data pointer and return the new variable
    Mat_VarFree(mat_var);
    return var;

}

VectorXd readVectorXd(mat_t *matfp, const std::string var_name, bool print_var) {

    // Get the variable's structure pointer and check validity
    matvar_t *mat_var = getVariable(matfp, var_name, print_var);

    // Read a vector, with automatic transpose to column vector always.
    VectorXd var;
    if (mat_var->dims[0] == 1) {
        if (print_var) std::cout << "Converting row vector '" << mat_var->name << "' to column vector" << std::endl;
        var = VectorXd(mat_var->dims[1], mat_var->dims[0]);
    } else {
        var = VectorXd(mat_var->dims[0], mat_var->dims[1]);
    }

    // Copy the data into the native Eigen types. However, we can also map to data already in
    // memory which avoids a copy. Read about mapping here:
    // http://eigen.tuxfamily.org/dox/group__TutorialMapClass.html
    // TODO consider mapping to reduce peak memory overhead
    double *var_d = (double *) mat_var->data;
    long int i = 0;
    for (i = 0; i < var.rows()*var.cols(); i++) {
        var[i] = var_d[i];
    }

    // Free the data pointer and return the new variable
    Mat_VarFree(mat_var);
    return var;

}


Array32d readArray32d(mat_t *matfp, const std::string var_name, bool print_var) {

    // Get the variable's structure pointer and check validity
    matvar_t *mat_var = getVariable(matfp, var_name, print_var);

    // Check for three elements always
    if ((mat_var->dims[0] != 3) || (mat_var->dims[1] != 2)) {
        std::string msg = "Variable '" + var_name + "' must be of size 3 x 2";
        throw std::invalid_argument(msg);
    }

    // Declare and size the array
    Eigen::Array<double, 3, 2> var;
    var = Eigen::Array<double, 3, 2>();

    // Copy the data into the native Eigen types. However, we can also map to data already in
    // memory which avoids a copy. Read about mapping here:
    // http://eigen.tuxfamily.org/dox/group__TutorialMapClass.html
    // TODO consider mapping to reduce peak memory overhead
    double *var_d = (double *) mat_var->data;
    long int i = 0;
    for (i = 0; i < var.rows()*var.cols(); i++) {
        var(i) = var_d[i];
    }

    // Free the data pointer and return the new variable
    Mat_VarFree(mat_var);
    return var;

}


ArrayXXd readArrayXXd(mat_t *matfp, const std::string var_name, bool print_var) {

    // Get the variable's structure pointer and check validity
    matvar_t *mat_var = getVariable(matfp, var_name, print_var);

    // Read an array. We always read as eigen arrays, not matrices, as these are far more flexible.
    // Can easily cast to matrix type later if linear algebra functionality is needed.

    ArrayXXd var;
    var = ArrayXXd(mat_var->dims[0],  mat_var->dims[1]);

    // Copy the data into the native Eigen types. However, we can also map to data already in
    // memory which avoids a copy. Read about mapping here:
    // http://eigen.tuxfamily.org/dox/group__TutorialMapClass.html
    // TODO consider mapping to reduce peak memory overhead
    double *var_d = (double *) mat_var->data;
    long int i = 0;
    long int j = 0;
    long int ind = 0;
    for (j = 0; j < var.cols(); j++) {
        for (i = 0; i < var.rows(); i++) {
            var(i, j) = var_d[ind];
            ind = ind+1;
        }
    }

    // Free the data pointer and return the new variable
    Mat_VarFree(mat_var);
    return var;

}

Tensor<double, 3> readTensor3d(mat_t *matfp, const std::string var_name, bool print_var) {

    // Get the variable's structure pointer and check validity for rank 3
    matvar_t *mat_var = getVariable(matfp, var_name, print_var, 3);

    // Read the tensor
    // TODO this code typecasts from size_t to long int, meaning possible data loss for very large arrays. Add a check.
    Eigen::Index dim0 = mat_var->dims[0];
    Eigen::Index dim1 = mat_var->dims[1];
    Eigen::Index dim2 = mat_var->dims[2];
    // std::cout  << "Initialising Tensor of size: " << dim0 << ", " << dim1 << ", "  << dim2 << std::endl;
    Eigen::Tensor<double, 3> var = Eigen::Tensor<double, 3>(dim0, dim1, dim2);
    var.setZero();

    // Copy the data into the native Eigen types. However, we can also map to data already in
    // memory which avoids a copy. Read about mapping here:
    // http://eigen.tuxfamily.org/dox/group__TutorialMapClass.html
    // TODO consider mapping to reduce peak memory overhead
    double *var_d = (double *) mat_var->data;
    long int i = 0;
    long int j = 0;
    long int k = 0;
    long int ind = 0;
    for (k = 0; k < var.dimensions()[2]; k++) {
        for (j = 0; j < var.dimensions()[1]; j++) {
            for (i = 0; i < var.dimensions()[0]; i++) {
                var(i, j, k) = var_d[ind];
                ind = ind + 1;
            }
        }
    }

    // Free the data pointer and return the new variable
    Mat_VarFree(mat_var);
    return var;

}

double readDouble(mat_t *matfp, const std::string var_name, bool print_var){

    // Get the variable's structure pointer and check validity
    matvar_t *mat_var = getVariable(matfp, var_name, print_var);

    // Read a single element double
    double *var_d = (double *) mat_var->data;
    double var = var_d[0];

    // Free the data pointer and return the new variable
    Mat_VarFree(mat_var);
    return var;
}


} /* end namespace */


#endif // ES_FLOW_VARIABLE_READERS_H
