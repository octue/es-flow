/*
 * variable_writers.h Helper functions for writing data to mat files from Eigen arrays and vectors.
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

#ifndef ES_FLOW_VARIABLE_WRITERS_H
#define ES_FLOW_VARIABLE_WRITERS_H

#include <iostream>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

#include "matio.h"

using Eigen::Array;
using Eigen::Array3d;
using Eigen::ArrayXXd;
using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::Tensor;
using Eigen::Dynamic;


typedef Eigen::Array<double, 3, 2> Array32d;


// TODO the Eigen variable writers are getting very unwieldy. Template them!


namespace es {


void writeString(mat_t *mat_fp, const std::string &var_name, std::string &var) {

    // Create the variable
    size_t dims[2] = {1, var.length()};
    matvar_t *mat_var = Mat_VarCreate(
        var_name.c_str(),
        MAT_C_CHAR,
        MAT_T_UTF8,
        2,
        dims,
        static_cast<void*>(const_cast<char*>(var.data())),
        0
    );
    if (nullptr == mat_var ) {
        throw std::runtime_error("Unable to write variable '" + var_name + "' to file.");
    }

    // Write the data and free the pointer
    Mat_VarWrite(mat_fp, mat_var, MAT_COMPRESSION_NONE);
    Mat_VarFree(mat_var);

}


void writeArray3d(mat_t *mat_fp, const std::string &var_name, const Eigen::Array3d &var) {

    // Create the variable
    size_t dims[2] = {3, 1};
    matvar_t *mat_var = Mat_VarCreate(
        var_name.c_str(),
        MAT_C_DOUBLE,
        MAT_T_DOUBLE,
        2,
        dims,
        static_cast<void*>(const_cast<double*>(var.data())),
        0
    );
    if (nullptr == mat_var ) {
        throw std::runtime_error("Unable to write variable '" + var_name + "' to file.");
    }

    // Write the data and free the pointer
    Mat_VarWrite(mat_fp, mat_var, MAT_COMPRESSION_NONE);
    Mat_VarFree(mat_var);

}


void writeArrayXd(mat_t *mat_fp, const std::string &var_name, const Eigen::ArrayXd &var) {

    // Create the variable
    size_t dims[2];
    dims[0] = var.rows();
    dims[1] = 1;
    matvar_t *mat_var = Mat_VarCreate(
        var_name.c_str(),
        MAT_C_DOUBLE,
        MAT_T_DOUBLE,
        2,
        dims,
        static_cast<void*>(const_cast<double*>(var.data())),
        0
    );
    if (nullptr == mat_var ) {
        throw std::runtime_error("Unable to write variable '" + var_name + "' to file.");
    }

    // Write the data and free the pointer
    Mat_VarWrite(mat_fp, mat_var, MAT_COMPRESSION_NONE);
    Mat_VarFree(mat_var);

}


void writeArrayXXd(mat_t *mat_fp, const std::string &var_name, const Eigen::ArrayXXd &var) {

    // Create the variable
    size_t dims[2];
    dims[0] = var.rows();
    dims[1] = var.cols();
    matvar_t *mat_var = Mat_VarCreate(
        var_name.c_str(),
        MAT_C_DOUBLE,
        MAT_T_DOUBLE,
        2,
        dims,
        static_cast<void*>(const_cast<double*>(var.data())),
        0
    );
    if (nullptr == mat_var ) {
        throw std::runtime_error("Unable to write variable '" + var_name + "' to file.");
    }

    // Write the data and free the pointer
    Mat_VarWrite(mat_fp, mat_var, MAT_COMPRESSION_NONE);
    Mat_VarFree(mat_var);

}


void writeArray32d(mat_t *mat_fp, const std::string &var_name, const Eigen::Array<double, 3, 2> &var) {

    // Create the variable
    size_t dims[2] = {3, 2};
    matvar_t *mat_var = Mat_VarCreate(
        var_name.c_str(),
        MAT_C_DOUBLE,
        MAT_T_DOUBLE,
        2,
        dims,
        static_cast<void*>(const_cast<double*>(var.data())),
        0
    );
    if (nullptr == mat_var ) {
        throw std::runtime_error("Unable to write variable '" + var_name + "' to file.");
    }

    // Write the data and free the pointer
    Mat_VarWrite(mat_fp, mat_var, MAT_COMPRESSION_NONE);
    Mat_VarFree(mat_var);

}


void writeTensor3d(mat_t *mat_fp, const std::string &var_name, const Tensor<double, 3> &var) {

    // Create the variable
    size_t dims[3];
    dims[0] = var.dimensions()[0];
    dims[1] = var.dimensions()[1];
    dims[2] = var.dimensions()[2];
    matvar_t *mat_var = Mat_VarCreate(
        var_name.c_str(),
        MAT_C_DOUBLE,
        MAT_T_DOUBLE,
        3,
        dims,
        static_cast<void*>(const_cast<double*>(var.data())),
        0
    );
    if (nullptr == mat_var ) {
        throw std::runtime_error("Unable to write variable '" + var_name + "' to file.");
    }

    // Write the data and free the pointer
    Mat_VarWrite(mat_fp, mat_var, MAT_COMPRESSION_NONE);
    Mat_VarFree(mat_var);

}


void writeDouble(mat_t *mat_fp, const std::string &var_name, const double var){

    // Create the variable
    size_t dims[2] = {1, 1};
    matvar_t *mat_var = Mat_VarCreate(
        var_name.c_str(),
        MAT_C_DOUBLE,
        MAT_T_DOUBLE,
        3,
        dims,
        static_cast<void*>(const_cast<double*>(&var)),
        0
    );
    if (nullptr == mat_var ) {
        throw std::runtime_error("Unable to write variable '" + var_name + "' to file.");
    }

    // Write the data and free the pointer
    Mat_VarWrite(mat_fp, mat_var, MAT_COMPRESSION_NONE);
    Mat_VarFree(mat_var);

}


} /* end namespace */


#endif // ES_FLOW_VARIABLE_WRITERS_H
