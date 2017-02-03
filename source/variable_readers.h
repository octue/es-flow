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

#ifndef ES_FLOW_VARIABLE_READERS_H
#define ES_FLOW_VARIABLE_READERS_H

#include <string>
#include "matio.h"
#include <eigen3/Eigen/Core>

using Eigen::Array;
using Eigen::ArrayXXd;
using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::Dynamic;

namespace es {

    void checkVariableType(matvar_t *mat_var, int matvar_type) {

        // Throw error if the requested data type does not match the saved data type
        if (mat_var->class_type != matvar_type) {
            std::string msg = "Error reading mat file (" + std::string(mat_var->name) + " incorrect variable type)";
            throw std::invalid_argument(msg);
        }

    }

    matvar_t * getVariable(mat_t *matfp, const std::string var_name, bool print_var = true) {

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

        // Crap out if rank > 2
        if (mat_var->rank > 2) {
            std::string msg = "Rank of variable " + var_name + " exceeds 2";
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

#endif //ES_FLOW_VARIABLE_READERS_H