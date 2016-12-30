/*
 * data_types.h Header defining data structures and validation methods
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
 * Copyright (c) 2016 Ocean Array Systems. All Rights Reserved.
 *
 */

#ifndef ES_FLOW_DATA_TYPES_H
#define ES_FLOW_DATA_TYPES_H

#include <string>
#include "matio.h"
#include <eigen3/Eigen/Core>

using Eigen::Array;
using Eigen::Matrix;
using Eigen::Vector3d;

namespace es {

    void checkVariableType(matvar_t *mat_var, int matvar_type) {

        // Throw error if the requested data type does not match the saved data type
        if (mat_var->class_type != matvar_type) {
            std::string msg = "Error reading mat file (" + std::string(mat_var->name) + " incorrect variable type)";
            throw std::invalid_argument(msg);
        }

    }

    void readVariable(mat_t *matfp, const std::string var_name, auto var, bool print_var) {

        // Get the variable's structure pointer and check validity
        matvar_t *mat_var = Mat_VarRead(matfp, var_name);
        if (mat_var == NULL) {
            std::string msg = "Error reading mat file (most likely incorrect file type: " + var_name + " variable is missing)";
            throw std::invalid_argument(msg);
        }


        // Get the mat_var_type required for the input
        std::cout << sprintf("%s\n",mat_var->name) << std::endl;
        std::cout << "Reading file type: " << typeid(var).name() << std::endl;
        switch (typeid(var).name()) {
            case "string" :
                    checkVariableType(mat_var, MAT_VAR_TYPE);
                    var = std::string((const char*)var->data, var->dims[1]);
                break;
            case "Array" :
                break;
        }

        // Crap out if rank > 2


        // Optionally print the variable
        if (print_var) Mat_VarPrint(var, true);

        Mat_VarFree(matvar);
        matvar = NULL;
    }
    }

    class BasicLidar {
        const std::string type = "lidar_basic";
        Array<double, Dynamic, 1> t;
        Vector<double, Dynamic, 1> z;
        Array<double, Dynamic, Dynamic> u;
        Array<double, Dynamic, Dynamic> v;
        Array<double, Dynamic, Dynamic> w;
        Vector3d position;
        double half_angle;
        struct {
            std::string t;
            std::string z;
            std::string u;
            std::string v;
            std::string w;
            std::string position;
            std::string half_angle;
        }units;

        void read(mat_t *matfp, bool print_var = true) {

            readVariable(matfp, 't', *t, print_var);
            readVariable(matfp, 'z', *z, print_var);
            readVariable(matfp, 'u', *u, print_var);
            readVariable(matfp, 'v', *v, print_var);
            readVariable(matfp, 'w', *w, print_var);
            readVariable(matfp, 'position', *postion, print_var);
            readVariable(matfp, 'half_angle', *half_angle, print_var);
            readVariable(matfp, 'units', *units, print_var);

        }

    };

} /* end namespace */

#endif //ES_FLOW_DATA_TYPES_H
