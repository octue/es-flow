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
#include "variable_readers.h"

using Eigen::Array;
using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::Dynamic;

namespace es {

    class BasicLidar {
    public:
        const std::string type = "lidar_basic";
        VectorXd t;
        VectorXd z;
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

            // Straightforward reads
            t = readVectorXd(matfp, "t", print_var);
            z = readVectorXd(matfp, "z", print_var);
            u = readArrayXXd(matfp, "u", print_var);
            v = readArrayXXd(matfp, "v", print_var);
            w = readArrayXXd(matfp, "w", print_var);
            half_angle = readDouble(matfp, "half_angle", print_var);

            // Handle initialisation of position as a two element vector, zero padded (elevation = 0) and as a three
            // element vector.
            VectorXd pos = readVectorXd(matfp, "position", print_var);
            if (pos.size() == 2) {
                position = Vector3d(pos(0), pos(1), 0.0);
            }else {
                position = Vector3d(pos(0), pos(1), pos(2));
            }

            // TODO read in and tests on units structure

        }

    };

} /* end namespace */

#endif //ES_FLOW_DATA_TYPES_H
