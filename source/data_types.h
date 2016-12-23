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
#include <Eigen/Core>

using Eigen::Array
using Eigen::Matrix
using Eigen::Vector

namespace es {

    class BasicLidar {
        const std::string type = "lidar_basic";
        Array<double, Dynamic, 1> t;
        Vector<double, Dynamic, 1> z;
        Array<double, Dynamic, Dynamic> u;
        Array<double, Dynamic, Dynamic> v;
        Array<double, Dynamic, Dynamic> w;
        Vector<double, 3, 1> position;
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
    };

} /* end namespace */

#endif //ES_FLOW_DATA_TYPES_H
