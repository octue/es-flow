/*
 * READERS.CPP Reader class and specialised subclasses for OAS data files
 *
 * References:
 *
 *   [1]
 *
 * Future Improvements:
 *
 *   [1] Add further readers to expand capable file types
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
 * Copyright (c) 2016-17 Ocean Array Systems, All Rights Reserved.
 */

#include "readers.h"
#include "matio.h"
#include <stdexcept>

Reader::Reader(const std::string &file) : file(file) {

    // Open the MAT file for reading and save the pointer
    matfp = Mat_Open(file.c_str(), MAT_ACC_RDONLY);
    if ( matfp == NULL ) {
        std::string msg = "Error reading MAT file: ";
        throw std::invalid_argument(msg + file);
    }

    // Get the file type
    matvar_t *type_var = Mat_VarReadInfo(matfp, "type");
    if ( type_var == NULL ) {
        throw std::invalid_argument(
                "Error reading mat file (most likely not an OAS standard file format - 'type' variable is missing)");
    } else {
        if(type_var->class_type != MAT_C_CHAR) {
            throw std::invalid_argument("Error reading mat file ('type' variable must be a string)");
        }
        Mat_VarPrint(type_var, true);
        file_type = std::string((char*)type_var->data, type_var->dims[0]);

    }

}

Reader::~Reader() {

    // Close the file on destruction of the Reader
    Mat_Close(matfp);

}

void Reader::checkTimeseries(enum timeseries_check_level) {

    // Check the integrity of the timeseries, automatically determine dt and frequency where possible

}

LidarReader::LidarReader(const std::string &file) : Reader(file) {

    // Check the file type
    if (file_type != "lidar_basic"){
        std::string msg = "Data file of incorrect type: Expected 'lidar_basic', got: ";
        throw std::invalid_argument(msg + file_type);
    }

}
