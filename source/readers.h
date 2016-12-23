//
// Created by Thomas Clark on 20/12/2016.
//

#ifndef ES_FLOW_READERS_H
#define ES_FLOW_READERS_H

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdbool.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include "matio.h"
#include "glog/logging.h"
#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>
#include <stdexcept>

namespace es {

    // Level of data validation applied to timeseries. Exception is thrown if timeseries checks fail
    enum timeseries_check_level {
        PRESENT             = 0,    // Checks that data is present and of correct type. Extracts start and end timestamps.
        INCREASING          = 1,    // As PRESENT + checks that the timestamp always increases.
        MONOTONIC           = 2,    // As INCREASING + checks that timestamp uniformly increases, allows data skip (e.g. instrument turned off then back on later). Extracts sampling frequency.
        STRICTLY_MONOTONIC  = 3,    // AS INCREASING + checks that timestamp uniformly increases, data skipping causes errors. Extracts sampling frequency.
    };

    template <typename DataType>
    class Reader {
    public:

        Reader(const std::string &file);

        ~Reader();

        void checkFileType();

        std::string logString() const;

        void read(<DataType> *data);

        void readWindow(const int index, <DataType> *data);

        void checkTimeseries(enum timeseries_check_level);

        double getWindowDuration() const;

        void setWindowDuration(double windowDuration);

        int getWindowSize() const;

        void setWindowSize(int windowSize);

    protected:

        std::string file;
        mat_t *matfp;
        std::string file_type;
        int windowSize;
        double windowDuration;
        double dt;
        std::vector<double> t;

    };

    template <class DataType>
    Reader<DataType>::Reader(const std::string &file) : file(file) {

        // Open the MAT file for reading and save the pointer
        matfp = Mat_Open(file.c_str(), MAT_ACC_RDONLY);
        if (matfp == NULL) {
            std::string msg = "Error reading MAT file: ";
            throw std::invalid_argument(msg + file);
        }

    }

    template <class DataType>
    Reader<DataType>::~Reader() {

        // Close the file on destruction of the Reader
        Mat_Close(matfp);

    }

    template <class DataType>
    void Reader<DataType>::checkFileType(){

        // Get the file type from the reserved 'type' variable in the mat file
        matvar_t *type_var = Mat_VarRead(matfp, "type");
        if (type_var == NULL) {
            throw std::invalid_argument(
                "Error reading mat file (most likely not an OAS standard file format - 'type' variable is missing)");
        } else {
            if (type_var->class_type != MAT_C_CHAR) {
                throw std::invalid_argument("Error reading mat file ('type' variable must be a character array)");
            }
            // Mat_VarPrint(type_var, true);
            file_type = std::string((const char*)type_var->data, type_var->dims[1]);
        }

    }

    template <class DataType>
    void Reader<DataType>::checkTimeseries(enum timeseries_check_level) {

        // Check the integrity of the timeseries, automatically determine dt and frequency where possible

    }

    template <class DataType>
    double Reader<DataType>::getWindowDuration() const {
        return windowDuration;
    }

    template <class DataType>
    void Reader<DataType>::setWindowDuration(double windowDuration) {
        Reader<DataType>::windowDuration = windowDuration;
    }

    template <class DataType>
    int Reader::getWindowSize() const {
        return windowSize;
    }

    template <class DataType>
    void Reader<DataType>::setWindowSize(int windowSize) {
        Reader<DataType>::windowSize = windowSize;
    }

    template <class DataType>
    std::string Reader<DataType>::logString() const {
        return std::string("Object Reader for type ") + file_type + std::string(", attached to file file ") + file;
    }

    template <class DataType>
    Reader<DataType>::read(<DataType> *data) {

    }

    template <class DataType>
    Reader<DataType>::readWindow(const int index, <DataType> *data) {

    }

    // Represent Reader classes in logs or ostream
    template <class DataType>
    ::std::ostream& operator<<(::std::ostream& os, const Reader<DataType>& reader) {
    // Represent in logs or ostream
    return os << reader.logString();
}

} /* namespace es */



#endif //ES_FLOW_READERS_H
