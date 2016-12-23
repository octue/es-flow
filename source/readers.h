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

namespace es {

    // Level of data validation applied to timeseries. Exception is thrown if timeseries checks fail
    enum timeseries_check_level {
        PRESENT             = 0,    // Checks that data is present and of correct type. Extracts start and end timestamps.
        INCREASING          = 1,    // As PRESENT + checks that the timestamp always increases.
        MONOTONIC           = 2,    // As INCREASING + checks that timestamp uniformly increases, allows data skip (e.g. instrument turned off then back on later). Extracts sampling frequency.
        STRICTLY_MONOTONIC  = 3,    // AS INCREASING + checks that timestamp uniformly increases, data skipping causes errors. Extracts sampling frequency.
    };

    class Reader {
    public:
        Reader(const std::string &file);

        ~Reader();

        void checkTimeseries(enum timeseries_check_level);

        void checkFileType();

        std::string logString() const;

    protected:
        std::string file;
        mat_t *matfp;
        std::string file_type;

    };

    // Represent Reader classes in logs or ostream
    ::std::ostream& operator<<(::std::ostream& os, const Reader& reader);

    class LidarReader : public Reader {
    public:
        LidarReader(const std::string &file);

        void read();

        double getWindowDuration() const;

        void setWindowDuration(double windowDuration);

    protected:
        double dt;
        std::vector<double> t;
        int windowSize;
        double windowDuration;

    public:
        int getWindowSize() const;

        void setWindowSize(int windowSize);

    };

} /* namespace es */



#endif //ES_FLOW_READERS_H
