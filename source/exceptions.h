/*
 * exceptions.h Customised exceptions for appropriate and fine grained error handling
 *
 * Author:              Tom Clark  (thclark@github)
 *
 * Copyright (c) 2017-9 Octue Ltd. All Rights Reserved.
 *
 */

#ifndef ES_FLOW_EXCEPTIONS_H
#define ES_FLOW_EXCEPTIONS_H


#include <exception>


namespace es {


struct NotImplementedException : public std::exception {
    std::string message = "Not yet implemented";
    const char *what() const throw() {
        return message.c_str();
    }
};

struct InvalidEddyTypeException : public std::exception {
    std::string message = "Unknown eddy type. Eddy type string must be one of 'A', 'B1', 'B2', 'B3' or 'B4'.";
    const char *what() const throw() {
        return message.c_str();
    }

};

}

#endif // ES_FLOW_EXCEPTIONS_H
