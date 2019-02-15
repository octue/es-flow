/*
 * main.cpp Command line executable for es-flow
 *
 * Author:              Tom Clark  (thclark @ github)
 *
 * Copyright (c) 2016-9 Octue Ltd. All Rights Reserved.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdbool.h>
#include "glog/logging.h"
#include "cxxopts.hpp"


int main(int argc, char* argv[]) {

    try {

        // Handle the input parsing and create the program help page

        // Set the program name (for --help option display) and default option behaviour
        cxxopts::Options options("es-flow", "EnvironmentSTUDIO flow library wrapper");
        bool logging = false;

        // Define the command line options, arranged visually (in the --help output) in two groups:
        options.add_options()
                ("l,log-file",
                 "Switch on logging, optionally specify the directory to save logfiles.",
                 cxxopts::value<std::string>()->implicit_value("logs"), "FILE")
                ("h,help",
                 "Display program help.",
                 cxxopts::value<bool>(), "BOOL");

        options.add_options("Input / output file")
                ("i,input-file", "Name of the input file (read only).", cxxopts::value<std::string>(), "FILE")
                ("o,output-file",
                 "Name of the output results file. Warning - this file will be overwritten if it already exists.",
                 cxxopts::value<std::string>(), "FILE");

        // Parse the input options
        options.parse(argc, argv);

        if (options.count("help")) {
            std::cout << options.help({"", "Input / output file"}) << std::endl;
            exit(0);
        }

        if (options.count("l")) {
            logging = true;
        }

        if (logging) {
            FLAGS_logtostderr = false;
            FLAGS_minloglevel = 0;
            FLAGS_log_dir = options["l"].as<std::string>();
            std::cout << "Logging to: " << FLAGS_log_dir << std::endl;
            google::InitGoogleLogging(argv[0]);
        }


    } catch (const cxxopts::OptionException &e) {
        std::cout << "Error parsing options: " << e.what() << std::endl;
        exit(1);

    } catch (...) {
        // Handle logging of general exceptions
        auto eptr = std::current_exception();
        try {
            // This deliberately rethrows the caught exception in the current scope, so we can log it
            if (eptr) {
                std::rethrow_exception(eptr);
            }
        } catch(const std::exception& e) {
            std::cout << "Caught exception: " << e.what() << std::endl;
        }
        exit(1);

    }

    exit(0);

}
