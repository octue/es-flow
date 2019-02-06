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
#include "matio.h"
#include "ceres/ceres.h"
#include "glog/logging.h"
#include "mkl_dfti.h"
#include "cxxopts.hpp"

#include <Eigen/Core>
#include <unsupported/Eigen/FFT>


using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;

// A templated cost functor that implements the residual r = 10 -
// x. The method operator() is templated so that we can then use an
// automatic differentiation wrapper around it to generate its
// derivatives.
struct CostFunctor {
    template <typename T> bool operator()(const T* const x, T* residual) const {
        residual[0] = T(10.0) - x[0];
        return true;
    }
};


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

        // Demonstrate MAT file read/write

        // Create the output file
        mat_t *matfp;
        matfp = Mat_CreateVer(options["o"].as<std::string>().c_str(), NULL, MAT_FT_MAT73);
        if ( NULL == matfp ) {
            std::string msg = "Error creating MAT file: ";
            throw  msg + options["o"].as<std::string>();
        }

        // We haven't written any data - just close the file
        Mat_Close(matfp);

        std::cout << "MATIO TEST COMPLETE" << std::endl;

        // Run the ceres solver

        // Define the variable to solve for with its initial value. It will be mutated in place by the solver.
        double x = 0.5;
        const double initial_x = x;

        // Build the problem.
        Problem problem;

        // Set up the only cost function (also known as residual). This uses
        // auto-differentiation to obtain the derivative (jacobian).
        CostFunction* cost_function = new AutoDiffCostFunction<CostFunctor, 1, 1>(new CostFunctor);
        problem.AddResidualBlock(cost_function, NULL, &x);

        // Run the solver!
        Solver::Options ceroptions;
        ceroptions.minimizer_progress_to_stdout = true;
        Solver::Summary summary;
        ceres::Solve(ceroptions, &problem, &summary);
        std::cout << summary.BriefReport() << "\n";
        std::cout << "x : " << initial_x << " -> " << x << "\n";

//        // Run an example FFT
//
//        // Make meaningless data and a size vector
//        float xf[200][100];
//        MKL_LONG len[2] = {200, 100};
//
//        // Create a decriptor, which is a pattern for what operation an FFT will undertake
//        DFTI_DESCRIPTOR_HANDLE fft;
//        DftiCreateDescriptor ( &fft, DFTI_SINGLE, DFTI_REAL, 2, len );
//        DftiCommitDescriptor ( fft );
//
//        // Compute a forward transform, in-place on the data
//        DftiComputeForward ( fft, xf );
//
//        // Free the descriptor
//        DftiFreeDescriptor ( &fft );
//
//        std::cout << "FFT COMPLETE\n";

        size_t dim_x = 28, dim_y = 126;
        Eigen::FFT<float> fft;
        Eigen::MatrixXf in = Eigen::MatrixXf::Random(dim_x, dim_y);
        Eigen::MatrixXcf out;
        out.setZero(dim_x, dim_y);

        for (int k = 0; k < in.rows(); k++) {
            Eigen::VectorXcf tmpOut(dim_x);
            fft.fwd(tmpOut, in.row(k));
            out.row(k) = tmpOut;
        }

        for (int k = 0; k < in.cols(); k++) {
            Eigen::VectorXcf tmpOut(dim_y);
            fft.fwd(tmpOut, out.col(k));
            out.col(k) = tmpOut;
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
