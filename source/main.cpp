#include <stdlib.h>
#include <stdio.h>
#include "matio.h"
#include "ceres/ceres.h"
#include "glog/logging.h"
#include "mkl_dfti.h"


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

int main(int argc, char** argv) {
    google::InitGoogleLogging(argv[0]);
    // The variable to solve for with its initial value. It will be
    // mutated in place by the solver.
    double x = 0.5;
    const double initial_x = x;
    // Build the problem.
    Problem problem;
    // Set up the only cost function (also known as residual). This uses
    // auto-differentiation to obtain the derivative (jacobian).
    CostFunction* cost_function =
            new AutoDiffCostFunction<CostFunctor, 1, 1>(new CostFunctor);
    problem.AddResidualBlock(cost_function, NULL, &x);
    // Run the solver!
    Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    std::cout << summary.BriefReport() << "\n";
    std::cout << "x : " << initial_x
              << " -> " << x << "\n";

    // run an FFT

    float xf[200][100];
    DFTI_DESCRIPTOR_HANDLE fft;
    MKL_LONG len[2] = {200, 100};
    // initialize x
    DftiCreateDescriptor ( &fft, DFTI_SINGLE, DFTI_REAL, 2, len );
    DftiCommitDescriptor ( fft );
    DftiComputeForward ( fft, xf );
    DftiFreeDescriptor ( &fft );

    std::cout << "FFT COMPLETE\n";


    mat_t *matfp;

    matfp = Mat_CreateVer("matfile73.mat", NULL, MAT_FT_MAT73);
    if ( NULL == matfp ) {
        fprintf(stderr,"Error creating MAT file \"matfile73.mat\"!\n");
        return EXIT_FAILURE;
    }
    Mat_Close(matfp);

    std::cout << "MATIO TEST COMPLETE";



    return 0;
}
