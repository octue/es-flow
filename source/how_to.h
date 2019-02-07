/*
 * how_to.h  Chunks of code that show us how to do things
 *
 *  May not be complete, may not compile. Not included in any build targets. Just copy, paste, adapt to our own uses.
 *
 * Author:              Tom Clark  (thclark @ github)
 *
 * Copyright (c) 2019 Octue Ltd. All Rights Reserved.
 *
 */


/***********************************************************************************************************************
 * HOW TO INTEGRATE USING THE NUMERICALINTEGRATION LIBRARY
 **********************************************************************************************************************/
//
#include <Eigen/Dense>
#include <Eigen/Core>
#include <iostream>
#include "NumericalIntegration.h"

/* Integrator for the velocity deficit
 *  We consider the example from:
 *
 *       http://www.gnu.org/software/gsl/manual/html_node/Numerical-integration-examples.html
 *
 *       int_0^1 x^{-1/2} log(x) dx = -4
 *
 *  The integrator expects the user to provide a functor as shown below.
 */

template<typename Scalar>
class IntegrandExampleFunctor
{
public:
    IntegrandExampleFunctor(const Scalar alpha):m_alpha(alpha)
    {
        assert(alpha>0);
    }

    Scalar operator()(const Scalar x) const
    {
        assert(x>0);
        return log(m_alpha*x) / sqrt(x);
    }

    void setAlpha(const Scalar alpha)
    {
        m_alpha = alpha;
    }
private:
    Scalar m_alpha;
};

double do_something(const double arg)
{
    // Define the scalar
    typedef double Scalar;

    // Define the functor
    Scalar alpha=1.;
    IntegrandExampleFunctor<Scalar> inFctr(alpha);

    //define the integrator
    Eigen::Integrator<Scalar> eigIntgtor(200);

    //define a quadrature rule
    Eigen::Integrator<Scalar>::QuadratureRule quadratureRule = Eigen::Integrator<Scalar>::GaussKronrod61;

    //define the desired absolute and relative errors
    Scalar desAbsErr = Scalar(0.);
    Scalar desRelErr = Eigen::NumTraits<Scalar>::epsilon() * 50.;

    //integrate
    Scalar result = eigIntgtor.quadratureAdaptive(inFctr, Scalar(0.),Scalar(1.), desAbsErr, desRelErr, quadratureRule);

    //expected result
    Scalar expected = Scalar(-4.);

    //print output
    size_t outputPrecision  = 18;
    std::cout<<std::fixed;
    std::cout<<"result          = "<<std::setprecision(outputPrecision)<<result<<std::endl;
    std::cout<<"exact result    = "<<std::setprecision(outputPrecision)<<expected<<std::endl;
    std::cout<<"actual error    = "<<std::setprecision(outputPrecision)<<(expected-result)<<std::endl;

    return 0.0;
}

TEST_F(MyTest, test_integrator_example) {

// Test a basic integrator out
double arg = 0.0;
double res = do_something(arg);

}



/***********************************************************************************************************************
 * HOW TO FIT USING CERES-SOLVER
 **********************************************************************************************************************/


//#include <stdlib.h>
//#include <stdio.h>
//#include <unistd.h>
//#include <stdbool.h>
//#include "matio.h"
#include "ceres/ceres.h"

#include <Eigen/Core>


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

void my_fitting_func() {

    // Run the ceres solver

    // Define the variable to solve for with its initial value. It will be mutated in place by the solver.
    double x = 0.5;
    const double initial_x = x;

    // Build the problem.
    Problem problem;

    // Set up the only cost function (also known as residual). This uses
    // auto-differentiation to obtain the derivative (jacobian).
    CostFunction *cost_function = new AutoDiffCostFunction<CostFunctor, 1, 1>(new CostFunctor);
    problem.AddResidualBlock(cost_function, NULL, &x);

    // Run the solver!
    Solver::Options ceroptions;
    ceroptions.minimizer_progress_to_stdout = true;
    Solver::Summary summary;
    Solve(ceroptions, &problem, &summary);
    std::cout << summary.BriefReport() << "\n";
    std::cout << "x : " << initial_x << " -> " << x << "\n";

}


/***********************************************************************************************************************
 * HOW TO DO AN FFT USING EIGEN
 **********************************************************************************************************************/


#include <Eigen/Core>
#include <unsupported/Eigen/FFT>


void my_fft_func() {
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
}


/***********************************************************************************************************************
 * HOW TO DO AN FFT USING MKL DIRECTLY
 **********************************************************************************************************************/

#include "mkl_dfti.h"

void my_mkl_fft() {
    // Run an example FFT

    // Make meaningless data and a size vector
    float xf[200][100];
    MKL_LONG len[2] = {200, 100};

    // Create a decriptor, which is a pattern for what operation an FFT will undertake
    DFTI_DESCRIPTOR_HANDLE fft;
    DftiCreateDescriptor (&fft, DFTI_SINGLE, DFTI_REAL, 2, len);
    DftiCommitDescriptor(fft);

    // Compute a forward transform, in-place on the data
    DftiComputeForward(fft, xf);

    // Free the descriptor
    DftiFreeDescriptor(&fft);

    std::cout << "FFT COMPLETE\n";
}


/***********************************************************************************************************************
 * HOW TO PLOT REYNOLDS STRESS PROFILES, EDDY SIGNATURES AND STRENGTH / SCALE DISTRIBUTIONS WITH CPPLOT
 **********************************************************************************************************************/

#include "cpplot.h"

void my_rs_plotting_func() {

    // Plot the Reynolds Stress profiles (comes from get_spectra)
    cpplot::Figure figa = cpplot::Figure();
    for (auto i = 0; i < 6; i++) {
        cpplot::ScatterPlot p = cpplot::ScatterPlot();
        p.x = Eigen::VectorXd::LinSpaced(data.reynolds_stress_a.rows(), 1, data.reynolds_stress_a.rows());
        p.y = data.reynolds_stress_a.col(i).matrix();
        figa.add(p);
    }
    figa.write("test_t2w_rij_a.json");
    //    legend({'R13A'; 'R13B'})
    //    xlabel('\lambda_E')

    cpplot::Figure figb = cpplot::Figure();
    for (auto i = 0; i < 6; i++) {
        cpplot::ScatterPlot p = cpplot::ScatterPlot();
        p.x = Eigen::VectorXd::LinSpaced(data.reynolds_stress_b.rows(), 1, data.reynolds_stress_b.rows());
        p.y = data.reynolds_stress_b.col(i).matrix();
        figb.add(p);
    }
    figb.write("test_t2w_rij_b.json");
    //    legend({'R13A'; 'R13B'})
    //    xlabel('\lambda_E')


    // Plot the eddy signatures
    cpplot::Figure fig2 = cpplot::Figure();
    cpplot::ScatterPlot ja = cpplot::ScatterPlot();
    ja.x = Eigen::VectorXd::LinSpaced(j13a.rows(), 1, j13a.rows());
    ja.y = j13a.matrix();
    fig2.add(ja);
    cpplot::ScatterPlot jb = cpplot::ScatterPlot();
    jb.x = ja.x;
    jb.y = j13b.matrix();
    fig2.add(jb);
    fig2.write("test_t2w_j13ab.json");
    //    legend({'J13A'; 'J13B'})

    // Plot strength and scale distributions
    cpplot::Figure fig3 = cpplot::Figure();
    cpplot::ScatterPlot twa = cpplot::ScatterPlot();
    twa.x = Eigen::VectorXd::LinSpaced(minus_t2wa.rows(), 1, minus_t2wa.rows());
    twa.y = minus_t2wa.matrix();
    fig3.add(twa);
    cpplot::ScatterPlot twb = cpplot::ScatterPlot();
    twb.x = twa.x;
    twb.y = minus_t2wb.matrix();
    fig3.add(twb);
    fig3.write("test_t2w_t2wab.json");
    //    legend({'T^2\omegaA'; 'T^2\omegaB'})

}


