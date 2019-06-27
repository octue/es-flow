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


// Remove template specialisation from doc (causes duplicate) @cond

/***********************************************************************************************************************
 * HOW TO DO PIECEWISE CUBIC HERMITE INTERPOLATION
 **********************************************************************************************************************/


// TODO This is a work in progress. It gives a pchip spline, but not one which is shape preserving. Need gradient
//  conditions at the maxima, minima and endpoints.
#include <iostream>
#include <algorithm>
#include <math.h>
#include <Eigen/Core>
#include <Eigen/Dense>
double evaluate(const double xi, const Eigen::VectorXd &x, const ceres::CubicInterpolator<ceres::Grid1D<double> > &interpolator){

    const Eigen::Index k = x.rows();
    Eigen::Index n_rows = x.rows()-1;
    Eigen::VectorXd x_percent;
    Eigen::VectorXd diff = x.bottomRows(n_rows) - x.topRows(n_rows);
//    if ((diff <= 0.0).any()) {
//        throw std::invalid_argument("Input values x must be strictly increasing");
//    }
    cumsum(x_percent, diff);
    x_percent = x_percent / x_percent(n_rows-1);
    double x_max = x.maxCoeff();
    double x_min = x.minCoeff();

    // Out of range xi values are constrained to the endpoints
    double bounded_target = std::max(std::min(xi, x_max), x_min);

    // Run a binary search for where the node is in the grid. This makes the algorithm O(log(N)) to evaluate one location
// Thanks to: https://stackoverflow.com/questions/6553970/find-the-first-element-in-a-sorted-array-that-is-greater-than-the-target
    Eigen::Index low = 0, high = k; // numElems is the size of the array i.e arr.size()
    while (low != high) {
        Eigen::Index mid = (low + high) / 2; // Or a fancy way to avoid int overflow
        if (x[mid] <= bounded_target) {
            /* This index, and everything below it, must not be the first element
             * greater than what we're looking for because this element is no greater
             * than the element.
             */
            low = mid + 1;
        }
        else {
            /* This element is at least as large as the element, so anything after it can't
             * be the first element that's at least as large.
             */
            high = mid;
        }
    }
    low = low - 1;

    /* Now, low and high surround to the element in question. */

    std::cout << low << " " << high << std::endl;

    double node_low = x[low];
    double node_high = x[high];

    std::cout << "node_low " << node_low << std::endl;
    std::cout << "node_high " << node_low << std::endl;
    std::cout << "low " << low << std::endl;
    std::cout << "high " << high << std::endl;


    double proportion = (bounded_target - node_low) / (node_high - node_low);
    double xi_mapped = double(low) + proportion;

    std::cout << "proportion " << proportion << std::endl;
    std::cout << "xi_mapped " << xi_mapped << std::endl;
    double f, dfdx;

    interpolator.Evaluate(xi_mapped, &f, &dfdx);
    std::cout << "evaluated: " << f << std::endl;

    return f;
}

double evaluateDirect(const double xi, const ceres::CubicInterpolator<ceres::Grid1D<double> > &interpolator){
    // This doesn't work because it doesn't map the proportion of the way through correctly
    double f, dfdx;

    interpolator.Evaluate(xi, &f, &dfdx);
    std::cout << "evaluated: " << f << std::endl;

    return f;
}

TEST_F(InterpTest, test_pchip_interp_double) {

Eigen::VectorXd x(5);
Eigen::VectorXd y(5);
Eigen::VectorXd xi = Eigen::VectorXd::LinSpaced(100, 0, 5);
Eigen::VectorXd yi(100);
x << 2, 4.5, 8, 9, 10;
y << 6, 1, 10, 12, 19;

const Eigen::Index k = y.rows();
ceres::Grid1D<double> grid(y.data(), 0, k);
// TODO the interpolator, unlike MATLAB's routine, isn't shape preserving, nor does it adjust for non-monotonic x.
//  So this works as an interpolant, but is pretty horrid.
ceres::CubicInterpolator<ceres::Grid1D<double> > interpolator(grid);

double x_max = x.maxCoeff();
double x_min = x.minCoeff();
for (Eigen::Index i=0; i<xi.rows(); i++) {
//yi[i] = evaluateDirect(xi[i], interpolator);
//xi[i] = xi[i] * ((x_max - x_min) / k) + x_min;
yi[i] = evaluate(xi[i], x, interpolator);
}


Figure fig = Figure();
ScatterPlot p_xy = ScatterPlot();
p_xy.x = x;
p_xy.y = y;
p_xy.name = "xy";
fig.add(p_xy);
ScatterPlot p_xiyi = ScatterPlot();
p_xiyi.x = xi;
p_xiyi.y = yi;
p_xiyi.name = "xiyi";
fig.add(p_xiyi);

// Add axis labeling
Layout lay = Layout("pchip check");
lay.xTitle("$x$");
lay.yTitle("$y$");
fig.setLayout(lay);

// Write figures
fig.write("check_pchip.json");

}

// TODO Once done, you can wrap this up into an interpolant much like in the interp.h file:



#include "ceres/ceres.h"
#include "ceres/cubic_interpolation.h"
#include "glog/logging.h"
#include "utilities/cumsum.h"

/** Helper to build an interpolator, avoiding overly complex initialisation list
 *
 * @param x_vec
 * @param y_vec
 * @return
 */
ceres::CubicInterpolator<ceres::Grid1D<double> > getInterpolator(Eigen::VectorXd const &x_vec, Eigen::VectorXd const &y_vec) {
    const Eigen::Index k = y_vec.rows();
    // TODO file an issue over at ceres-solver. Grid1D is deeply unhelpful. At the very least it should accept Eigen::Index instead of int.
    std::cout << "k: " << k << std::endl;
    ceres::Grid1D<double> grid(y_vec.data(), 0, k);
    ceres::CubicInterpolator<ceres::Grid1D<double> > interpolator(grid);
    return interpolator;

}

/** @brief Performs cubic interpolation using a piecewise cubic hermite interpolant.
 *
 * Written in response to:
 *    https://stackoverflow.com/questions/56789920/how-do-i-use-cerescubicinterpolator-with-data-not-on-a-uniform-grid
 *
 * Based on:
 *    https://github.com/ceres-solver/ceres-solver/blob/master/examples/sampled_function/sampled_function.cc
 *
 * Also consider:
 *    https://blog.demofox.org/2015/08/08/cubic-hermite-interpolation/
 *
 * This is much friendlier than CubicSplineInterpolant, whose basis function is a bezier spline, since there's no
 * wildly oscillatory behaviour where the input nodes are jagged or noisy.
 *
 * Initialise and use an interpolant as follows:
 *
 *    Eigen::VectorXd xvals(3);
 *    Eigen::VectorXd yvals(xvals.rows());
 *
 *    xvals << 0, 15, 30;
 *    yvals << 0, 12, 17;
 *
 *    PiecewiseCubicHermiteInterpolant s(xvals, yvals);
 *
 */
class PiecewiseCubicHermiteInterpolant {
public:
    // TODO template for arrays
    PiecewiseCubicHermiteInterpolant(Eigen::VectorXd const &x_vec, Eigen::VectorXd const &y_vec)
        : k_(y_vec.rows()),
          interpolator_(getInterpolator(x_vec, y_vec)) {
        // TODO assert that x is strictly ascending or strictly descending
        std::cout << "HEREbf" <<std::endl;
        Eigen::Index n_rows = x_vec.rows();
        Eigen::VectorXd x_percent;
        cumsum(x_percent, x_vec.bottomRows(n_rows-1) - x_vec.topRows(n_rows-1));
    }


    /** @brief Evaluate interpolant at values xi whose yi values are unknown
     *
     * Performs piecewise cubic interpolation using the ceres CubicInterpolator
     *
     *@param xi double (or eigen VectorXd) of target x values for the interpolation
     *@return double interpolated value yi corresponding to location xi
     */
    double operator()(double xi) const {
        // x values need to be scaled down in extraction as well.
        std::cout << "evaluating" << std::endl;
        //
        double f, dfdx;
        interpolator_.Evaluate(xi, &f, &dfdx);
        std::cout << "evaluated" << std::endl;

        return f;
    }

    Eigen::VectorXd operator()(Eigen::VectorXd const &xi_vec) {
        Eigen::VectorXd yi_vec(xi_vec.rows());
        yi_vec = xi_vec.unaryExpr([this](double xi) { return this->operator()(xi); }).transpose();
        return yi_vec;
    }

private:


    /** @brief Scales x values to the domain [0, 1] on which the interpolant was defined
     *
     *@param x double (or eigen VectorXd)
     *@return double (or eigen VectorXd) the scaled value(s)
     */
//    double scaled_value(double x) const {
//        return (x - x_min) / (x_max - x_min);
//    }

//    Eigen::RowVectorXd scaled_values(Eigen::VectorXd const &x_vec) const {
//        return x_vec.unaryExpr([this](double x) { return scaled_value(x); }).transpose();
//    }

    // Number of samples in the input data
    Eigen::Index k_;

    // Interpolator, initialised with the grid of values
    const ceres::CubicInterpolator<ceres::Grid1D<double> > interpolator_;
};







/***********************************************************************************************************************
 * HOW TO DIFFERENTIATE USING EIGEN::AUTODIFF
 **********************************************************************************************************************/


#include <unsupported/Eigen/AutoDiff>
/**
 * Application of automatic differentiation with one variable, using a templated function which can be called
 * with active scalars or normal floats/doubles.
 */
template <typename T>
T myfun2(const double a, T const & b){
    T c = pow(sin(a),2.) + pow(cos(b),2.) + 1.;
    return c;
}

void test_scalar() {
    std::cout << "== test_scalar() ==" << std::endl;
    double a;
    a = 0.3;
    typedef Eigen::AutoDiffScalar<Eigen::VectorXd> AScalar;
    AScalar Ab;
    Ab.value() = 0.5;
    Ab.derivatives() = Eigen::VectorXd::Unit(1,0);
    AScalar Ac = myfun2(a,Ab);
    std::cout << "Result: " << Ac.value() << std::endl;
    std::cout << "Gradient: " << Ac.derivatives().transpose() << std::endl;
}


/**
 * Typical application of automatic differentiation in two dimensions
 * A templated function which can be called with active scalars
 * or normal floats/doubles.
 */
template <typename T>
T myfun(T const & a, T const & b){
	T c = pow(sin(a),2.) + pow(cos(b),2.) + 1.;
	return c;
}

void test_scalar(){
	std::cout << "== test_scalar() ==" << std::endl;
	// use with normal floats
	double a,b;
	a = 0.3;
	b = 0.5;
	double c = myfun(a,b);
	std::cout << "Result: " << c << std::endl;

	// use with AutoDiffScalar
	typedef Eigen::AutoDiffScalar<Eigen::VectorXd> AScalar;
	AScalar Aa,Ab;
	Aa.value() = 0.3;
	Ab.value() = 0.5;
	Aa.derivatives() = Eigen::VectorXd::Unit(2,0); // This is a unit vector [1, 0]
	Ab.derivatives() = Eigen::VectorXd::Unit(2,1); // This is a unit vector [0, 1]
	AScalar Ac = myfun(Aa,Ab);
	std::cout << "Result: " << Ac.value() << std::endl;
	std::cout << "Gradient: " << Ac.derivatives().transpose() << std::endl;
}


/***********************************************************************************************************************
 * HOW TO INTEGRATE USING THE NUMERICALINTEGRATION LIBRARY
 **********************************************************************************************************************/

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

/***********************************************************************************************************************
 * HOW TO DO MAT FILE READ WRITE
 **********************************************************************************************************************/

#include "matio.h"

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




// @endcond
