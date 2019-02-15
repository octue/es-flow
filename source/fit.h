/*
 * Fit.h Use ceres-solver to fit mean speed and spectra to measured data
 *
 * Author:              Tom Clark  (thclark @ github)
 *
 * Copyright (c) 2019 Octue Ltd. All Rights Reserved.
 *
 */

#ifndef ES_FLOW_FIT_H
#define ES_FLOW_FIT_H

#include <stdio.h>
#include "ceres/ceres.h"
#include <Eigen/Core>

#include "relations/velocity.h"
#include "definitions.h"

using ceres::AutoDiffCostFunction;
using ceres::DynamicAutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;
using ceres::CauchyLoss;


typedef Eigen::Array<bool, 5, 1> Array5b;
typedef Eigen::Array<double, 5, 1> Array5d;


namespace es {


/** @brief Cost functor for fitting power law speed profiles.
 *
 * Templated for automatic differentiation.
 *
 * Internal use only (see ``fit_power_law_speed``).
 *
 */
struct PowerLawSpeedResidual {
	PowerLawSpeedResidual(double z, double u, double z_ref=1.0, double u_ref=1.0) : z_(z), u_(u), u_ref_(u_ref), z_ref_(z_ref) {}
	template <typename T> bool operator()(const T* const alpha, T* residual) const {

		residual[0] = u_ - power_law_speed(T(z_), u_ref_, z_ref_, alpha[0]);
		return true;
	}
private:
    // Observations for a sample
    const double u_;
    const double z_;
    const double u_ref_;
    const double z_ref_;
};


/** @brief Determine best fit of power law to input speed profile.
 *
 * @param[in] z Vertical locations in m (or normalised, if z_ref = 1.0)
 * @param[in] u
 * @param[in] z_ref Optional reference height (default 1.0) by which the input z is normalised
 * @param[in] u_ref Optional reference velocity (default 1.0) by which the input u is normalised
 * @return
 */
double fit_power_law_speed(const Eigen::ArrayXd &z, const Eigen::ArrayXd &u, const double z_ref=1.0, const double u_ref=1.0) {

	// Define the variable to solve for with its initial value. It will be mutated in place by the solver.
	double alpha = 0.3;

	// TODO Assert that z.size() == u.size()

	// Build the problem
	Problem problem;
    Solver::Summary summary;
    Solver::Options options;

    options.max_num_iterations = 400;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = true;
//std::cout << options << std::endl;
	// Set up the only cost function (also known as residual), using cauchy loss function for
	// robust fitting and auto-differentiation to obtain the derivative (jacobian)
	for (auto i = 0; i < z.size(); i++) {
	    CostFunction *cost_function = new AutoDiffCostFunction<PowerLawSpeedResidual, 1, 1>(
	        new PowerLawSpeedResidual(z[i], u[i], z_ref, u_ref));
        problem.AddResidualBlock(cost_function, new CauchyLoss(0.5), &alpha);
	}

	// Run the solver
    Solve(options, &problem, &summary);
    std::cout << summary.BriefReport() << "\n";
    std::cout << "Initial alpha: " << 0.3 << "\n";
    std::cout << "Final   alpha: " << alpha << "\n";

    return alpha;
}


/** @brief Cost functor for fitting lewkowicz speed profiles.
 *
 * Implements operator() as required by ceres-solver. Allows some parameters to be fixed.
 *
 * Internal use only (see ``fit_lewkowicz_speed``).
 */
struct LewkowiczSpeedResidual {

    /** Construct the cost functor.
     *
     * Parameters, where supplied in input ``initial_params`` and masked in input ``fixed_params`` are in this order:
     *    - ``pi_coles = params[0]``
     *    - ``kappa = params[1]``
     *    - ``u_inf = params[2]``
     *    - ``shear_ratio = params[3]``
     *    - ``delta_c = params[4]``
     *
     * @param[in] z Observation data point, vertical coordinate in m
     * @param[in] u Observation data point, horizontal speed in m/s
     * @param[in] fixed_params Eigen::Array<bool, 5, 1>, logical mask, true where the parameter is fixed, false where it is allowed to vary
     * @param[in] initial_params Eigen::Array<double, 5, 1> containing initial parameter values. These are used where fixed parameters are specified
     */
    LewkowiczSpeedResidual(double z, double u, const Array5b &fixed_params, const Array5d &initial_params) : z_(z), u_(u), fixed_params_(fixed_params), initial_params_(initial_params) { }

    template <typename T> bool operator()(T const* const* parameters, T* residual) const {

        // Mask out fixed and variable parameters
        // TODO Once Eigen 3.4.x comes out, reduce this code to use the parameter mask as a logical index.
        //  See this issue which will allow such indexing in eigen: http://eigen.tuxfamily.org/bz/show_bug.cgi?id=329#c27
        int param_ctr = 0;
        std::vector<T> params(5);
        for (auto i = 0; i < 5; i++) {
            if (fixed_params_(i)) {
                params[i] = T(initial_params_(i));
            } else {
                params[i] = parameters[param_ctr][0];
                param_ctr += 1;
            }
        }
        T spd = lewkowicz_speed(T(z_), params[0], params[1], params[2], params[3], params[4]);
        residual[0] = u_ - spd;
        return true;
    }
private:
    const double u_;
    const double z_;
    Array5b fixed_params_;
    Array5d initial_params_;
};


/** @brief Determine best fit of Lewkowicz relation to input speed profile.
 *
 * Robustly fits Lewkowicz speed relation @f$ u = f(z, \Pi, \kappa, U_{inf}, S, \delta_{c}) @f$ to input
 * data @f$ u, z @f$.
 *
 * An initial guess is required, and parameters can be fixed to help constrain the fitting process. See also
 * fit_lewkowicz_speed(const Eigen::ArrayXd &z, const Eigen::ArrayXd &u, bool print_report = true) for providing a set
 * of default values and fixed parameters.
 *
 * Parameters are as follows:
 *    - ``pi_coles = params[0]``
 *    - ``kappa = params[1]``
 *    - ``u_inf = params[2]``
 *    - ``shear_ratio = params[3]``
 *    - ``delta_c = params[4]``
 *
 * The fitting process uses a Levenberg-Marquadt solver (provided by ``ceres-solver``) with Automatic differentiation
 * (for improved numeric stability and convergence over numerical differentiation approaches), and a Cauchy Loss
 * function of parameter 0.5 (for robustness agains outlying data entries).
 *
 * @param[in] z Vertical locations in m (or normalised, if z_ref = 1.0)
 * @param[in] u Observed speeds at corresponding vertical locations, in m/s
 * @param[in] fixed_params Eigen::Array<bool, 5, 1>, logical mask, true where the parameter is fixed, false where it is allowed to vary
 * @param[in] initial_params Eigen::Array<double, 5, 1> containing initial parameter values. These are used where fixed parameters are specified
 * @param[in] print_report bool, default true. If true, the solver prints iterations summary and full final convergence and solution report (for debugging and validation purposes).
 *
 * @return Array5d containing fitted values ``[pi_coles, kappa, u_inf, shear_ratio, delta_c]``
 */
Array5d fit_lewkowicz_speed(const Eigen::ArrayXd &z, const Eigen::ArrayXd &u, const Array5b &fixed_params, const Array5d &initial_params, bool print_report = true) {

    // Initialise a results array, and get a vector of pointers to the 'unfixed' parameters (the ones to optimise)
    Array5d final_params(initial_params);
    std::vector<double *> unfixed_param_ptrs;
    for (int p_ctr = 0; p_ctr <5; p_ctr++) {
        if (!fixed_params(p_ctr)) { unfixed_param_ptrs.push_back(final_params.data()+p_ctr); };
    }

    // Build the problem
    Problem problem;
    Solver::Summary summary;
    Solver::Options options;
    options.max_num_iterations = 400;
    options.minimizer_progress_to_stdout = print_report;

    for (auto i = 0; i < z.size(); i++) {

        // Set the stride such that the entire set of derivatives is computed at once (since we have maximum 5)
        auto cost_function = new DynamicAutoDiffCostFunction<LewkowiczSpeedResidual, 5>(
            new LewkowiczSpeedResidual(z[i], u[i], fixed_params, initial_params)
        );

        // Add N parameters, where N is the number of free parameters in the problem
        auto pc = 0;
        for (auto p_ctr = 0; p_ctr < 5; p_ctr++) {
            if (!fixed_params(p_ctr)) {
                cost_function->AddParameterBlock(1);
                pc += 1;
            }
        }
        cost_function->SetNumResiduals(1);
        problem.AddResidualBlock(cost_function, new CauchyLoss(0.5), unfixed_param_ptrs);
    }

    Solve(options, &problem, &summary);

    if (print_report) {
        std::cout << summary.FullReport() << std::endl;
        std::cout << "Initial values: " << initial_params.transpose() << std::endl;
        std::cout << "Final   values: " << final_params.transpose() << std::endl;
    }

    return final_params;
}


/** @brief Determine best fit of Lewkowicz relation to input speed profile.
 *
 * Robustly fits Lewkowicz speed relation @f$ u = f(z, \Pi, \kappa, U_{inf}, S, \delta_{c}) @f$ to input
 * data @f$ u, z @f$, using 'sensible' initial guesses (for Atmospheric Boundary Layer in Northern European conditions)
 * and fixed parameters to constrain the fitting process:
 *
 *     - von Karman constant @f$ \kappa @f$ is fixed according to the value in definitions.h.
 *     - Boundary layer thickness @f$ delta_{c} @f$ is fixed to a default value of ``1000 m``.
 *
 * @param[in] z Vertical locations in m (or normalised, if z_ref = 1.0)
 * @param[in] u Observed speeds at corresponding vertical locations, in m/s
 * @param[in] print_report bool, default true. If true, the solver prints iterations summary and full final convergence and solution report (for debugging and validation purposes).
 *
 * @return Array5d containing fitted values ``[pi_coles, kappa, u_inf, shear_ratio, delta_c]``
 */
Array5d fit_lewkowicz_speed(const Eigen::ArrayXd &z, const Eigen::ArrayXd &u, bool print_report = true) {

    // Set default initial parameters
    double pi_coles = 0.5;
    double kappa = KAPPA_VON_KARMAN;
    double shear_ratio = 20;
    double u_inf = u.maxCoeff();
    double delta_c = 1000;
    Array5d initial_params;
    initial_params << pi_coles, kappa, u_inf, shear_ratio, delta_c;

    // Set default fixed parameters (von karman constant and boundary layer thickness)
    Array5b fixed_params;
    fixed_params << false, true, false, false, true;

    return fit_lewkowicz_speed(z, u, fixed_params, initial_params, print_report);

}

} /* namespace es */


#endif // ES_FLOW_FIT_H
