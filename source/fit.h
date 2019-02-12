/*
 * Fit.h Use ceres-solver to fit mean speed and spectra to measured data
 *
 * Author:              Tom Clark  (thclark @ github)
 *
 * Copyright (c) 2019 Octue Ltd. All Rights Reserved.
 *
 */

#ifndef SOURCE_FIT_H_
#define SOURCE_FIT_H_

#include <stdio.h>
#include "ceres/ceres.h"
#include <Eigen/Core>

#include "relations/velocity.h"

using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;
using ceres::CauchyLoss;

namespace es {


/** @brief Cost functor for fitting power law speed profiles
 *
 * Suitable for automatic differentiation.
 *
 */
struct PowerLawSpeedResidual {
	PowerLawSpeedResidual(double z, double u, double z_ref=1.0, double u_ref=1.0) : z_(z), u_(u), u_ref_(u_ref), z_ref_(z_ref) {}
	template <typename T> bool operator()(const T* const alpha, T* residual) const {
//		T z_norm = z_ / z_ref_;
		residual[0] = u_ - pow((z_ / z_ref_), alpha[0]) * u_ref_;
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
    options.max_num_iterations = 25;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = true;

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


} /* namespace es */


#endif /* SOURCE_FIT_H_ */
