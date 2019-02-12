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



/** @brief Cost functor for fitting lewkowicz speed profiles.
 *
 * Implements operator() as required by ceres-solver. Allows some parameters to be fixed (see ``fit_lewkowicz_speed``).
 *
 */
struct LewkowiczSpeedResidual {
    /** Construct the cost functor.
     *
     * Parameters, where supplied in input ``initial_params`` and masked in input ``fix_param`` are in this order:
     * pi_coles = params[0]
     * kappa = params[1]
     * u_inf = params[2]
     * u_tau = params[3]
     * delta_c = params[4]
     *
     * @param z Observation data point, vertical coordinate in m
     * @param u Observation data point, horizontal speed in m/s
     * @param fix_param Eigen::Array<bool, 5, 1>, logical mask, true where the parameter is fixed, false where it is allowed to vary
     * @param initial_params Eigen::Array<double, 5, 1> containing initial parameter values. These are used where fixed parameters are specified
     */
    LewkowiczSpeedResidual(double z, double u, Array5b &fix_param, Array5d &initial_params) : z_(z), u_(u), fix_param_(fix_param), initial_params_(initial_params) { }

    template <typename T> bool operator()(T const* const* parameters, T* residual) const {

        // Mask out fixed and variable parameters
        // TODO Once Eigen 3.4.x comes out, reduce this code to use the parameter mask as a logical index.
        //  See this issue which will allow such indexing in eigen: http://eigen.tuxfamily.org/bz/show_bug.cgi?id=329#c27
        int param_ctr = 0;
        std::vector<T> params(5);
        for (auto i = 0; i < 5; i++) {
            if (fix_param_(i)) {
                params[i] = T(initial_params_(i));
            } else {
                params[i] = parameters[param_ctr][0];
                param_ctr += 1;
            }
        }
        T spd = lewkowicz_speed(T(z_), params[0], params[1], params[2], params[3], params[4]); // TODO can I spread these?
        residual[0] = u_ - spd;
        return true;
    }
private:
    const double u_;
    const double z_;
    Array5b fix_param_;
    Array5d initial_params_;
};


/** @brief Determine best fit of Lewkowicz relation to input speed profile.
 *
 * TODO Expose nondefault initial values and parameter fixing options on the API to this function.
 *
 * Present limitations:
 *      - von Karman constant Kappa is fixed according to the value in definitions.h.
 *      - Boundary layer thickness delta_c is fixed to a default value of 1000 m.
 *
 * @param[in] z Vertical locations in m (or normalised, if z_ref = 1.0)
 * @param[in] u Observed speeds at corresponding vertical locations, in m/s
 * @param[in] z_ref Optional reference height (default 1.0) by which the input z is normalised
 * @param[in] u_ref Optional reference velocity (default 1.0) by which the input u is normalised
 * @return
 */
Array5d fit_lewkowicz_speed(const Eigen::ArrayXd &z, const Eigen::ArrayXd &u) {

    // TODO Assert that z.size() == u.size()

    double pi_coles = 0.75;
    double kappa = KAPPA_VON_KARMAN;
    double shear_ratio = 10;
    double u_inf = u.maxCoeff(); // TODO make this less sensitive to noise. Maybe median?
    double u_tau = u_inf / shear_ratio;
    double delta_c = 1000;

    Array5b fix_params;
    Array5d initial_params;
    fix_params << false, true, false, false, true;
    initial_params << pi_coles, kappa, u_inf, u_tau, delta_c;

    // Initialise a results array, and get a vector of pointers to the free parameters
    Array5d final_params(initial_params);
    std::vector<double *> unfixed_param_ptrs;
    for (int p_ctr = 0; p_ctr <5; p_ctr++) {
        if (!fix_params(p_ctr)) { unfixed_param_ptrs.push_back(final_params.data()+p_ctr); };
    }
    std::cout << "HERE2 " << unfixed_param_ptrs[0][0] << " "  << unfixed_param_ptrs[1][0] << " " << unfixed_param_ptrs[2][0] << " " <<std::endl;
    // Build the problem
    Problem problem;
    Solver::Summary summary;
    Solver::Options options;
    options.max_num_iterations = 25;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = true;
    
    // Set up the only cost function (also known as residual), using cauchy loss function for
    // robust fitting and auto-differentiation to obtain the jacobian
    for (auto i = 0; i < z.size(); i++) {

        // Set the stride such that the entire set of derivatives is computed at once (since we have maximum 5)
        DynamicAutoDiffCostFunction<LewkowiczSpeedResidual, 5>* cost_function =
            new DynamicAutoDiffCostFunction<LewkowiczSpeedResidual, 5>(
                new LewkowiczSpeedResidual(z[i], u[i], fix_params, initial_params)
            );
        // Add N parameters, where N is the number of free parameters in the problem
        auto pc = 0;
        for (auto p_ctr = 0; p_ctr < 5; p_ctr++) {
            if (!fix_params(p_ctr)) {
                cost_function->AddParameterBlock(1);
                pc += 1;
            }
        }
        std::cout << "added n params n = " << pc << std::endl;
        cost_function->SetNumResiduals(1);
        problem.AddResidualBlock(cost_function, new CauchyLoss(0.5), &pi_coles, &u_inf, &u_tau);
//        problem.AddResidualBlock(cost_function, new CauchyLoss(0.5), unfixed_param_ptrs);
    }

    // Run the solver
    Solve(options, &problem, &summary);
    std::cout << summary.BriefReport() << std::endl;
    std::cout << "Initial values: " << initial_params.transpose() << std::endl;
    final_params(0) = pi_coles;
    final_params(2) = u_inf;
    final_params(3) = u_tau;
    std::cout << "Final   values: " << final_params.transpose() << std::endl;

    return final_params;
}

} /* namespace es */


#endif /* SOURCE_FIT_H_ */
