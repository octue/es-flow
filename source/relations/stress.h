/*
 * stress.h Atmospheric Boundary Layer Reynolds Stress relations
 *
 * Author:              Tom Clark  (thclark @ github)
 *
 * Copyright (c) 2019 Octue Ltd. All Rights Reserved.
 *
 */

#ifndef ES_FLOW_STRESS_H_
#define ES_FLOW_STRESS_H_

#include <Eigen/Dense>
#include <Eigen/Core>
#include <iostream>
#include <iomanip>
#include <math.h>
#include "NumericalIntegration.h"
#include "profile.h"
#include "utilities/trapz.h"
#include "utilities/cumtrapz.h"

#include "cpplot.h"

using namespace utilities;
using namespace cpplot;


namespace es {


/** @brief Get the Coles wake distribution @f$ w_{c} @f$ from the wake parameter @f$ \Pi @f$.
 *
 */
template<typename T>
T get_coles_wake(const T eta, const double pi_coles) {
    T wc;
    T eta_pow_2 = pow(eta, 2.0);
    wc = (2.0 * eta_pow_2 * (3.0 - (2.0 * eta)));
//        - (1.0 / pi_coles) * eta_pow_2 * (1.0 - eta) * (1.0 - (2.0 * eta));
    std::cout << "REMOVE THE mult BODGE" << std::endl;
    return wc;
}

// Remove template specialisation from doc (causes duplicate) @cond
Eigen::ArrayXd get_coles_wake(const Eigen::ArrayXd &eta, const double pi_coles) {
    Eigen::ArrayXd wc;
    Eigen::ArrayXd eta_pow_2;
    eta_pow_2 = eta.pow(2.0);
    wc = (2.0 * eta_pow_2 * (3.0 - (2.0 * eta)));
//        - (1.0 / pi_coles) * eta_pow_2 * (1.0 - eta) * (1.0 - (2.0 * eta));

    std::cout << "REMOVE THE mult BODGE" << std::endl;

    return wc;
}
// @endcond


/** @brief Get the integrand @f$ f @f$ used in computation of @f$ R_{13} @f$ in the modified Lewkowicz method.
 *
 * Uses equations 2 and 7 Perry and Marusic 1995 Part 1.
 *
 */
Eigen::ArrayXd get_f_integrand(const Eigen::ArrayXd &eta, const double kappa, const double pi_coles, const double shear_ratio) {

    Eigen::ArrayXd f;
    Eigen::ArrayXd ones;
    ones.setOnes(eta.size());
    f = (-1.0 / kappa) * eta.log()
        + (pi_coles/kappa) * get_coles_wake(ones, pi_coles)
        - (pi_coles/kappa) * get_coles_wake(eta, pi_coles);
    for (int k = 0; k < f.size(); k++) {
        if (std::isinf(f[k])) {
            f(k) = shear_ratio; // from eqs 2 and 7 P&M part 1 1995
        }
    }
    return f;
}
Eigen::ArrayXd get_f_integrand_vik(const Eigen::ArrayXd &eta, const double kappa, const double pi_coles, const double shear_ratio) {

    Eigen::ArrayXd f;
    Eigen::ArrayXd ones;
    ones.setOnes(eta.size());
    f = (eta.pow(3.0) - 1.0)/(3.0*kappa)
        -(1.0 / kappa) * eta.log()
        + 2.0*pi_coles*(1.0 - 3.0*eta.pow(2.) + 2.0*eta.pow(3.0))/kappa;
    std::cout << "USING VIK" <<std::endl;
    for (int k = 0; k < f.size(); k++) {
        if (std::isinf(f[k])) {
            f(k) = shear_ratio; // from eqs 2 and 7 P&M part 1 1995
        }
    }
    return f;
}



/** @brief Compute Horizontal-vertical Reynolds Stress R13 profile using modified Lewkowicz formulation.
 *
 * Gets Reynolds Stress profiles due to Type A and B eddies. Adopts the approach of Perry and Marusic 1995,
 * using the modified Lewkowicz formulation (Perry and Marusic eq.51).
 *
 * TODO - presently, the integrations here are sensitive to the chosen eta distribution. Refactor according to
 * issue #38
 *
 * TODO - Determine whether it will be an advantage to template this function so that it can be called with active scalars (allows use of autodiff), doubles/floats,
 * Eigen::Arrays (directly) or Eigen::VectorXds (via template specialisation) of eta values.
 *
 * # References
 *  [1] Perry AE and Marusic I (1995) A wall-wake model for turbulent boundary
 *      layers. Part 1. Extension of the attached eddy hypothesis J Fluid Mech
 *      vol 298 pp 361-388
 *
 * # Future Improvements
 *
 *   [1] Optional different mean profile formulations including account for the
 *       free surface
 *
 *   [2] Support directional variation with height; i.e. compatible with mean
 *       profile formulations using U(y)
 *
 *   [3] Added formulation for contribution of smaller Type C eddies
 *
 *   [4] For the given wall formulation, can f be calculated more efficiently?
 *       See schlichting and gersten p. 593
 *
 * @param[in] beta Clauser parameter @f$ \Beta @f$, representing acceleration/decelaration of the boundary layer
 * @param[in] eta Nondimensional vertical coordinates at which you want to get stress @f$ R_{13} @f$. Values must ascend but not necessarily be monotonic.
 * @param[in] kappa von Karman constant.
 * @param[in] pi_coles Coles wake parameter @f$ \Pi @f$
 * @param[in] shear_ratio Ratio between free-stream and skin friction velocities @f$ S = U_{inf}/U_{\tau} @f$
 * @param[in] zeta Scaled streamwise derivative @f$ \zeta @f$ of the Coles wake parameter @f$ \Pi @f$
 *
 */
void reynolds_stress_13(Eigen::ArrayXd &r13_a, Eigen::ArrayXd &r13_b, const double beta, const Eigen::ArrayXd &eta, const double kappa, const double pi_coles, const double shear_ratio, const double zeta){

    // Check the input has a valid range; the integration will be messed up otherwise.
//    if ((eta(eta.rows()-1) != 1.0) || (eta(0) != 0.0)) {
//        throw std::invalid_argument("Input eta must be defined in the range 0, 1 exactly");
//    }

    // Get f between eta = 0 and eta = 1 (bounds for C1 integration)
    Eigen::ArrayXd f = get_f_integrand_vik(eta, kappa, pi_coles, shear_ratio);
    const double d_pi = 0.01 * pi_coles;
    Eigen::ArrayXd f_plus  = get_f_integrand_vik(eta, kappa, (pi_coles + d_pi), shear_ratio);
    Eigen::ArrayXd f_minus  = get_f_integrand_vik(eta, kappa, (pi_coles - d_pi), shear_ratio);

    // TODO can we cast this returned value directly?
    Eigen::ArrayXd c1_tmp;
    trapz(c1_tmp, eta, f);
    double c1 = c1_tmp(0);

    // Do the integrations for ei with central differencing for numerical differentiation of e coefficients 5-7
    Eigen::ArrayXd e1;
    Eigen::ArrayXd e1_plus;
    Eigen::ArrayXd e1_minus;
    Eigen::ArrayXd e2;
    Eigen::ArrayXd e2_plus;
    Eigen::ArrayXd e2_minus;
    Eigen::ArrayXd e3;
    Eigen::ArrayXd e4;
    Eigen::ArrayXd e5;
    Eigen::ArrayXd e6;
    Eigen::ArrayXd e7;
    cumtrapz(e1, eta, f);
    cumtrapz(e1_plus, eta, f_plus);
    cumtrapz(e1_minus, eta, f_minus);
    cumtrapz(e2, eta, f.pow(2.0));
    cumtrapz(e2_plus, eta, f_plus.pow(2.0));
    cumtrapz(e2_minus, eta, f_minus.pow(2.0));
    e3 = f * e1;
    e4 = eta * f;
    e5 = (e1_plus - e1_minus) / (2.0 * d_pi);
    e6 = (e2_plus - e2_minus) / (2.0 * d_pi);
    e7 = f * e5;

    // Get A coefficients from equations A2 a-d
    //MARK correct
    Eigen::ArrayXd a1 = e2 - e3 + shear_ratio * (e4 - e1);
    //MARK correct
    Eigen::ArrayXd a2 = (-2.0 * e2) + e3 + (shear_ratio * e1);
    Eigen::ArrayXd a3 = e6 - e7 - (shear_ratio * e5);
    //MARK correct
    Eigen::ArrayXd a4 = (2.0 * e2) - e3 + (shear_ratio * e4) - (shear_ratio * 3.0 * e1);

    // B coefficients are simply evaluated at eta = 1 (eqn A6)
    //MARK correct
    double b1 = a1(eta.rows()-1);
    //MARK correct
    double b2 = a2(eta.rows()-1);
    double b3 = a3(eta.rows()-1);
    //MARK correct
    double b4 = a4(eta.rows()-1);

    // E1 from (eqn A4). Can't call it E1 due to name conflict with above.
    //MARK correct
    double e1_coeff = 1.0 / (kappa * shear_ratio + 1.0);

    // N from (eqn A5) using central differencing as before
    double wc_minus = get_coles_wake(1.0, pi_coles - d_pi);
    double wc_plus =  get_coles_wake(1.0, pi_coles + d_pi);
    double n = get_coles_wake(1.0, pi_coles) + pi_coles * (wc_plus - wc_minus) / (2.0 * d_pi);
    std::cout << "BODGE: MULTIPLYING N BY 1.17 SEEMS TO HELP FOR NO REASON (getting R13 analytic in non-equilibrium cases)"<< pi_coles << std::endl;
//    n *=1.17;

    // Compile f_i terms
    Eigen::ArrayXd f1;
    Eigen::ArrayXd f2;
    Eigen::ArrayXd f3;
    f1 = 1.0
        - a1 / (b1 + e1_coeff * b2)
        - e1_coeff * a2 / (b1 + e1_coeff * b2);

    f2 = (e1_coeff * n * a2 * b1
        + a3 * b1
        - e1_coeff * n * a1 * b2
        + e1_coeff * a3 * b2
        - a1 * b3
        - e1_coeff * a2 * b3) / (b1 + e1_coeff * b2);

    f3 = (e1_coeff * a2 * b1
        + a4 * b1
        - e1_coeff * a1 * b2
        + e1_coeff * a4 * b2
        - a1 * b4
        - e1_coeff * a2 * b4) / (b1 + e1_coeff * b2);

    // Convert f2 and f3 into g1 and g2, apply eq. 15
    Eigen::ArrayXd g1 = f2 / shear_ratio;
    Eigen::ArrayXd g2 = -f3 / (c1 * shear_ratio);
    Eigen::ArrayXd minus_r13 = f1 + g1 * zeta + g2 * beta;

    // As a validation check, reproduce figure
    Figure fig = Figure();
    ScatterPlot p1 = ScatterPlot();
    p1.x = eta;
    p1.y = f1;
    p1.name = "f1";
    ScatterPlot p2 = ScatterPlot();
    p2.x = eta;
    p2.y = g1*zeta;
    p2.name = "g1*zeta";
    ScatterPlot p3 = ScatterPlot();
    p3.x = eta;
    p3.y = g2*beta;
    p3.name = "g2*beta";
    ScatterPlot p4 = ScatterPlot();
    p4.x = eta;
    p4.y = minus_r13;
    p4.name = "tau / tau_0 (= -r13)";
    fig.add(p1);
    fig.add(p2);
    fig.add(p3);
    fig.add(p4);
    Layout lay = Layout("Check reproduction of Perry & Marusic 1995 Figure 1");
    lay.xTitle("$\\eta$");
    fig.setLayout(lay);
    fig.write("check_perry_marusic_fig_1");

    std::cout << "BODGE HERE!!!! zeta, beta" << zeta << ", " << beta << std::endl;

    // Get the component due to equilibrium sink flow (OLD version - see P&M eqn 53)
    // minusReynStressA = ones(size(eta)) - eta + eta.*log(eta);

    // Lewkowicz 1982 shear stress for equilibrium sink flow, Perry and Marusic eqn. 51
    Eigen::ArrayXd minus_r13_a = 1.0 - (60.0/59.0)*eta - (20.0/59.0)*eta.pow(3.0) + (45.0/59.0)*eta.pow(4.0) - (24.0/59.0)*eta.pow(5.0) + (60.0/59.0)*eta*eta.log();

    // Handle the log(0) singularity
    for (int i = 0; i < minus_r13_a.size(); i++) {
        if (std::isinf(minus_r13_a(i)) || std::isnan(minus_r13_a(i))) {
            minus_r13_a(i) = 1.0;
        }
    }

    // Output correctly signed reynolds stress components
    r13_a = -1.0 * minus_r13_a;
    r13_b = -1.0 * (minus_r13 - minus_r13_a);

}

} /* namespace es */

#endif /* ES_FLOW_STRESS_H_ */
