/*
 * velocity.h Atmospheric Boundary Layer velocity profile relations
 *
 * Author:              Tom Clark  (thclark @ github)
 *
 * Copyright (c) 2013-9 Octue Ltd. All Rights Reserved.
 *
 */

#ifndef SOURCE_RELATIONS_VELOCITY_H_
#define SOURCE_RELATIONS_VELOCITY_H_

#include <Eigen/Dense>
#include <Eigen/Core>
#include "profile.h"


namespace es {


/** @brief Compute speed profile according to the power law.
 *
 * Templated so that it can be called with active scalars (allows use of autodiff), doubles/floats,
 * Eigen::Arrays (directly) or Eigen::VectorXds (via template specialisation) of z values.
 *
 * Power law speed is computed as:
 * \f[ \frac{\overline{U}}{\overline{U}_{ref}} = \left(\frac{z}{z_{ref}}\right)^{\alpha} \f]
 *
 * @param[in]  z     Height(s) in m at which you want to get speed.
 * @param[in]  u_ref Reference speed in m/s
 * @param[in]  z_ref Reference height in m
 * @param[in]  alpha Power law exponent
 */
template <typename T>
T power_law_speed(T const & z, const double u_ref, const double z_ref, const double alpha){
    T z_norm = z / z_ref;
    T speed = pow(z_norm, alpha) * u_ref;
    return speed;
};

// Remove template specialisations from doc (causes duplicate) @cond
template <>
Eigen::VectorXd power_law_speed(Eigen::VectorXd const & z, const double u_ref, const double z_ref, const double alpha) {
    // Template specialisation for VectorXd type
    Eigen::VectorXd z_norm = z / z_ref;
    Eigen::VectorXd speed = pow(z_norm.array(), alpha) * u_ref;
    return speed;
};
// @endcond


/** @brief Compute speed profile according to the MOST law.
 *
 * Templated so that it can be called with active scalars (allows use of autodiff), doubles/floats,
 * Eigen::Arrays (directly) or Eigen::VectorXds (via template specialisation) of z values.
 *
 * MOST law speed is computed as:
 *
 * TODO
 *
 * @param[in]  z        Height value(s) at which you want to get speed (m)
 * @param[in]  kappa    von Karman constant
 * @param[in]  d        Zero plane offset distance (e.g. for forest canopies) (m)
 * @param[in]  z0       Roughness length (m)
 * @param[in]  L        Monin-Obukhov length (m)
 *
 */
template <typename T>
T most_law_speed(T const & z, const double kappa, const double d, const double z0, const double L){
    std::cout << "MOST Law not implemented yet" << std::endl;
    T speed;
    return speed;
}


// Remove template specialisation from doc (causes duplicate) @cond
template <>
Eigen::VectorXd most_law_speed(Eigen::VectorXd const & z, const double kappa, const double d, const double z0, const double L){
    std::cout << "MOST Law not implemented yet" << std::endl;
    Eigen::VectorXd speed;
    return speed;
};
// @endcond


/** @brief Compute speed profile according to Marusic and Jones' relations.
 *
 * Templated so that it can be called with active scalars (allows use of autodiff), doubles/floats,
 * Eigen::Arrays (directly) or Eigen::VectorXds (via template specialisation) of z values.
 *
 * Speed is computed as:
 * @f[ \frac{\overline{U}}{U_{\tau}} & = & \frac{1}{\kappa} \ln \left( \frac{z+z_0}{k_s} \right) + Br + \frac{\Pi_j}{\kappa} W_c[\eta, \Pi_j] \\ W_c[\eta, \Pi_j] & = & 2 \eta^2 \left( 3 - 2\eta \right) - \frac{1}{3\Pi_j}\eta^3 \\ \eta & = & \frac{z+z_0}{\delta + z_0} @f]
 *
 * which reduces to the defecit relation:
 * @f[ U_{D}^{*} = \frac{U_\infty-\overline{U}}{U_{\tau}} = -\frac{1}{\kappa} \ln \left( \eta \right) + \frac{1}{3\kappa} \left(\eta^3 - 1 \right) + 2 \frac{\Pi_j}{\kappa} \left(1 - 3\eta^2 + 2\eta^3 \right) @f]
 *
 * @param[in]  z        Height(s) in m at which you want to get speed.
 * @param[in]  pi_j     Jones' modification of the Coles wake factor. Double or AutoDiffScalar type acccepted.
 * @param[in]  kappa    von Karman constant
 * @param[in]  z_0      Roughness length - represents distance of hypothetical smooth wall from actual rough wall z0 = 0.25k_s (m)
 * @param[in]  delta    Boundary layer thickness (m)
 * @param[in]  u_inf    Speed of flow at z = delta (m/s)
 * @param[in]  u_tau    Shear / skin friction velocity (governed by ratio parameter S = u_inf / u_tau)
 *
 */
template <typename T_z, typename T_pi_j>
T_z marusic_jones_speed(T_z const & z, T_pi_j const pi_j, const double kappa, const double z_0,
                      const double delta, const double u_inf, const double u_tau){
    T_z eta = (z + z_0) / (delta + z_0);
    T_z eta_cubed = pow(eta, 3.0);
    T_z term1 = log(eta) / kappa;
    T_z term2 = (eta_cubed - 1.0) / (3.0 * kappa);
    T_z term3 = 2.0 * pi_j * (1.0 - pow(eta, 2.0) * 3.0 + eta_cubed * 2.0) / kappa;
    T_z u_deficit = term2 - term1 + term3;
    T_z speed = u_inf - u_deficit * u_tau;
    return speed;
};
// Remove template specialisation from doc (causes duplicate) @cond
template <typename T_pi_j>
Eigen::VectorXd marusic_jones_speed(Eigen::VectorXd const & z, T_pi_j const pi_j, const double kappa, const double z_0,
                             const double delta, const double u_inf, const double u_tau){
    // Template specialisation for VectorXd type
    Eigen::VectorXd eta = (z.array() + z_0) / (delta + z_0);
    Eigen::VectorXd eta_cubed = eta.array().cube();
    Eigen::VectorXd term1 = eta.array().log() / kappa;
    Eigen::VectorXd term2 = (eta_cubed.array() - 1.0) / (3.0 * kappa);
    Eigen::VectorXd term3 = 2.0 * pi_j * (1.0 - eta.array().square() * 3.0 + eta_cubed.array() * 2.0) / kappa;
    Eigen::VectorXd u_deficit = term2 - term1 + term3;
    Eigen::VectorXd speed = u_inf - u_deficit.array() * u_tau;
    return speed;
};
//@endcond


/** @brief Compute coles wake parameter.
 *
 * Used by Perry and Marusic 1995.
 *
 * Templated so that it can be called with active scalars (allows use of autodiff), doubles/floats,
 * Eigen::Arrays (directly) or Eigen::VectorXds (via template specialisation) of z values.
 *
 * @f[ W_c[\eta, \Pi_j] & = & 2 \eta^2 \left( 3 - 2\eta \right) - \frac{1}{3\Pi_j}\eta^3 \\ \eta & = & \frac{z+z_0}{\delta + z_0} @f]
 *
 * Translated from MATLAB:
 * \code
 * function wc = colesWake(eta, Pi)
 *     wc = 2*eta.^2.*(3-2*eta) - (1/Pi).*eta.^2.*(1-eta).*(1-2*eta);
 * end
 * \endcode
 *
 * @param[in]  eta          Nondimensional height values
 * @param[in]  capital_pi   The coles wake parameter Pi
 *
 */
template <typename T>
T coles_wake(T const & eta, const double capital_pi){
    T wc, eta_sqd;
    eta_sqd = pow(eta, 2.0);
    wc = 2.0 * eta_sqd * (3.0 - 2.0 * eta)
        - eta_sqd * (1.0 - eta) * (1.0 - 2.0*eta) / capital_pi;
    return wc;
};
// Remove template specialisation from doc (causes duplicate) @cond
template <>
Eigen::VectorXd coles_wake(Eigen::VectorXd const & eta, const double capital_pi){
    Eigen::VectorXd wc, eta_sqd;
    eta_sqd = eta.array().pow(2.0);
    wc = 2.0 * eta_sqd.array() * (3.0 - 2.0 * eta.array())
        - eta_sqd.array() * (1.0 - eta.array()) * (1.0 - 2.0*eta.array()) / capital_pi;
    return wc;
};
//@endcond


/** Compute Lewkowicz (1982) velocity profile.
 *
 * Used by Perry and Marusic 1995 (from eqs 2 and 7).
 *
 * Templated so that it can be called with active scalars (allows use of autodiff), doubles/floats,
 * Eigen::Arrays (directly) or Eigen::VectorXds (via template specialisation) of z values.
 *
 * Translated from MATLAB:
 * \code
 * function [f] = getf(eta, Pi, kappa, S)
 *     f = (-1/kappa)*log(eta) + (Pi/kappa)*colesWake(1,Pi)*ones(size(eta)) - (Pi/kappa)*colesWake(eta,Pi);
 *     f(isinf(f)) = S;
 * end
 * \endcode
 *
 * @param[in]  eta Nondimensional (eta = z/delta_c) height value or values
 * @param[in]  pi_coles The coles wake parameter Pi
 * @param[in]  kappa von Karman constant
 * @param[in]  u_inf Speed of flow at z = delta (m/s)
 * @param[in]  u_tau Shear / skin friction velocity (governed by ratio parameter shear_ratio = u_inf / u_tau)
 *
 */
template <typename T>
T lewkowicz_speed(T const & eta, const double pi_coles, const double kappa, const double u_inf, const double u_tau) {
    T f, speed;
    f = pi_coles * coles_wake(1.0, pi_coles) / kappa;
    f = f - log(eta) / kappa;
    f = f - pi_coles * coles_wake(eta, pi_coles);
    // TODO sort this out so it can be template compliant
    //if (std::isinf(f)) {
    //    f = u_inf/u_tau;
    //}
    speed = u_inf - f*u_tau;
    return speed;
};
// Remove template specialisation from doc (causes duplicate) @cond
template <>
Eigen::VectorXd lewkowicz_speed(Eigen::VectorXd const & eta, const double pi_coles, const double kappa, const double u_inf, const double u_tau){
    Eigen::VectorXd f, speed;
    Eigen::VectorXd term1 = eta.array().log() / (-1.0*kappa);
    double term2 = pi_coles * coles_wake(1.0, pi_coles) / kappa;
    Eigen::VectorXd term3 = pi_coles * coles_wake(eta, pi_coles) / kappa;
    // std::cout << "term 1 = [" << term1.transpose() << "]" << std::endl;
    // std::cout << "term 2 = [" << term2 << "]" << std::endl;
    // std::cout << "term 3 = [" << term3.transpose() << "]" << std::endl;
    f = term1.array() + term2 - term3.array();
    for (int k = 0; k < f.size(); k++) {
        if (std::isinf(f[k])) {
            f(k) = u_inf / u_tau;
        }
    }
    speed = u_inf - f.array()*u_tau;
    return speed;
};
// @endcond


} /* namespace es */


#endif /* SOURCE_RELATIONS_VELOCITY_H_ */
