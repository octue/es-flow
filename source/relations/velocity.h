/*
 * velocity.h Atmospheric Boundary Layer velocity profile relations
 *
 * Author:              Tom Clark  (thclark @ github)
 *
 * Copyright (c) 2013-9 Octue Ltd. All Rights Reserved.
 *
 */

#ifndef ES_FLOW_VELOCITY_H_
#define ES_FLOW_VELOCITY_H_

#include <Eigen/Dense>
#include <Eigen/Core>
#include "profile.h"
#include <iostream>


namespace es {

/** @brief Compute speed profile according to the power law.
 *
 * Templated so that it can be called with active scalars (allows use of autodiff), doubles/floats,
 * Eigen::Arrays (directly) or Eigen::VectorXds (via template specialisation) of z values (to obtain shear profile
 * u/dz) and/or alpha values (to obtain variational jacobean for fitting parameter alpha).
 *
 * Power law speed is computed as:
 * \f[ \frac{\overline{U}}{\overline{U}_{ref}} = \left(\frac{z}{z_{ref}}\right)^{\alpha} \f]
 *
 * @param[in]  z     Height(s) in m at which you want to get speed.
 * @param[in]  u_ref Reference speed in m/s
 * @param[in]  z_ref Reference height in m
 * @param[in]  alpha Power law exponent. Must be of same type as input z (allows autodifferentiation).
 */
template <typename T, typename Talf>
T power_law_speed(T const & z, const double u_ref, const double z_ref, Talf const & alpha){
    T z_norm = z / z_ref;
    T speed = pow(z_norm, alpha) * u_ref;
    return speed;
};

// Remove template specialisations from doc (causes duplicate) @cond
Eigen::VectorXd power_law_speed(Eigen::VectorXd const & z, const double u_ref, const double z_ref, const double alpha) {
    Eigen::VectorXd z_norm = z.array() / z_ref;
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
 * @param[in]  pi_j     Jones' modification of the Coles wake factor. Double or AutoDiffScalar type accepted.
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


/** @brief Compute corner-corrected coles wake parameter @f W_c @f.
 *
 * There are two alternative functional forms for the wake parameter. Perry and Marusic (1995) used the
 * Lewkowicz (1982) wake function, where the second (@f 1/\Pi... @f) term is artificially applied to ensure that
 * the velocity defect function has zero gradient at @f z= \delta_c @f:
 *
 * @f[ W_c[\eta, \Pi] & = & 2 \eta^2 \left( 3 - 2\eta \right) - \frac{1}{\Pi}\eta^2 \left( 1 - \eta \right) \left( 1 - 2\eta \right) @f]
 *
 * Note: Jones, Marusic and Perry 2001 refer to this correction as the 'corner' function and use an alternative
 * correction (recommended by Prof. Coles) to the pure wall flow (instead of the wake) which is not a function of Pi,
 * allowing reversion to the original (Coles 1956) wake parameter @f W @f:
 *
 * @f[ W[\eta, \Pi] & = & 2 \eta^2 \left( 3 - 2\eta \right) @f]
 *
 * This function implements the former by default!
 *
 * This function is templated so that it can be called with active scalars (allows use of autodiff), doubles/floats,
 * Eigen::Arrays (directly) or Eigen::VectorXds (the latter via template specialisation) of z values.
 *
 * @param[in]  eta          Nondimensional height values
 * @param[in]  pi_coles     The coles wake parameter Pi
 * @param[in]  corrected    Boolean flag, default true. If true, return Lewkowicz (1982) corrected wake function. If false, return the original uncorrected wake parameter.
 *
 */
template <typename T_z, typename T_param>
T_z coles_wake(T_z const &eta, T_param const &pi_coles, const bool corrected=true){
    T_z wake_param, eta_sqd;
    eta_sqd = pow(eta, 2.0);
    if (corrected) {
        wake_param = 2.0 * eta_sqd * (3.0 - 2.0 * eta)
            - eta_sqd * (1.0 - eta) * (1.0 - 2.0 * eta) / pi_coles;
    } else {
        wake_param = 2.0 * eta_sqd * (3.0 - 2.0 * eta);
    }
    return wake_param;
};

// Remove template specialisation from doc (causes duplicate) @cond
template <typename T_param>
Eigen::VectorXd coles_wake(Eigen::VectorXd const &eta, T_param const &pi_coles, const bool corrected=true){
    Eigen::VectorXd wake_param, eta_sqd;
    eta_sqd = eta.array().pow(2.0);
    if (corrected) {
        wake_param = 2.0 * eta_sqd.array() * (3.0 - 2.0 * eta.array())
            - eta_sqd.array() * (1.0 - eta.array()) * (1.0 - 2.0*eta.array()) / pi_coles;
    } else {
        wake_param = 2.0 * eta_sqd.array() * (3.0 - 2.0 * eta.array());
    }
    return wake_param;
};
//@endcond


// Do not document @cond
/* Template wrapper for std::isinf to kill off problems where ceres::Jet is used (autodifferentiation) instead of
 * a double. Call it is_dbl_inf rather than isinf to avoid accidental use.
 */
bool is_double_inf(double in) {
    return std::isinf(in);
};

template <typename T>
bool is_double_inf(T in) {
    return false;
};
// @endcond


/** @brief Get the velocity deficit integrand @f$ f @f$ used in computation of @f$ R_{13} @f$.
 *
 * TODO move into velocity.h, reuse in lewkowicz_speed and jones_speed functions;.
 * In their original derivation of the ADEM, Perry and Marusic 1995 eqn. 2 use the velocity distribution
 * of Lewkowicz 1982 as a basis to determine the shear stress profile. This led to an integrand for the mean velocity
 * deficit profile:
 *
 * @f[ f & = & \frac{U_{1} - \overline{U}}{U_{\tau}} \\ & = & - \frac{1}{\kappa} \ln \left( \eta \right) + \frac{\Pi}{\kappa} W_c[1, \Pi] - \frac{\Pi}{\kappa} W_c[\eta, \Pi], \\ W_c[\eta, \Pi] & = & 2 \eta^2 \left( 3 - 2\eta \right) - \frac{1}{\Pi}\eta^2 \left( 1 - \eta \right) \left( 1 - 2\eta \right) @f]
 *
 * Jones, Marusic and Perry (2001) eqn 1.9 uses an alternative formulation, thereby removing the non-physical dependency
 * of the wake factor on @f \Pi @f. This leads to the velocity deficit function:
 *
 * @f[ f & = & \frac{U_{1} - \overline{U}}{U_{\tau}} \\ & = & - \frac{1}{\kappa} \ln \left( \eta \right) + \frac{1}{3\kappa} \left(\eta^{3} - 1 \right) + 2 \frac{\Pi}{\kappa} \left( 1 - 3\eta^{2} + 2\eta^{3} \right) @f]
 *
 * This function is templated so that it can be called with active scalars (allows use of autodiff), doubles/floats,
 * Eigen::Arrays (directly) or Eigen::VectorXds (the latter via template specialisation) of z values.
 *
 * @param[in] eta           Nondimensional height values
 * @param[in] kappa         von Karman constant
 * @param[in] pi_coles      The coles wake parameter Pi
 * @param[in] shear_ratio   Shear ratio S
 * @param[in] lewkowicz     Boolean flag, default false. If true, use Lewkowicz (1982) velocity deficit funciton. If false, use Jones et at (2001).
 *
 * @return
 */
Eigen::ArrayXd deficit(const Eigen::ArrayXd &eta, const double kappa, const double pi_coles, const double shear_ratio, const bool lewkowicz=false) {
    Eigen::ArrayXd f;
    Eigen::ArrayXd ones;
    ones.setOnes(eta.size());
    f = -1.0 * eta.log()/ kappa;
    if (lewkowicz) {
        f = f + (pi_coles/kappa) * coles_wake(ones, pi_coles)
              - (pi_coles/kappa) * coles_wake(eta, pi_coles);
    } else {
        f = f + (eta.pow(3.0) - 1.0)/(3.0*kappa)
              + 2.0*pi_coles*(1.0 - 3.0*eta.pow(2.) + 2.0*eta.pow(3.0))/kappa;
    }
    for (int k = 0; k < f.size(); k++) {
        if (std::isinf(f[k])) {
            f(k) = shear_ratio;
        }
    }
    return f;
}


/** Compute Lewkowicz (1982) velocity profile.
 *
 * TODO refactor to base it on velocity deficit, keep code DRY
 * TODO Refactor to find usages as VectorXd and change them to array uses; remove the VectorXd template specialization
 *
 * Used by Perry and Marusic 1995 (from eqs 2 and 7).
 *
 * Templated so that it can be called with active scalars (allows use of autodiff), doubles/floats,
 * Eigen::Arrays (directly) or Eigen::VectorXds (via template specialisation) of z values.
 *
 * @param[in]  z height in m (or Nondimensional heights (eta = z/delta_c) where delta_c = 1.0)
 * @param[in]  pi_coles The coles wake parameter Pi
 * @param[in]  kappa von Karman constant
 * @param[in]  u_inf Speed of flow at z = delta (m/s)
 * @param[in]  shear_ratio Shear / skin friction velocity ratio (shear_ratio = u_inf / u_tau)
 * @param[in]  delta_c Boundary layer thickness in m, used to normalise z. Defaults to 1.0.
 *
 */
template <typename T_z, typename T_param>
T_z lewkowicz_speed(T_z const & z, T_param const & pi_coles, T_param const & kappa, T_param const & u_inf, T_param const & shear_ratio, T_param const &delta_c) {
    T_z f, speed, eta;
    eta = z / delta_c;
    T_param u_tau = u_inf / shear_ratio;
    T_z term1 = log(eta) / (-1.0*kappa);
    T_z term2 = pi_coles * coles_wake(T_z(1.0), pi_coles, true) / kappa;
    T_z term3 = pi_coles * coles_wake(eta, pi_coles, true) / kappa;
    f = term1 + term2 - term3;
    if (is_double_inf(f)) {
        f = u_inf / u_tau;
    };
    speed = u_inf - f * u_tau;
    return speed;
};

// Remove template specialisation from doc (causes duplicate) @cond
template <typename T_param>
Eigen::VectorXd lewkowicz_speed(Eigen::VectorXd const & z, T_param const &pi_coles, T_param const &kappa, T_param const &u_inf, T_param const &shear_ratio, T_param const &delta_c=1.0){
    Eigen::VectorXd f, speed, eta;
    eta = z.array() / delta_c;
    T_param u_tau = u_inf/shear_ratio;
    Eigen::VectorXd term1 = eta.array().log() / (-1.0*kappa);
    double term2 = pi_coles * coles_wake(1.0, pi_coles) / kappa;
    Eigen::VectorXd term3 = pi_coles * coles_wake(eta, pi_coles).array() / kappa;
    f = term1.array() + term2 - term3.array();
    for (int k = 0; k < f.size(); k++) {
        if (std::isinf(f[k])) {
            f(k) = u_inf / u_tau;
        }
    }
    speed = u_inf - f.array() * u_tau;
    return speed;
};

template <typename T_param>
Eigen::ArrayXd lewkowicz_speed(Eigen::ArrayXd const & z, T_param const &pi_coles, T_param const &kappa, T_param const &u_inf, T_param const &shear_ratio, T_param const &delta_c=1.0){
    Eigen::ArrayXd f, speed, eta;
    eta = z / delta_c;
    T_param u_tau = u_inf/shear_ratio;
    Eigen::ArrayXd term1 = eta.log() / (-1.0*kappa);
    double term2 = pi_coles * coles_wake(1.0, pi_coles) / kappa;
    Eigen::ArrayXd term3 = pi_coles * coles_wake(eta, pi_coles) / kappa;
    f = term1 + term2 - term3;
    for (int k = 0; k < f.size(); k++) {
        if (std::isinf(f[k])) {
            f(k) = u_inf / u_tau;
        }
    }
    speed = u_inf - f * u_tau;
    return speed;
};
// @endcond


} /* namespace es */


#endif /* ES_FLOW_VELOCITY_H_ */
