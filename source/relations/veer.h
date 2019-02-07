/*
 * veer.h Atmospheric Boundary Layer Veer relations
 *
 * Author:              Tom Clark  (thclark @ github)
 *
 * Copyright (c) 2015-9 Octue Ltd. All Rights Reserved.
 *
 */

#ifndef SOURCE_RELATIONS_VELOCITY_H_
#define SOURCE_RELATIONS_VELOCITY_H_

#include <Eigen/Dense>
#include <Eigen/Core>
#include "definitions.h"
#include "profile.h"


namespace es {


/** @brief Computes the left hand side of the veer relations given u(z) or v(z).
 *
 * Templated so that it can be called with active scalars (allows use of autodiff), doubles/floats,
 * Eigen::Arrays (directly) or Eigen::VectorXds (via template specialisation) of u values.
 *
 * The veer relations are:
 * @f[ 2 |\mathbf{\Omega}| sin(\phi) (\overline{v}_g - \overline{v}) & = & \frac{\partial}{\partial z} \bigg( \nu\frac{\partial\overline{u}}{\partial z} - \overline{u'w'}\bigg) \\ 2 |\mathbf{\Omega}| sin(\phi) (\overline{u}_g - \overline{u}) & = & \frac{\partial}{\partial z} \bigg(\nu\frac{\partial\overline{v}}{\partial z} - \overline{v'w'}\bigg) @f]
 *
 * @param[in]  ui       Mean velocity component at a given height (m/s)
 * @param[in]  ui_g     Mean geostrophic velocity component outside the atmospheric boundary layer (m/s)
 * @param[in]  phi      Latitude (degrees)
 */
template <typename T>
T veer_lhs(T const & ui, const double ui_g, const double phi){
    T lhs = 2.0*omega_world*sind(phi)*(u_g - u);
    return lhs;
}


/** @brief Computes the right hand side of the veer relations given u(z) or v(z).
 *
 * See `veer_lhs()` for the full relation.
 *
 * Templated so that it can be called with active scalars (allows use of autodiff), doubles/floats,
 * Eigen::Arrays (directly) or Eigen::VectorXds (via template specialisation) of u values.
 *
 * @param[in]  ui       Mean velocity component at a given height (m/s){
 * @param[in]  uiu3_bar Mean cross term of unsteady velocity components \f$ \overline{u_i'u_3'}} \f$ (m^2/s^2)
 * @param[in]  nu       Kinematic viscosity of the fluid (m^2/s)
 */
template <typename T>
T veer_rhs(T & ui, T & uiu3_bar, const double nu){
    T lhs, mixing_term;
    mixing_term = nu*ui.getZDerivative() + uiu3_bar;
    rhs = mixing_term.getZDerivative();
    return rhs;
}

} /* namespace es */

#endif /* SOURCE_RELATIONS_VELOCITY_H_ */
