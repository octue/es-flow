/*
 * velocity.h
 *
 */

#ifndef SOURCE_RELATIONS_VELOCITY_H_
#define SOURCE_RELATIONS_VELOCITY_H_

#include <Eigen/Dense>
#include <Eigen/Core>
#include "profile.h"

using namespace Eigen;

namespace es {


    /// Compute speed profile according to the power law.
    /**
     * Templated so that it can be called with active scalars (allows use of autodiff), doubles/floats,
     * Eigen::Arrays (directly) or Eigen::VectorXds (via template specialisation) of z values.
     *
     * Power law speed is computed as:
     *
     * \f{eqnarray*}{\frac{\overline{U}}{\overline{U}_{ref}} & = & \left(\frac{z}{z_{ref}}\right)^{\alpha} \f}
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
	}

	template <>
	VectorXd power_law_speed(VectorXd const & z, const double u_ref, const double z_ref, const double alpha) {
		// Template specialisation for VectorXd type
		VectorXd z_norm = z / z_ref;
		VectorXd speed = pow(z_norm.array(), alpha) * u_ref;
		return speed;
	};

    /// Compute speed profile according to the MOST law.
    /**
     * Templated so that it can be called with active scalars (allows use of autodiff), doubles/floats,
     * Eigen::Arrays (directly) or Eigen::VectorXds (via template specialisation) of z values.
     *
     * MOST law speed is computed as:
     *
     * TODO
     *
     * @param[in]  z        Height(s) in m at which you want to get speed.
     * @param[in]  kappa    von Karman constant
     * @param[in]  d        Zero plane offset distance (e.g. for forest canopies) (m)
     * @param[in]  z0       Roughness length (m)
     * @param[in]  L        Monin-Obukhov length (m)
     *
     */
    template <typename T>
    T most_law_speed(T const & z, const double kappa, const double d, const double z0, const double L){

    /* Compute speed profile according to the MOST law.
     */
		std::cout << "MOST Law not implemented yet" << std::endl;
		T speed;
		return speed;
    }

    /// Compute speed profile according to Maursic's and Jones' relations.
    /**
     * Templated so that it can be called with active scalars (allows use of autodiff), doubles/floats,
     * Eigen::Arrays (directly) or Eigen::VectorXds (via template specialisation) of z values.
     *
     * Speed is computed as:
     *
     * TODO
     *
     * @param[in]  z        Height(s) in m at which you want to get speed.
     * @param[in]  pi_j     Jones' modification of the Coles wake factor
     * @param[in]  kappa    von Karman constant
     * @param[in]  z_0       Roughness length - represents distance of hypothetical smooth wall from actual rough wall z0 = 0.25k_s (m)
     * @param[in]  delta    Boundary layer thickness (m)
     * @param[in]  u_inf    Speed of flow at z = delta (m/s)
     * @param[in]  u_tau    Shear / skin friction velocity (governed by ratio parameter S = u_inf / u_tau)
     *
     */
	template <typename T>
	T marusic_jones_speed(T const & z, const double pi_j, const double kappa, const double z_0,
                          const double delta, const double u_inf, const double u_tau){
        T eta = (z + z_0) / (delta + z_0);
        T eta_cubed = pow(eta, 3.0);
        T term1 = log(eta) / kappa;
        T term2 = (eta_cubed - 1.0) / (3.0 * kappa);
        T term3 = 2.0 * pi_j * (1.0 - pow(eta, 2.0) * 3.0 + eta_cubed * 2.0) / kappa;
        T u_deficit = term2 - term1 + term3;
        T speed = u_inf - u_deficit * u_tau;
        return speed;
    }

	template <>
	VectorXd marusic_jones_speed(VectorXd const & z, const double pi_j, const double kappa, const double z_0,
                                 const double delta, const double u_inf, const double u_tau){
		// Template specialisation for VectorXd type
        VectorXd eta = (z.array() + z_0) / (delta + z_0);
        VectorXd eta_cubed = eta.array().cube();
        VectorXd term1 = eta.array().log() / kappa;
        VectorXd term2 = (eta_cubed.array() - 1.0) / (3.0 * kappa);
        VectorXd term3 = 2.0 * pi_j * (1.0 - eta.array().square() * 3.0 + eta_cubed.array() * 2.0) / kappa;
        VectorXd u_deficit = term2 - term1 + term3;
        VectorXd speed = u_inf - u_deficit.array() * u_tau;
        return speed;
	};



} /* namespace es */

#endif /* SOURCE_RELATIONS_VELOCITY_H_ */
