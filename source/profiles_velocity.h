/*
 * profiles_velocity.h
 *
 *  Created on: 29 Nov 2016
 *      Author: thc29
 */

#ifndef SOURCE_PROFILES_VELOCITY_H_
#define SOURCE_PROFILES_VELOCITY_H_

#include <Eigen/Dense>
#include <Eigen/Core>
#include "profile.h"

using namespace Eigen;

namespace es {

	template <typename T>
	T power_law_speed(T const & z, const double u_ref, const double z_ref, const double alpha){
		/// Compute speed profile according to the power law.
        /**
		 *
		 * u_ref Reference speed in m/s
		 * z_ref Reference height in m
		 * alpha Power law exponent
		 *
		 * Templated so that it can be called with active scalars (allows use of autodiff), doubles/floats,
		 * Eigen::Arrays (directly) or Eigen::VectorXds (via template specialisation) of z values
		 */
		T z_norm = z / z_ref;
		T speed = pow(z_norm, alpha) * u_ref;
		//    std::cout << speed << std::endl;
		return speed;
	}

	template <>
	VectorXd power_law_speed(VectorXd const & z, const double u_ref, const double z_ref, const double alpha) {
		// Template specialisation for VectorXd type
		VectorXd z_norm = z / z_ref;
		VectorXd speed = pow(z_norm.array(), alpha) * u_ref;
		return speed;
	};

    template <typename T>
    T most_law_speed(T const & z, const double kappa, const double d, const double z0, const double L){

    /* Compute speed profile according to the MOST law.
     * von karman constant
     * zero plane offset distance (e.g. for forest canopies)
     * roughness length
     * Monin-Obukhov length
     */
		std::cout << "MOST Law not implemented yet" << std::endl;
		T speed;
		return speed;
    }

	template <typename T>
	T marusic_jones_speed(T const & z, const double u_ref, const double z_ref, const double alpha){
		/* Compute speed profile according to marusic/jones.
		 *
		 * u_ref Reference speed in m/s
		 * z_ref Reference height in m
		 * alpha Power law exponent
		 *
		 * Templated so that it can be called with active scalars (allows use of autodiff), doubles/floats,
		 * Eigen::Arrays (directly) or Eigen::VectorXds (via template specialisation) of z values
		 */
		T z_norm = z / z_ref;
		T speed = pow(z_norm, alpha) * u_ref;
		//    std::cout << speed << std::endl;
		return speed;
	}

	template <>
	VectorXd marusic_jones_speed(VectorXd const & z, const double u_ref, const double z_ref, const double alpha) {
		// Template specialisation for VectorXd type
		VectorXd z_norm = z / z_ref;
		VectorXd speed = pow(z_norm.array(), alpha) * u_ref;
		return speed;
	};


	class VelocityProfile: public Profile<double> {
	public:
		VelocityProfile();
		virtual ~VelocityProfile();

		VectorXd AutoDiff();
	};

} /* namespace es */

#endif /* SOURCE_PROFILES_VELOCITY_H_ */
