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



    /// Compute Reynolds Stress R13 profile according to Perry & Marusic 1995
    /**
     * Templated so that it can be called with active scalars (allows use of autodiff), doubles/floats,
     * Eigen::Arrays (directly) or Eigen::VectorXds (via template specialisation) of z values.
     *
     * Power law speed is computed as:
     * \f{eqnarray*}{\frac{\overline{U}}{\overline{U}_{ref}} & = & \left(\frac{z}{z_{ref}}\right)^{\alpha} \f}
     *
     * @param[in]  z     Height(s) in m at which you want to get speed.
     * @param[in]  u_ref Reference speed in m/s
     * @param[in]  z_ref Reference height in m
     * @param[in]  alpha Power law exponent
     */
	template <typename T>
	T r13(T const & z, const double u_ref, const double z_ref, const double alpha){
		T r13;
		return r13;
	}

	template <>
	VectorXd r13(VectorXd const & z, const double u_ref, const double z_ref, const double alpha) {
		// Template specialisation for VectorXd type
        VectorXd r13;
		return r13;
	};



} /* namespace es */

#endif /* SOURCE_RELATIONS_VELOCITY_H_ */
