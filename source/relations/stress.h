/*
 * velocity.h
 *
 */

#ifndef SOURCE_RELATIONS_STRESS_H_
#define SOURCE_RELATIONS_STRESS_H_

#include <Eigen/Dense>
#include <Eigen/Core>
//#include <NumericalIntegration.h>
#include "profile.h"
#include <iostream>
#include <iomanip>
#include <unsupported/Eigen/NumericalIntegration>

using namespace Eigen;

namespace es {


    /// Integrator for the velocity deficit
    /**
     *  We consider the example from:
     *
     *       http://www.gnu.org/software/gsl/manual/html_node/Numerical-integration-examples.html
     *
     *       int_0^1 x^{-1/2} log(x) dx = -4
     *
     *  The integrator expects the user to provide a functor as shown below.
     */
    template<typename Scalar>
    class IntegrandExampleFunctor
    {
    public:
        IntegrandExampleFunctor(const Scalar alpha):m_alpha(alpha)
        {
            assert(alpha>0);
        }

        Scalar operator()(const Scalar x) const
        {
            assert(x>0);
            return log(m_alpha*x) / sqrt(x);
        }

        void setAlpha(const Scalar alpha)
        {
            m_alpha = alpha;
        }
    private:
        Scalar m_alpha;
    };


    double do_something(const double arg)
    {
        // Define the scalar
        typedef double Scalar;

        // Define the functor
        Scalar alpha=1.;
        IntegrandExampleFunctor<Scalar> inFctr(alpha);

        //define the integrator
        Eigen::Integrator<Scalar> eigIntgtor(200);

        //define a quadrature rule
        Eigen::Integrator<Scalar>::QuadratureRule quadratureRule = Eigen::Integrator<Scalar>::GaussKronrod61;

        //define the desired absolute and relative errors
        Scalar desAbsErr = Scalar(0.);
        Scalar desRelErr = Eigen::NumTraits<Scalar>::epsilon() * 50.;

        //integrate
        Scalar result = eigIntgtor.quadratureAdaptive(inFctr, Scalar(0.),Scalar(1.), desAbsErr, desRelErr, quadratureRule);

        //expected result
        Scalar expected = Scalar(-4.);

        //print output
        size_t outputPrecision  = 18;
        std::cout<<std::fixed;
        std::cout<<"result          = "<<std::setprecision(outputPrecision)<<result<<std::endl;
        std::cout<<"exact result    = "<<std::setprecision(outputPrecision)<<expected<<std::endl;
        std::cout<<"actual error    = "<<std::setprecision(outputPrecision)<<(expected-result)<<std::endl;

        return 0.0;
    }



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
