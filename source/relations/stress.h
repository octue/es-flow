/*
 * stress.h Atmospheric Boundary Layer Reynolds Stress relations
 *
 * Author:              Tom Clark  (thclark @ github)
 *
 * Copyright (c) 2019 Octue Ltd. All Rights Reserved.
 *
 */

#ifndef SOURCE_RELATIONS_STRESS_H_
#define SOURCE_RELATIONS_STRESS_H_

#include <Eigen/Dense>
#include <Eigen/Core>
#include <iostream>
#include <iomanip>
#include <math.h>
#include "NumericalIntegration.h"
#include "profile.h"
#include "utilities/trapz.h"
#include "utilities/cumtrapz.h"

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


    /** @brief Get the Coles wake distribution @f$ w_{c} @f$ from the wake parameter @f$ \Pi @f$
     *
     */
    template<typename T>
    T get_coles_wake(const T eta, const double pi_coles) {
        T wc;
        T eta_pow_2 = pow(eta, 2.0);
        wc = (2.0 * eta_pow_2 * (3.0 - (2.0 * eta)));
        wc -= (1.0 / pi_coles) * eta_pow_2 * (1.0 - eta) * (1.0 - (2.0 * eta));
        return wc;
    }
    Eigen::ArrayXd get_coles_wake(const Eigen::ArrayXd &eta, const double pi_coles) {
        Eigen::ArrayXd wc;
        Eigen::ArrayXd eta_pow_2;
        eta_pow_2 = eta.pow(2.0);
        wc = (2.0 * eta_pow_2 * (3.0 - (2.0 * eta)));
        wc -= (1.0 / pi_coles) * eta_pow_2 * (1.0 - eta) * (1.0 - (2.0 * eta));
        return wc;
    }



/** @brief Get the integrand @f$ f @f$ used in computation of @f$ R_{13} @f$ in the modified Lewkowicz method
     *
     * Uses equations 2 and 7 Perry and Marusic 1995 Part 1
     *
     */
    Eigen::ArrayXd get_f_integrand(const Eigen::ArrayXd eta, const double kappa, const double pi_coles, const double shear_ratio) {

        Eigen::ArrayXd f;
        Eigen::ArrayXd ones;
        ones.setOnes(eta.size());
        f = (-1.0 / kappa) * eta.log()
            + (pi_coles/kappa) * get_coles_wake(ones, pi_coles) * ones
            - (pi_coles/kappa) * get_coles_wake(eta, pi_coles);
        for (int k = 0; k < f.size(); k++) {
            if (std::isinf(f[k])) {
                f(k) = shear_ratio; // from eqs 2 and 7 P&M part 1 1995
            }
        }

        return f;
    }


    /** @brief Compute Horizontal-vertical Reynolds Stress R13 profile using modified Lewkowicz formulation
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

        // Get f between eta = 0 and eta = 1 (bounds for C1 integration)
        Eigen::ArrayXd f = get_f_integrand(eta, kappa, pi_coles, shear_ratio);
        const double d_pi = 0.01 * pi_coles;
        Eigen::ArrayXd f_plus  = get_f_integrand(eta, kappa, (pi_coles + d_pi), shear_ratio);
        Eigen::ArrayXd f_minus  = get_f_integrand(eta, kappa, (pi_coles - d_pi), shear_ratio);

        // TODO can we cast this returned value directly?
        Eigen::ArrayXd c1_tmp;
        trapz(c1_tmp, eta, f);
        const double c1 = c1_tmp(0);

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
        Eigen::ArrayXd a1 = e2 - e3 + shear_ratio * (e4 - e1);
        Eigen::ArrayXd a2 = (-2.0 * e2) + e3 + (shear_ratio * e1);
        Eigen::ArrayXd a3 = e6 - e7 - (shear_ratio * e5);
        Eigen::ArrayXd a4 = (2.0 * e2) - e3 + (shear_ratio * e4) - (shear_ratio * 3.0 * e1);

        // B coefficients are simply evaluated at eta = 1 (eqn A6)
        // TODO eta must be defined up to 1! check this
        double b1 = a1(eta.rows()-1);
        double b2 = a2(eta.rows()-1);
        double b3 = a3(eta.rows()-1);
        double b4 = a4(eta.rows()-1);

        // E1 from (eqn A4). Can't call it E1 due to name conflict with above.
        double e1_coeff = 1.0 / (kappa * shear_ratio + 1.0);

        // N from (eqn A5) using central differencing as before
        double wc_minus = get_coles_wake(1.0, pi_coles - d_pi);
        double wc_plus =  get_coles_wake(1.0, pi_coles + d_pi);
        double n = get_coles_wake(1.0, pi_coles) + pi_coles * (wc_plus - wc_minus) / (2.0 * d_pi);

        // Compile f_i terms
        Eigen::ArrayXd f1;
        Eigen::ArrayXd f2;
        Eigen::ArrayXd f3;
        f1 = 1
            - a1 / (b1 + e1_coeff * b2)
            - e1_coeff * a2 / (b1 + e1_coeff * b2);

        f2 = (e1_coeff * n * a2 * b1
            + a3 * b1
            - e1_coeff * n * a1 * b2
            + e1_coeff * a3 * b2
            - a1 * b3
            - e1_coeff * a2 * b3) / (b1 + e1_coeff * b2);

        f3 = (e1_coeff * a2 * b1 + a4 * b1 - e1_coeff * a1 * b2 + e1_coeff * a4 * b2 - a1 * b4 - e1_coeff * a2 * b4) / (b1 + e1_coeff * b2);

        // Convert f2 and f3 into g1 and g2 ready for application of eq. 15
        Eigen::ArrayXd g1 = f2 / shear_ratio;
        Eigen::ArrayXd g2 = -f3 / (c1 * shear_ratio);

        // Top 3 boxes of figure 18
        Eigen::ArrayXd minus_r13 = f1 + g1 * zeta + g2 * beta;

        // Get the component due to equilibrium sink flow (OLD version - see P&M eqns 51,53)
        // minusReynStressA = ones(size(eta)) - eta + eta.*log(eta);

        // Lewkowicz 1982 shear stress for equilibrium sink flow, Perry and Marusic eqn. 51
        Eigen::ArrayXd minus_r13_a = 1.0 - (60.0/59.0)*eta - (20.0/59.0)*eta.pow(3.0) + (45.0/59.0)*eta.pow(4.0) - (24.0/59.0)*eta.pow(5.0) + (60.0/59.0)*eta*eta.log();

        // Handle the log(0) singularity
        for (int i = 0; i < minus_r13_a.size(); i++) {
            if (std::isinf(minus_r13_a(i))) {
                minus_r13_a(i) = 1.0;
            }
        }

        // Output correctly signed reynolds stress components
        r13_a = -1.0 * minus_r13_a;
        r13_b = -1.0 * (minus_r13 - minus_r13_a);

	}

} /* namespace es */

#endif /* SOURCE_RELATIONS_VELOCITY_H_ */
