/*
 * profiles_velocity.cpp
 *
 *  Created on: 29 Nov 2016
 *      Author: thc29
 */

#include <Eigen/Dense>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <unsupported/Eigen/AutoDiff>
#include "constants.h"
#include "profiles_velocity.h"

using namespace Eigen;

namespace es {





    /* Use profile base class constructor and destructor
    VelocityProfile::VelocityProfile() {
    // TODO Auto-generated constructor stub
    }

    VelocityProfile::~VelocityProfile() {
    // TODO Auto-generated destructor stub
    }
    */

    VectorXd VelocityProfile::AutoDiff() {
        // Differentiate the velocity profile wrt z

	    typedef Eigen::AutoDiffScalar<Eigen::VectorXd> AScalar;
	    // AScalar stores a scalar and a derivative vector.

	    // Instantiate an AutoDiffScalar variable with a normal Scalar
	    double s = 0.3;
	    AScalar As(s);

        // Get the value from the Instance
        std::cout << "value: " << As.value() << std::endl;

        // The derivative vector is As.derivatives();

        // Resize the derivative vector to the number of dependent variables
        As.derivatives().resize(2);

        // Set the initial derivative vector
        As.derivatives() = Eigen::VectorXd::Unit(2,0);
        std::cout << "Derivative vector : " << As.derivatives().transpose() << std::endl;

        // Instantiate another AScalar
        AScalar Ab(4);
        Ab.derivatives() = Eigen::VectorXd::Unit(2,1);

        // Do the most simple calculation
        AScalar Ac = As * Ab;

        std::cout << "Result/Ac.value()" << Ac.value() << std::endl;
        std::cout << "Gradient: " << Ac.derivatives().transpose() << std::endl;

        Eigen::VectorXd a(1);
        return a;
    }

} /* namespace es */
