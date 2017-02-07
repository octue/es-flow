
#include "profile.h"
#include "math.h"

// TODO move out to general purpose libraries
double sind(double angle) {
    double angle_radians = angle * M_PI / 180.0f;
    return sin(angle_radians) * M_PI / 180.0f;
}
double tand(double angle) {
    double angle_radians = angle * M_PI / 180.0f;
    return tan(angle_radians) * M_PI / 180.0f;
}

namespace es {

    Bins::Bins(std::vector<double> z): z_ctrs(z) {
        // Construct from a spacing with default bin size from central difference of the spacing
        int i;
        n_bins = z_ctrs.size();
        dx = std::vector<double>(n_bins);
        dy = std::vector<double>(n_bins);
        dz = std::vector<double>(n_bins);
        dz[0] = fabs(z[1]-z[0]);
        for (i = 1; i < n_bins-1; i++) {
            dz[i] = 0.5*fabs(z[i+1]-z[i-1]);
        }
        dz[n_bins-1] = fabs(z[n_bins-1]-z[n_bins-2]);
        for (i = 0; i < n_bins-1; i++) {
            dx[i] = dz[i];
            dy[i] = dz[i];
        }
    }

    Bins::Bins(std::vector<double> z, std::vector<double> d_x, std::vector<double> d_y, std::vector<double> d_z): Bins(z) {
        // Construct from a spacing and bin size in x, y, z
        dx = d_x;
        dy = d_y;
        dz = d_z;
        if ((dx.size() != n_bins) ||  (dy.size() != n_bins) || (dz.size() != n_bins)){
            throw std::range_error("Length of dx, dy or dz does not match the number of bins");
        }
    }

    Bins::Bins(std::vector<double> z, double half_angle_degrees): Bins(z) {
        // Construct from spacing, with bin height from central differencing of the spacing and bin widths using a half angle
        int i;
        double tha;
        tha = tand(half_angle_degrees);
        dx = std::vector<double>(z_ctrs.size());
        dy = std::vector<double>(z_ctrs.size());
        for (i = 0; i < n_bins - 1; i++) {
            dx[i] = 2.0 * z[i] * tha;
            dy[i] = dx[i];
        }
    }

    Bins::Bins(std::vector<double> z, double half_angle_degrees, std::vector<double> d_z): Bins(z, half_angle_degrees) {
        // Construct from spacing, with bin height given and bin widths using a half angle
        int i;
        dz = d_z;
        if (dz.size() != n_bins){
            throw std::range_error("Length of dx, dy or dz does not match the number of bins");
        }
    }

    Bins::~Bins() {
        //Destructor
    }

    ::std::ostream& operator<<(::std::ostream& os, const Bins& bins) {
        // Represent in logs or ostream
        return os << "debug statement for bins class";
    }


//    /* Use profile base class constructor and destructor
//    VelocityProfile::VelocityProfile() {
//    // TODO Auto-generated constructor stub
//    }
//
//    VelocityProfile::~VelocityProfile() {
//    // TODO Auto-generated destructor stub
//    }
//    */
//
//    VectorXd VelocityProfile::AutoDiff() {
//    // Differentiate the velocity profile wrt z
//
//    typedef Eigen::AutoDiffScalar<Eigen::VectorXd> AScalar;
//    // AScalar stores a scalar and a derivative vector.
//
//    // Instantiate an AutoDiffScalar variable with a normal Scalar
//    double s = 0.3;
//    AScalar As(s);
//
//    // Get the value from the Instance
//    std::cout << "value: " << As.value() << std::endl;
//
//    // The derivative vector is As.derivatives();
//
//    // Resize the derivative vector to the number of dependent variables
//    As.derivatives().resize(2);
//
//    // Set the initial derivative vector
//    As.derivatives() = Eigen::VectorXd::Unit(2,0);
//    std::cout << "Derivative vector : " << As.derivatives().transpose() << std::endl;
//
//    // Instantiate another AScalar
//    AScalar Ab(4);
//    Ab.derivatives() = Eigen::VectorXd::Unit(2,1);
//
//    // Do the most simple calculation
//    AScalar Ac = As * Ab;
//
//    std::cout << "Result/Ac.value()" << Ac.value() << std::endl;
//    std::cout << "Gradient: " << Ac.derivatives().transpose() << std::endl;
//
//    Eigen::VectorXd a(1);
//    return a;
//}

} /* namespace es */
