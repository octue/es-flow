
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

    Bins::Bins(int n, double *z): n_bins(n), z_ctrs(z) {
        // Construct from a spacing with default bin size from central differencing the spacing
        int i;
        dz[0] = fabs(z[1]-z[0]);
        for (i = 1; i < n-1; i++) {
            dz[i] = 0.5*fabs(z[i+1]-z[i-1]);
        }
        dz[n-1] = fabs(z[n-1]-z[n-2]);
        for (i = 0; i < n-1; i++) {
            dx[i] = dz[i];
            dy[i] = dz[i];
        }
    }

    Bins::Bins(int n, double *z, double *d_x, double *d_y, double *d_z): Bins(n, z) {
        // Construct from a spacing and bin size in x, y, z
        dx = d_x;
        dy = d_y;
        dz = d_z;
    }

    Bins::Bins(int n, double *z, double half_angle_degrees): Bins(n, z) {
        // Construct from spacing, with bin height from central differencing of the spacing and bin widths using a half angle
        int i;
        double tha;
        tha = tand(half_angle_degrees);
        for (i = 0; i < n - 1; i++) {
            dx[i] = 2.0 * z[i] * tha;
            dy[i] = dx[i];
        }
    }

    Bins::Bins(int n, double *z, double half_angle_degrees, double *d_z): Bins(n, z, half_angle_degrees) {
        // Construct from spacing, with bin height given and bin widths using a half angle
        int i;
        double tha;
        tha = tand(half_angle_degrees);
        double *dx = NULL;
        double *dy = NULL;
        double *dz = NULL;
        dx = new double[n];
        dy = new double[n];
        dz = new double[n];
        for (i = 0; i < n - 1; i++) {
            dx[i] = 2.0 * z[i] * tha;
            dy[i] = dx[i];
            dz[i] = d_z[i];
        }
    }

    Bins::~Bins() {
        //Destructor
        delete [] dx;
        delete [] dy;
        delete [] dz;
    }

    template <class ProfileType>
    Profile<ProfileType>::Profile(const Bins &bins) : bins(bins) {
        // Construct from just bins, with default position {0.0, 0.0, 0.0}

        // Initialise position to the default
        position[0] = 0.0;
        position[1] = 0.0;
        position[2] = 0.0;

    }

    template <class ProfileType>
    Profile<ProfileType>::Profile(const Bins &bins, double x, double y, double z): bins(bins) {
        // Construct from bins and a position

        // Initialise position to the default 0,0,0
        position[0] = x;
        position[1] = y;
        position[2] = z;

    }

    template <class ProfileType>
    Profile<ProfileType>::~Profile() {
        // Destructor
    }

    template <class ProfileType>
    void Profile<ProfileType>::setValues(const std::vector<ProfileType> &values) {
        Profile::values = values;
    }

    template <class ProfileType>
    const std::vector<ProfileType> &Profile<ProfileType>::getValues() const {
        return values;
    }

} /* namespace es */
