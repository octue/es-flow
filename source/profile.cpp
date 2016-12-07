
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

    Bins::Bins(int n, double *z, double *dx, double *dy, double *dz): Bins(n, z), dx(dx), dy(dy), dz(dz) {
        // Construct from a spacing and bin size in x, y, z
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

    Bins::Bins(int n, double *z, double half_angle_degrees, double *dz): Bins(n, z, half_angle_degrees), dz(dz) {
        // Construct from spacing, with bin height given and bin widths using a half angle
        int i;
        double tha;
        tha = tand(half_angle_degrees);
        for (i = 0; i < n - 1; i++) {
            dx[i] = 2.0 * z[i] * tha;
            dy[i] = dx[i];
        }
    }

    Bins::~Bins() {
        //Destructor
    }

    Profile::Profile(Bins b): bins(b) {
        // Construct from just bins, with default position {0.0, 0.0, 0.0}

        // Initialise position to the default
        position[0] = 0.0;
        position[1] = 0.0;
        position[2] = 0.0;

    }

    // Construct from bins and a position
    Profile::Profile(Bins b, double x, double y, double z): bins(b) {

        // Initialise position to the default 0,0,0
        position[0] = x;
        position[1] = y;
        position[2] = z;

    }

    Profile::~Profile() {
        // Destructor
    }

    template <class ProfileType>
    void Profile::setValues(const std::vector<ProfileType> &values) {
        Profile::values = values;
    }

    template <class ProfileType>
    const std::vector<ProfileType> &Profile::getValues() const {
        return values;
    }

} /* namespace es */
