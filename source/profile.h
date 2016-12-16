#ifndef SOURCE_PROFILE_H
#define SOURCE_PROFILE_H

#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>
#include <stdexcept>

namespace es {

    class Bins {
    public:

        // Construct from a spacing with bin sizes set by central differencing of spacing
        Bins(int n, double *z);

        // Construct from spacing with bin sizes specified
        Bins(int n, double *z, double *dx, double *dy, double *dz);

        // Construct from spacing and a half angle in degrees, dz set by central differencing of spacing
        Bins(int n, double *z, double half_angle_degrees);

        // Construct from spacing and a half angle in degrees, with dz specified
        Bins(int n, double *z, double half_angle_degrees, double *dz);

        // Destroy
        ~Bins();

        // Vectors of z locations and bin dimensions
        int n_bins;
        double *z_ctrs;
        double *dx;
        double *dy;
        double *dz;

    };

    template <class ProfileType>
    class Profile {
    private:
        std::vector<ProfileType> values;

    public:

        // Construct using bins only (zero values and position)
        Profile(const Bins &bins);

        // Construct using bins and a global origin location
        Profile(const Bins &bins, double x, double y, double z);

        // Destroy
        virtual ~Profile();

        // Get the values, whatever type they might be
//        ProfileType *getValues() const;
        const std::vector<ProfileType> &getValues() const;

        void setValues(const std::vector<ProfileType> &values);

        // Properties
        double position[3];
        Bins bins;


    };

} /* namespace es */

#endif //SOURCE_PROFILE_H
