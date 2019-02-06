/*
 * profile.h Profile handling for (e.g.) Velocity, Reynolds Stress and Spectral Tensor Profile management
 *
 * Function definitions are included here, rather than the cpp file, because of the templating constraint. See http://stackoverflow.com/questions/495021/why-can-templates-only-be-implemented-in-the-header-file
 *
 * Author:                   Tom Clark (thclark @ github)
 *
 * Copyright (c) 2016-9 Octue Ltd. All Rights Reserved.
 *
 */

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
        Bins(std::vector<double> z);

        // Construct from spacing with bin sizes specified
        Bins(std::vector<double> z, std::vector<double> dx, std::vector<double> dy, std::vector<double> dz);

        // Construct from spacing and a half angle in degrees, dz set by central differencing of spacing
        Bins(std::vector<double> z, double half_angle_degrees);

        // Construct from spacing and a half angle in degrees, with dz specified
        Bins(std::vector<double> z, double half_angle_degrees, std::vector<double> dz);

        // Destroy
        ~Bins();

        // Vectors of z locations and bin dimensions
        unsigned long n_bins;
        std::vector<double> z_ctrs;
        std::vector<double> dx;
        std::vector<double> dy;
        std::vector<double> dz;

    };

    // Represent Bins class in logs or ostream
    ::std::ostream& operator<<(::std::ostream& os, const Bins& bins);

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
        const std::vector<ProfileType> &getValues() const;

        void setValues(const std::vector<ProfileType> &values);

        // Properties
        double position[3];
        Bins bins;

    };

    // Represent Profiles class in logs or ostream
    template <class ProfileType>
    ::std::ostream& operator<<(::std::ostream& os, const Profile<ProfileType>& profile);

    template <class ProfileType>
    Profile<ProfileType>::Profile(const Bins &bins) : bins(bins) {
        // Construct from just bins, with default position {0.0, 0.0, 0.0}
        position[0] = 0.0;
        position[1] = 0.0;
        position[2] = 0.0;
    }

    template <class ProfileType>
    Profile<ProfileType>::~Profile() {
        // Destructor
    }

    template <class ProfileType>
    Profile<ProfileType>::Profile(const Bins &bins, double x, double y, double z): bins(bins) {
        // Construct from bins and a position
        position[0] = x;
        position[1] = y;
        position[2] = z;
    }

    template <class ProfileType>
    void Profile<ProfileType>::setValues(const std::vector<ProfileType> &values) {
        if (values.size() != bins.n_bins) {
            throw std::out_of_range("size of vector 'values' does not equal the number of bins for this profile");
        }
        Profile::values = values;
    }

    template <class ProfileType>
    const std::vector<ProfileType> &Profile<ProfileType>::getValues() const {
        return values;
    }

//    template <class ProfileType>
//    const std::vector<ProfileType> &Profile<ProfileType>::getZDerivative() const {
//        // Get the first derivative with respect to Z location using Eigen's Autodiff
//        // Some examples can be found in eigen/unsupported/doc/examples/AutoDiff.cpp
//
//
//        return derivatives;
//    }

    template <class ProfileType>
    ::std::ostream& operator<<(::std::ostream& os, const Profile<ProfileType>& profile) {
        // Represent in logs or ostream
        return os << "debug statement for profile class";
    }



    class VelocityProfile: public Profile<double> {
    public:
        VelocityProfile();
        virtual ~VelocityProfile();

//        VectorXd AutoDiff();
    };

} /* namespace es */






#endif //SOURCE_PROFILE_H
