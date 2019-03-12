/*
 * biot_savart.h A parallelised biot-savart function for use with the Attached-Detached Eddy Method
 *
 * Author:              Tom Clark  (thclark @ github)
 *
 * Copyright (c) 2014-2019 Octue Ltd. All Rights Reserved.
 *
 */

#ifndef ES_FLOW_BIOT_SAVART_H
#define ES_FLOW_BIOT_SAVART_H

#include <tbb/tbb.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/atomic.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <iostream>
#include <Eigen/Dense>


using namespace Eigen;
using namespace std;


#define pi          3.14159265358979323846264;
#define fourpi      12.5663706143592172953851;
#define invfourpi	0.079577471545948;


namespace es {


// @cond - Don't document these convenience inline functions or the kernel

// Inverse magnitude of a Vector3d element
static inline double InvMagnitude(Vector3d vec){
    return 1 / vec.norm();
}


// Inverse magnitude of a Vector3d element handling /0 errors to return 0 always
static inline double InvMagnitudeZero(Vector3d vec)
{
    double a = vec.norm();
    if (a == 0.) {
        return 0.;
    } else {
        return 1 / a;
    }
}


// Kernel that applies induction from one vortex line to one control point
static inline Vector3d BiotElementalKernel(Vector3d control_locations, Vector3d start_nodes, Vector3d end_nodes, double gamma, double effective_core_radius_squared)
{

    Vector3d l12 = end_nodes - start_nodes;

    Vector3d r1 = control_locations - start_nodes;

    Vector3d r2 = control_locations - end_nodes;

    Vector3d cross1 = l12.cross(r1);

    // Used invmagnitudezero to solve the problem where the cross product is {0,0,0} due to point p being on the line defined by l12 (r1 collinear with l12)
    Vector3d induction;
    induction = InvMagnitudeZero(cross1) * cross1;
    double mag_l12 = InvMagnitude(l12);
    double mag_r1 = InvMagnitudeZero(r1);
    double mag_r2 = InvMagnitudeZero(r2);
    double costheta_1 = l12.dot(r1) * (mag_l12 * mag_r1);
    double costheta_2 = l12.dot(r2) * (mag_l12 * mag_r2);
    double h = (cross1.norm()) * (mag_l12);
    double magu = gamma * h * (costheta_1 - costheta_2) * invfourpi;
    magu /= sqrt(pow(effective_core_radius_squared, 2) + pow(h, 4));
    induction = induction * magu;

    return induction;

};


// @endcond


/** @brief O(N^2) 'Naive' application of the biot savart law to determine the action of many vortex lines on many control points
 *
 * @param[in] startNodes [3 x n] Start nodes of vortex lines in x, y, z coords
 * @param[in] endNodes [3 x n] End nodes of vortex lines in x, y, z coords
 * @param[in] locations [3 x p] Locations in x, y, z coords at which induction is required
 * @param[in] gamma [n x 1] Circulation of each vortex line
 * @param[in] rcEffSqd [n x 1] Square of the effective vortex core radius for each vortex line
 *
 *  Outputs:
 *
 *    induction[3 x p]     double
 *    Induction in u, v, w directions of the input vortex lines at the input control points
 *
 */
Matrix3Xd NaiveBiotSavart(Matrix3Xd startNodes, Matrix3Xd endNodes, Matrix3Xd locations, VectorXd gamma, VectorXd effective_core_radius_squared){

    assert (startNodes.cols() == gamma.size());
    assert (startNodes.cols() == endNodes.cols());
    assert (startNodes.cols() == effective_core_radius_squared.size());
    Index start = 0;
    Index n_control_point = locations.cols();
    Index n_vortex = gamma.size();
    Index n = n_control_point * n_vortex;

    Matrix3Xd induction(3, n_control_point);
    induction.setZero();

    // Parallel version
    tbb::affinity_partitioner partitioner1;
    tbb::affinity_partitioner partitioner2;

    tbb::parallel_for(start, n_control_point, [&](Index control_i){

        Vector3d p = locations.col(control_i);

        Vector3d induction_accumulator;
        induction_accumulator.setZero();
        tbb::spin_mutex induction_mutex;

        tbb::parallel_for(start, n_vortex, [&](Index vortex_i){
            Vector3d A1 = startNodes.col(vortex_i);
            Vector3d B1 = endNodes.col(vortex_i);
            double a = gamma(vortex_i);
            double rc4 = effective_core_radius_squared(vortex_i);
            Vector3d induction_elemental = BiotElementalKernel(p, A1, B1, a, rc4);
            {
                tbb::spin_mutex::scoped_lock lock(induction_mutex);
                induction_accumulator += induction_elemental;
            }
        }, partitioner1);
        induction.block(0, control_i, 3, 1) = induction_accumulator;

    }, partitioner2);

    /* Serial version
    for (long control_i = 0; control_i < n_control_point; control_i++)
    {
        Vector3d p = locations.block(0,control_i,3,1);
        long vortex_i;
        for (vortex_i = 0; vortex_i < n_vortex; vortex_i++) {
            Vector3d A1 = startNodes.block(0, vortex_i, 3, 1);
            Vector3d B1 = endNodes.block(0, vortex_i, 3, 1);
            double a = gamma(vortex_i);
            double rc4 = effective_core_radius_squared(vortex_i);
            induction.block(0, control_i, 3, 1) = induction.block(0, control_i, 3, 1) + BiotElementalKernel(p, A1, B1, a, rc4);
        };

    }
    */

    return induction;
};

};

#endif
