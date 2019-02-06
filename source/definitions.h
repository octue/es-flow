/*
 * definitions.h Constants, definitions and identities
 *
 * Author:              Tom Clark  (thclark @ github)
 *
 * Copyright (c) 2019 Octue Ltd. All Rights Reserved.
 *
 */

#ifndef SOURCE_DEFINITIONS_H
#define SOURCE_DEFINITIONS_H

// Rotational speed of the world in rad/s
double omega_world = 7.2921159e-05;

// Degrees based trigonometry
#define sind(x) (sin(fmod((x), 360.0) * M_PI / 180.0))
#define cosd(x) (cos(fmod((x), 360.0) * M_PI / 180.0))
#define asind(x) (asin(x) * 180.0 / M_PI)
#define acosd(x) (acos(x) * 180.0 / M_PI)
#define tand(x) (tan(x * M_PI / 180.0) * M_PI / 180.0)

#endif //SOURCE_DEFINITIONS_H
