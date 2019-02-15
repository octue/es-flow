/*
 * definitions.h Constants, definitions and identities
 *
 * Author:              Tom Clark  (thclark @ github)
 *
 * Copyright (c) 2019 Octue Ltd. All Rights Reserved.
 *
 */

#ifndef ES_FLOW_DEFINITIONS_H
#define ES_FLOW_DEFINITIONS_H

// Rotational speed of the world in rad/s
#define OMEGA_WORLD 7.2921159e-05

// von Karman constant
#define KAPPA_VON_KARMAN 0.41

// Degrees based trigonometry
#define sind(x) (sin(fmod((x), 360.0) * M_PI / 180.0))
#define cosd(x) (cos(fmod((x), 360.0) * M_PI / 180.0))
#define asind(x) (asin(x) * 180.0 / M_PI)
#define acosd(x) (acos(x) * 180.0 / M_PI)
#define tand(x) (tan(x * M_PI / 180.0) * M_PI / 180.0)

#endif // ES_FLOW_DEFINITIONS_H
