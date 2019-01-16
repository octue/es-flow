/*
 * constants.h
 *
 * References:
 *
 *   [1]
 *
 * Future Improvements:
 *
 *   [1] 
 *
 * Author:                   T. Clark
 * Work address:             Ocean Array Systems Ltd
 *                           Hauser Forum
 *                           3 Charles Babbage Road
 *                           Cambridge
 *                           CB3 0GT
 * Email:                    tom.clark@oceanarraysystems.com
 * Website:                  www.oceanarraysystems.com
 *
 * Copyright (c) 2017 Ocean Array Systems. All Rights Reserved.
 *
 */

#ifndef SOURCE_CONSTANTS_H
#define SOURCE_CONSTANTS_H

// Rotational speed of the world in rad/s
double omega_world = 7.2921159e-05;

// Degrees based trigonometry
#define sind(x) (sin(fmod((x),360) * M_PI / 180))
#define cosd(x) (cos(fmod((x),360) * M_PI / 180))
#define asind(x) (asin(x) * 180 / M_PI)
#define acosd(x) (acos(x) * 180 / M_PI)

#endif //SOURCE_CONSTANTS_H
