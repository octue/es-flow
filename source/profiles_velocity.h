/*
 * profiles_velocity.h
 *
 *  Created on: 29 Nov 2016
 *      Author: thc29
 */

#ifndef SOURCE_PROFILES_VELOCITY_H_
#define SOURCE_PROFILES_VELOCITY_H_

#include "profile.h"

namespace es {

	class VelocityProfile: public Profile {
	public:
		VelocityProfile();
		virtual ~StressProfile();
	};

} /* namespace es */

#endif /* SOURCE_PROFILES_VELOCITY_H_ */
