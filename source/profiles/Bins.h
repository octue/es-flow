/*
 * Bins.h
 *
 *  Created on: 29 Nov 2016
 *      Author: thc29
 */

#ifndef SOURCE_BASIC_PROFILES_BINS_H_
#define SOURCE_BASIC_PROFILES_BINS_H_

namespace es {

class Bin {
public:
	double x;
	double y;
	double z;
	double dx;
	double dy;
	double dz;
	
	Bin();
	virtual ~Bin();
};

} /* namespace es */

#endif /* SOURCE_BASIC_PROFILES_BINS_H_ */
