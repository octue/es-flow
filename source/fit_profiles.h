/*
 * Fit.h
 *
 *  Created on: 29 Nov 2016
 *      Author: thc29
 */

#ifndef SOURCE_FIT_H_
#define SOURCE_FIT_H_

namespace es {

class Fit {
public:
	Fit();
	virtual ~Fit();
	int Fitness() {
		return 1;
	};
};

} /* namespace es */

#endif /* SOURCE_FIT_H_ */
