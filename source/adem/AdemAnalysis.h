/*
 * AdemAnalysis.h
 *
 *  Created on: 29 Nov 2016
 *      Author: thc29
 */

#ifndef SOURCE_ADEM_ADEMANALYSIS_H_
#define SOURCE_ADEM_ADEMANALYSIS_H_

namespace es {

class AdemAnalysis {
public:
	AdemAnalysis();
	virtual ~AdemAnalysis();
};

class WindAdemAnalysis: public AdemAnalysis {
public:
	WindAdemAnalysis();
	virtual ~WindAdemAnalysis();
};

class TidalAdemAnalysis: public AdemAnalysis {
public:
	TidalAdemAnalysis();
	virtual ~TidalAdemAnalysis();
};

} /* namespace es */

#endif /* SOURCE_ADEM_ADEMANALYSIS_H_ */
