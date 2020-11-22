#ifndef AD872_H_
#define AD872_H_

#include "ReadCircuit.h"

class AD872 : public ReadCircuit{
public:

	AD872(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell);

    double rdPulseWidth;
	void CalculateLatency(double numRead);
	void CalculatePower(double numRead);
};


#endif
