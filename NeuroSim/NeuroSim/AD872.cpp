#include "AD872.h"
#include <cmath>
#include <iostream>
#include "constant.h"
#include "formula.h"
using namespace std;

AD872::AD872(const InputParameter& _inputParameter, const Technology& _tech, const MemCell& _cell): ReadCircuit(_inputParameter, _tech, _cell) {
	initialized = false;
}
void AD872::CalculateLatency(double numRead) {
	if (!initialized) {
		cout << "[AD872] Error: Require initialization first!" << endl;
	} else {
		readLatency = 0;

		readLatency += 1.2/(1e9); //AD872 is 10M Hz frequency.
		readLatency += 1/clkFreq;
		//printf("numRead = %.4e\n",numRead);
		readLatency *= numRead;
        //printf("rdLatency=%.9e  readPulseWidth=%.9e   voltageIntThreshold=%.9e  maxIntPer=%.9e  \n",readLatency,cell.readPulseWidth,voltageIntThreshold,maxNumIntPerCycle);
	}
}
void AD872::CalculatePower(double numRead) {
	if (!initialized) {
		cout << "[ReadCircuit] Error: Require initialization first!" << endl;
	} else {
		leakage = 0;
		readDynamicEnergy = 0;
		
		// Leakage (DFF) (rough calculation)
		leakage += CalculateGateLeakage(INV, 1, widthDffInvN, widthDffInvP, inputParameter.temperature, tech) * tech.vdd * 8 * numDff;
		// Leakage (Read circuit body)
		if (mode == CMOS) {
			// Analytical result
			leakage += CalculateGateLeakage(INV, 1, widthNmos1, widthPmos1, inputParameter.temperature, tech) * tech.vdd;
			leakage += CalculateGateLeakage(INV, 1, widthNmos2, 0, inputParameter.temperature, tech) * tech.vdd;
			leakage += CalculateGateLeakage(INV, 1, widthNmos3, widthPmos3, inputParameter.temperature, tech) * tech.vdd;
			leakage += CalculateGateLeakage(INV, 1, widthNmos4, widthPmos4, inputParameter.temperature, tech) * tech.vdd;
			leakage += CalculateGateLeakage(INV, 1, widthNmos5, widthPmos5, inputParameter.temperature, tech) * tech.vdd;
			leakage += CalculateGateLeakage(INV, 1, widthNmos6, 0, inputParameter.temperature, tech) * tech.vdd;
			leakage += CalculateGateLeakage(INV, 1, widthNmos7, 0, inputParameter.temperature, tech) * tech.vdd;
			leakage += CalculateGateLeakage(INV, 1, widthNmos8, widthPmos8, inputParameter.temperature, tech) * tech.vdd;
			// Buffer
			leakage += CalculateGateLeakage(INV, 1, widthInvN, widthInvP, inputParameter.temperature, tech) * tech.vdd * 2;

			// SPICE result	(65nm tech node)
			//leakage += 104.9e-6;

		} else {	// mode==OSCILLATION, only one INV
			// Analytical result
			//leakage += CalculateGateLeakage(INV, 1, widthInvN, widthInvP, inputParameter.temperature, tech) * tech.vdd;

			// SPICE result (65nm tech node)
			leakage += 35.84e-9;
		}

		leakage *= numReadCol;

        readDynamicEnergy = 312.5e-6*1e-7;//result from SPICE
		readDynamicEnergy *= numReadCol;
		readDynamicEnergy *= numRead;

		if (!readLatency) {
			//cout << "[ReadCircuit] Error: Need to calculate read latency first" << endl;
		} else {
			readPower = readDynamicEnergy/readLatency;
		}
	}
}
