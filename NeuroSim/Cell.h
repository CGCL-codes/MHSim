/*******************************************************************************
* Copyright (c) 2015-2017
* School of Electrical, Computer and Energy Engineering, Arizona State University
* PI: Prof. Shimeng Yu
* All rights reserved.
*   
* This source code is part of NeuroSim - a device-circuit-algorithm framework to benchmark 
* neuro-inspired architectures with synaptic devices(e.g., SRAM and emerging non-volatile memory). 
* Copyright of the model is maintained by the developers, and the model is distributed under 
* the terms of the Creative Commons Attribution-NonCommercial 4.0 International Public License 
* http://creativecommons.org/licenses/by-nc/4.0/legalcode.
* The source code is free and you can redistribute and/or modify it
* by providing that the following conditions are met:
*   
*  1) Redistributions of source code must retain the above copyright notice,
*     this list of conditions and the following disclaimer. 
*   
*  2) Redistributions in binary form must reproduce the above copyright notice,
*     this list of conditions and the following disclaimer in the documentation
*     and/or other materials provided with the distribution.
*   
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
* SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
* 
* Developer list: 
*   Pai-Yu Chen     Email: pchen72 at asu dot edu 
*                     
*   Xiaochen Peng   Email: xpeng15 at asu dot edu
********************************************************************************/

#ifndef CELL_H_
#define CELL_H_

#include <random>
#include <vector>

class Cell {
public:
	int x, y;	// Cell location: x (column) and y (row) start from index 0
	static double heightInFeatureSize, widthInFeatureSize;	// Cell height/width in terms of feature size (F)
	static double area;	// Cell area (m^2)
	virtual ~Cell() {}	// Add a virtual function to enable dynamic_cast
};

class eNVM: public Cell {
public:
	static double readVoltage;	// On-chip read voltage (Vr) (V)
	static double readPulseWidth;	// Read pulse width (s) (will be determined by ADC)
	static double readEnergy;	// Dynamic variable for calculation of read energy (J)
	static double writeVoltageLTP;	// Write voltage (V) for LTP or weight increase
	static double writeVoltageLTD;	// Write voltage (V) for LTD or weight decrease
	static double writePulseWidthLTP;	// Write pulse width (s) of LTP or weight increase
	static double writePulseWidthLTD;	// Write pulse width (s) of LTD or weight decrease
	static double writeEnergy;	// Dynamic variable for calculation of write energy (J)
	double conductance;	// Current conductance (S) (Dynamic variable) at on-chip Vr (different than the Vr in the reported measurement data)
	double conductancePrev;	// Previous conductance (S) (Dynamic variable) at on-chip Vr (different than the Vr in the reported measurement data)
	static double maxConductance;	// Maximum cell conductance (S)
	static double minConductance;	// Minimum cell conductance (S)
	static double avgMaxConductance;   // Average maximum cell conductance (S)
	static double avgMinConductance;   // Average minimum cell conductance (S)
	static bool cmosAccess;	// True: Pseudo-crossbar (1T1R), false: cross-point
    static bool isSTTMRAM; // if it is a STTMRAM device
	static bool FeFET;			// True: FeFET structure (Pseudo-crossbar only, should be cmosAccess=1)
	static double resistanceAccess;	// The resistance of transistor (Ohm) in Pseudo-crossbar array when turned ON
	static bool nonlinearIV;	// Consider I-V nonlinearity or not (Currently this option is for cross-point array. It is hard to have this option in pseudo-crossbar since it has an access transistor and the transistor's resistance can be comparable to RRAM's resistance after considering the nonlinearity. In this case, we have to iteratively find both the resistance and Vw across RRAM.)
	static bool readNoise;	// Consider read noise or not
	static double sigmaReadNoise;	// Sigma of read noise in gaussian distribution
	static double NL;	// Nonlinearity in write scheme (the current ratio between Vw and Vw/2), assuming for the LTP side


	static double conductanceAtVwLTP;		// Conductance at the LTP write voltage
	static double conductanceAtVwLTD;		// Conductance at the LTD write voltage
	static double conductanceAtHalfVwLTP;	// Conductance at 1/2 LTP write voltage
	static double conductanceAtHalfVwLTD;	// Conductance at 1/2 LTD write voltage
	static bool conductanceRangeVar;	// Consider variation of conductance range or not
	static double maxConductanceVar;	// Sigma of maxConductance variation (S)
	static double minConductanceVar;	// Sigma of minConductance variation (S)


    static std::mt19937 localGen;
    static bool init;


	//std::normal_distribution<double> *gaussian_dist4;	// Normal distribution object
	//std::normal_distribution<double> *gaussian_dist5;	// Normal distribution object
	virtual ~eNVM(){
//		printf("111 \n");
//		if(gaussian_dist != NULL)
//		{
//			delete gaussian_dist;
//			gaussian_dist = NULL;
//		}
//		if(gaussian_dist2 != NULL)
//		{
//			delete gaussian_dist2;
//			gaussian_dist2 = NULL;
//		}
//		if(gaussian_dist3 != NULL)
//		{
//			delete gaussian_dist3;
//			gaussian_dist3 = NULL;
//		}
//		if(gaussian_dist4 != NULL)
//		{
//			delete gaussian_dist4;
//			gaussian_dist4 = NULL;
//		}
//		if(gaussian_dist5 != NULL)
//		{
//			delete gaussian_dist5;
//			gaussian_dist5 = NULL;
//		}
//		if(gaussian_dist_maxConductance != NULL)
//		{
//			delete gaussian_dist_maxConductance;
//			gaussian_dist_maxConductance = NULL;
//		}
//		if(gaussian_dist_minConductance != NULL)
//		{
//			delete gaussian_dist_minConductance;
//			gaussian_dist_minConductance = NULL;
//		}

	}
};


class SRAM: public Cell {
public:
	SRAM(int x, int y);
	int bit;	// Stored bit (1 or 0) (dynamic variable)
	int bitPrev;	// Previous bit
	double widthSRAMCellNMOS;	// Pull-down NMOS width in terms offeature size (F)
	double widthSRAMCellPMOS;	// Pull-up PMOS width in terms of feature size (F)
	double widthAccessCMOS;		// Access transistor width in terms of feature size (F)
	double minSenseVoltage;		// Minimum voltage difference (V) for sensing
	double readEnergy;			// Dynamic variable for calculation of read energy (J) 
	double writeEnergy;			// Dynamic variable for calculation of write energy (J)
	double readEnergySRAMCell;	// Read energy (J) per SRAM cell (currently not used, it is included in the peripheral circuits of SRAM array in NeuroSim)
	double writeEnergySRAMCell;	// Write energy (J) per SRAM cell (will be obtained from NeuroSim)
	double Read(){return 0;}	// Currently not used
	void Write(){}	// Currently not used
};

class AnalogNVM: public eNVM {
public:
	static int maxNumLevelLTP;	// Maximum number of conductance states during LTP or weight increase
	static int maxNumLevelLTD;	// Maximum number of conductance states during LTD or weight decrease
	static int numPulse;   // Number of write pulses used in the most recent write operation (Positive number: LTP, Negative number: LTD) (dynamic variable)
	static double writeLatencyLTP;	// Write latency of a cell during LTP or weight increase (different cells use different # write pulses, thus latency values are different). writeLatency will be calculated for each cell first, and then replaced by the maximum one in the batch write.
	static double writeLatencyLTD;	// Write latency of a cell during LTD or weight decrease (different cells use different # write pulses, thus latency values are different). writeLatency will be calculated for each cell first, and then replaced by the maximum one in the batch write.
	static bool FeFET;			// True: FeFET structure (Pseudo-crossbar only, should be cmosAccess=1)
	static double gateCapFeFET;	// Gate Capacitance of FeFET (F)
	/* Non-identical write pulse scheme */
	static bool nonIdenticalPulse;	// Use non-identical pulse scheme in weight update or not (put the parameter here due to the access from Train.cpp)
	static double VinitLTP;    // Initial write voltage for LTP or weight increase (V)
	static double VstepLTP;    // Write voltage step for LTP or weight increase (V)
	static double VinitLTD;    // Initial write voltage for LTD or weight decrease (V)
	static double VstepLTD;    // Write voltage step for LTD or weight decrease (V)
	static double PWinitLTP;   // Initial write pulse width for LTP or weight increase (s)
	static double PWstepLTP;   // Write pulse width for LTP or weight increase (s)
	static double PWinitLTD;   // Initial write pulse width for LTD or weight decrease (s)
	static double PWstepLTD;   // Write pulse width for LTD or weight decrease (s)
	static double writeVoltageSquareSum;   // Sum of V^2 of non-identical pulses (for weight update energy calculation in subcircuits)

	virtual double Read(double voltage) = 0;
	virtual void Write(double deltaWeightNormalized, double weight, double minWeight, double maxWeight) = 0;
	double GetMaxReadCurrent(){
      if(cmosAccess)
          return readVoltage * 1/(1/avgMaxConductance+resistanceAccess);
      else
          return readVoltage * avgMaxConductance;}
	double GetMinReadCurrent(){
      if(cmosAccess)
          return readVoltage * 1/(1/avgMinConductance+resistanceAccess);
      else
          return readVoltage * avgMinConductance;}
	void WriteEnergyCalculation(double wireCapCol);
};

class DigitalNVM: public eNVM {
public:
    static std::normal_distribution<double> *gaussian_dist;	// Normal distribution object
    static std::normal_distribution<double> *gaussian_dist_maxConductance;	// Normal distribution object
    static bool gd_init;
    static std::normal_distribution<double> *gaussian_dist_minConductance;	// Normal distribution object
    DigitalNVM(int x, int y);
	int bit;	// Stored bit (1 or 0) (dynamic variable), for internel check only and not be used for read
	int bitPrev;	// Previous bit
	double refCurrent;	// Reference current for S/A
	double Read(double voltage);	// Return read current (A)
       // modified below
    bool isSTTMRAM;  // if it is STTMRAM, then, we can relax the cell area
    bool parallelRead; // if it is a parallel readout for STT-MRAM
	void Write(int bitNew, double wireCapCol);
};

class IdealDevice: public AnalogNVM {
public:
	IdealDevice(int x, int y);
	double Read(double voltage);	// Return read current (A)
	void Write(double deltaWeightNormalized, double weight, double minWeight, double maxWeight);
    static std::normal_distribution<double> *gaussian_dist;	// Normal distribution object
    static bool gd_init;
    static std::normal_distribution<double> *gaussian_dist_maxConductance;	// Normal distribution object
    static std::normal_distribution<double> *gaussian_dist_minConductance;	// Normal distribution object
    ~IdealDevice(){
//        if(gaussian_dist != NULL) {
//            delete gaussian_dist;
//        }
//        if(gaussian_dist_maxConductance!=NULL)
//            delete gaussian_dist_maxConductance;
//        if(gaussian_dist_minConductance!=NULL)
//            delete gaussian_dist_minConductance;
//        gaussian_dist = NULL;
//        gaussian_dist_maxConductance = NULL;
//        gaussian_dist_minConductance = NULL;
	}
};

class RealDevice: public AnalogNVM {
public:
	bool nonlinearWrite;	// Consider weight update nonlinearity or not
	static double xPulse;		// Conductance state in terms of the pulse number (doesn't need to be integer)
	static double NL_LTP;		// LTP nonlinearity
	static double NL_LTD;		// LTD nonlinearity
	static double paramALTP;	// Parameter A for LTP nonlinearity
	static double paramBLTP;	// Parameter B for LTP nonlinearity
	static double paramALTD;	// Parameter A for LTD nonlinearity
	static double paramBLTD;	// Parameter B for LTD nonlinearity
	static double sigmaDtoD;	// Sigma of device-to-device variation on weight update nonliearity baseline
	static double sigmaCtoC;	// Sigma of cycle-to-cycle variation on weight update
    static std::normal_distribution<double> *gaussian_dist;	// Normal distribution object
    static std::normal_distribution<double> *gaussian_dist2;	// Normal distribution object
    static std::normal_distribution<double> *gaussian_dist3;	// Normal distribution object
    static std::normal_distribution<double> *gaussian_dist_maxConductance;	// Normal distribution object
    static std::normal_distribution<double> *gaussian_dist_minConductance;	// Normal distribution object
    static bool gd_init;
	RealDevice(int x, int y);
	double Read(double voltage);	// Return read current (A)
	void Write(double deltaWeightNormalized, double weight, double minWeight, double maxWeight);
	~RealDevice(){
//        if(gaussian_dist != NULL) {
//            delete gaussian_dist;
//        }
//        if(gaussian_dist_maxConductance!=NULL)
//            delete gaussian_dist_maxConductance;
//        if(gaussian_dist_minConductance!=NULL)
//            delete gaussian_dist_minConductance;
//        if(gaussian_dist3!=NULL)
//            delete gaussian_dist3;
//        if(gaussian_dist2!=NULL)
//            delete gaussian_dist2;
//        gaussian_dist = NULL;
//        gaussian_dist_maxConductance = NULL;
//        gaussian_dist_minConductance = NULL;
//        gaussian_dist3 = NULL;
//        gaussian_dist2 = NULL;
	}
};


class MeasuredDevice: public AnalogNVM {
public:
	bool nonlinearWrite;	// Consider weight update nonlinearity or not
	bool symLTPandLTD;	// True: use LTP conductance data for LTD
	double xPulse;		// Conductance state in terms of the pulse number (doesn't need to be integer)
	std::vector<double> dataConductanceLTP;	// LTP conductance data at different pulse number
	std::vector<double> dataConductanceLTD;	// LTD conductance data at different pulse number
    std::normal_distribution<double> *gaussian_dist;
	MeasuredDevice(int x, int y);
	double Read(double voltage);	// Return read current (A)
	void Write(double deltaWeightNormalized, double weight, double minWeight, double maxWeight);
	~MeasuredDevice(){
//		printf("222");
//		printf("release gaussian1! \n");
		if(gaussian_dist != NULL)
		{
//			printf("release gaussian2! \n");
			delete gaussian_dist;
			gaussian_dist = NULL;
		}
		std::vector<double> x;
		dataConductanceLTD.swap(x);
		dataConductanceLTP.swap(x);

	}
};

#endif
