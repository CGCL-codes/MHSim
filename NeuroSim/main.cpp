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

#include <cstdio>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <random>
#include <vector>
#include "Cell.h"
#include "Array.h"
#include "formula.h"
#include "NeuroSim.h"
//#include "Param.h"
#include "IO.h"
//#include "Train.h"
//#include "Test.h"
#include "Mapping.h"
//#include "Definition.h"
#include <math.h>
#include <sys/time.h>

#include "neuro/Memristor_CU.hpp"
//#include <mcheck.h>
extern "C" {
#include <cblas.h>
}


using namespace std;


int main() {

//	mtrace();
    std::mt19937 gen;
	gen.seed(0);








	//memristor_mm(unsigned int M, unsigned int N, unsigned int K,
	//				float alpha, float *A, float *B, float beta, float *C);
	printf("RealDevice size = %d\n",sizeof(RealDevice));
    //unsigned int M = 4096, N = 5, K = 9216;
	
    unsigned int M = 40, N = 5, K = 92;
    float *A = new float[M*K], *B = new float[K*N], *C = new float[M*N], *D = new float[M*N], *F = new float[M*N], *G = new float[M*N];


    for(int i = 0 ; i < M*K; i++)
	{
		A[i] = i;
    }
	for(int i = 0 ; i < K*N; i++)
	{
		B[i] = -2*i;
    }
	for(int i = 0 ; i < M*N; i++)
	{
		C[i] = 0.;
		D[i] = 0.;
		F[i] = 0.;
		G[i] = 0.;
	}

    struct timeval ss,s,e,mmmm;
    struct timezone tz;
    gettimeofday(&ss,&tz);

	clock_t sc,mc,ec;
    sc = clock();

    gettimeofday(&ss,&tz);
    Memristor_CU<RealDevice> *mcu = new Memristor_CU<RealDevice>(M,K,A,WMI,2,128,128,4);

    gettimeofday(&s,&tz);
    mcu->memristor_mm(B,N,F,1.,0.);
    gettimeofday(&mmmm,&tz);
    //mcu->ReadLatency();
    //mcu->Energy(B,N);


	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1., A, K, B,
	      N, 0., D, N);

    double writeLatency = 0.;
    for (int k = 0; k < K; ++k) {
        writeLatency += mcu->getWriteLatency(k);
    }
    printf("ReadLatency = %.5e   readEnergy = %.5e    writeLatency = %.5e  writeEnergy = %.5e\n",mcu->getLatency(),mcu->getReadEnergy(),writeLatency,mcu->getWriteEnergy());


//	delete param;
//	delete[] A;
//	delete[] B;
	delete[] D;
	delete[] C;
	delete[] F;
	delete[] G;

	delete mcu;
    	//Array *arrayIH1 = new Array(100, 400, 100);
    	//arrayIH1->GetMaxCellReadCurrent(0,0);

	/* Load in MNIST data */
//	ReadTrainingDataFromFile("patch60000_train.txt", "label60000_train.txt");
//	ReadTestingDataFromFile("patch10000_test.txt", "label10000_test.txt");
//
//	/* Initialization of synaptic array from input to hidden layer */
//	//arrayIH->Initialization<IdealDevice>();
//	arrayIH->Initialization<RealDevice>();
//	//arrayIH->Initialization<MeasuredDevice>();
//	//arrayIH->Initialization<SRAM>(param->numWeightBit);
//	//arrayIH->Initialization<DigitalNVM>(param->numWeightBit,true);
//
//
//	/* Initialization of synaptic array from hidden to output layer */
//	//arrayHO->Initialization<IdealDevice>();
//	arrayHO->Initialization<RealDevice>();
//	//arrayHO->Initialization<MeasuredDevice>();
//	//arrayHO->Initialization<SRAM>(param->numWeightBit);
//	//arrayHO->Initialization<DigitalNVM>(param->numWeightBit,true);
//
//
//	/* Initialization of NeuroSim synaptic cores */
//	param->relaxArrayCellWidth = 0;
//	NeuroSimSubArrayInitialize(subArrayIH, arrayIH, inputParameterIH, techIH, cellIH);
//	param->relaxArrayCellWidth = 1;
//	NeuroSimSubArrayInitialize(subArrayHO, arrayHO, inputParameterHO, techHO, cellHO);
//	/* Calculate synaptic core area */
//	NeuroSimSubArrayArea(subArrayIH);
//	NeuroSimSubArrayArea(subArrayHO);
//
//	/* Calculate synaptic core standby leakage power */
//	NeuroSimSubArrayLeakagePower(subArrayIH);
//	NeuroSimSubArrayLeakagePower(subArrayHO);
//
//	/* Initialize the neuron peripheries */
//	NeuroSimNeuronInitialize(subArrayIH, inputParameterIH, techIH, cellIH, adderIH, muxIH, muxDecoderIH, dffIH, subtractorIH);
//	NeuroSimNeuronInitialize(subArrayHO, inputParameterHO, techHO, cellHO, adderHO, muxHO, muxDecoderHO, dffHO, subtractorHO);
//	/* Calculate the area and standby leakage power of neuron peripheries below subArrayIH */
//	double heightNeuronIH, widthNeuronIH;
//	NeuroSimNeuronArea(subArrayIH, adderIH, muxIH, muxDecoderIH, dffIH, subtractorIH, &heightNeuronIH, &widthNeuronIH);
//	double leakageNeuronIH = NeuroSimNeuronLeakagePower(subArrayIH, adderIH, muxIH, muxDecoderIH, dffIH, subtractorIH);
//	/* Calculate the area and standby leakage power of neuron peripheries below subArrayHO */
//	double heightNeuronHO, widthNeuronHO;
//	NeuroSimNeuronArea(subArrayHO, adderHO, muxHO, muxDecoderHO, dffHO, subtractorHO, &heightNeuronHO, &widthNeuronHO);
//	double leakageNeuronHO = NeuroSimNeuronLeakagePower(subArrayHO, adderHO, muxHO, muxDecoderHO, dffHO, subtractorHO);
//
//	/* Print the area of synaptic core and neuron peripheries */
//	double totalSubArrayArea = subArrayIH->usedArea + subArrayHO->usedArea;
//	double totalNeuronAreaIH = adderIH.area + muxIH.area + muxDecoderIH.area + dffIH.area + subtractorIH.area;
//	double totalNeuronAreaHO = adderHO.area + muxHO.area + muxDecoderHO.area + dffHO.area + subtractorHO.area;
//	printf("Total SubArray (synaptic core) area=%.4e m^2\n", totalSubArrayArea);
//	printf("Total Neuron (neuron peripheries) area=%.4e m^2\n", totalNeuronAreaIH + totalNeuronAreaHO);
//	printf("Total area=%.4e m^2\n", totalSubArrayArea + totalNeuronAreaIH + totalNeuronAreaHO);
//
//	/* Print the standby leakage power of synaptic core and neuron peripheries */
//	printf("Leakage power of subArrayIH is : %.4e W\n", subArrayIH->leakage);
//	printf("Leakage power of subArrayHO is : %.4e W\n", subArrayHO->leakage);
//	printf("Leakage power of NeuronIH is : %.4e W\n", leakageNeuronIH);
//	printf("Leakage power of NeuronHO is : %.4e W\n", leakageNeuronHO);
//	printf("Total leakage power of subArray is : %.4e W\n", subArrayIH->leakage + subArrayHO->leakage);
//	printf("Total leakage power of Neuron is : %.4e W\n", leakageNeuronIH + leakageNeuronHO);
//
//	/* Initialize weights and map weights to conductances for hardware implementation */
//	WeightInitialize();
//	if (param->useHardwareInTraining) { WeightToConductance(); }
//
//	srand(0);	// Pseudorandom number seed
//
//	ofstream mywriteoutfile;
//	mywriteoutfile.open("my_log.csv");
//
//	for (int i=1; i<=param->totalNumEpochs/param->interNumEpochs; i++) {
//        //cout << "Training Epoch : " << i << endl;
//		Train(param->numTrainImagesPerEpoch, param->interNumEpochs,param->optimization_type);
//		if (!param->useHardwareInTraining && param->useHardwareInTestingFF) { WeightToConductance(); }
//		Validate();
//		mywriteoutfile << i*param->interNumEpochs << ", " << (double)correct/param->numMnistTestImages*100 << endl;
//
//		printf("Accuracy at %d epochs is : %.2f%\n", i*param->interNumEpochs, (double)correct/param->numMnistTestImages*100);
//		printf("\tRead latency=%.4e s\n", subArrayIH->readLatency + subArrayHO->readLatency);
//		printf("\tWrite latency=%.4e s\n", subArrayIH->writeLatency + subArrayHO->writeLatency);
//		printf("\tRead energy=%.4e J\n", arrayIH->readEnergy + subArrayIH->readDynamicEnergy + arrayHO->readEnergy + subArrayHO->readDynamicEnergy);
//		printf("\tWrite energy=%.4e J\n", arrayIH->writeEnergy + subArrayIH->writeDynamicEnergy + arrayHO->writeEnergy + subArrayHO->writeDynamicEnergy);
//	}
//	printf("\n");
	return 0;
}


