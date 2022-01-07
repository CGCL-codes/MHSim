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

#include "GEMM/Memristor_CU.hpp"
//#include <mcheck.h>
extern "C" {
#include <cblas.h>
}


using namespace std;


int main() {
//	mtrace();
    std::mt19937 gen;
	gen.seed(0);
    unsigned int M = 20, N = 159, K = 25;
    //unsigned int M = 4, N = 4, K = 3;
    float *A = new float[M*K], *B = new float[K*N], *C = new float[M*N], *D = new float[M*N], *E = new float[M*N], *F = new float[M*N], *G = new float[M*N];
    char comma;
    for(int i = 0 ; i < M*K; i++)
	{
        A[i] = (i+1)/K;
    }
	for(int i = 0 ; i < K*N; i++)
	{
        B[i] = i/K;
    }
	for(int i = 0 ; i < M*N; i++)
	{
		C[i] = 0.;
		D[i] = 0.;
		F[i] = 0.;
        E[i] = 0.;
		G[i] = 0.;
	}
    float *A1 = new float[M*K], *B1 = new float[K*N], *B2 = new float[K*N];
    for(int i = 0; i < K; i++)
        for(int j = 0; j < N; j++)
            B2[i*N + j] = B[i + j*K];
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1., A, K, B2,
	      N, 0., C, N);
    Memristor_CU<RealDevice> *mcu = new Memristor_CU<RealDevice>(M,K,A,WMI,2,128,128,8);
    Memristor_CU<RealDevice> *mcu1 = new Memristor_CU<RealDevice>(N,K,B2,IMW,2,128,128,8);
    mcu->memristor_mm(B2,N,F,1.,0.,16);
    mcu1->memristor_mm(A,M,F,1.,0.,16);
    mcu->memristor_mm(B,N,G,1.,0.,16);
    mcu->ReadLatency(16);
	mcu->Energy(B,N);
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, 1., A, K, B,
	      K, 0., D, N);
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1., A, K, B,
	      N, 0., E, N);

    for(int ii = 0 ; ii < 9; ii++){
        printf("F = %.4e  D = %.4e  C = %.4e E = %.4e  G = %.4e \n", F[ii], D[ii], C[ii], E[ii], G[ii]);
    }

    double writeLatency = 0.;
    for (int k = 0; k < K; ++k) {
        writeLatency += mcu->getWriteLatency(k);
    }
    printf("ReadLatency = %.5e   readEnergy = %.5e    writeLatency = %.5e  writeEnergy = %.5e\n",mcu->getLatency(),mcu->getReadEnergy(),writeLatency,mcu->getWriteEnergy());


	delete[] A;
	delete[] B;
    delete[] A1;
    delete[] B1;
    delete[] B2;
	delete[] D;
	delete[] C;
	delete[] F;
    delete[] E;
	delete[] G;

	delete mcu;
	return 0;
}


