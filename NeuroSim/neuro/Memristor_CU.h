/*
 * Memristor_CU.h
 *
 *  Created on: 2019-8-26
 *      Author: jhxu
 */

#ifndef NEURO_MEMRISTOR_CU_HPP_
#define NEURO_MEMRISTOR_CU_HPP_
#include "../Array.h"
#include "../NeuroSim/SubArray.h"
#include "../Cell.h"
#include <math.h>
#include "../Param.h"
//#include "../Definition.h"
#include <stdint.h>


typedef enum{
    IMW, //input * weight
    WMI //weight * input
} OperationType;


template<class memoryType>
class Memristor_CU{
public:

	Memristor_CU(uint32_t M, uint32_t K, const float *weight, OperationType ot = WMI, short bits = 2, short row_size = 128, short col_size = 128, short cells = 12 ){
		this->ot = ot;
		this->M = M;
		this->K = K;
		tech = new Technology();
		inputParameter = new InputParameter();
		inputParameter->deviceRoadmap = HP;	// HP: high performance, LSTP: low power
		inputParameter->temperature = 301;	// Temperature (K)
		inputParameter->processNode = 14;	//processNode = 32 Technology node (nm)aasx
		tech->Initialize(inputParameter->processNode, inputParameter->deviceRoadmap);
		cell = new MemCell();
        this->row_size = row_size;
        this->col_size = col_size;
		maxExponentA = new unsigned char*[(int)(K/row_size)+1];
		for(int i = 0 ; i < (int)(K/row_size)+1; i++)
        {
		    maxExponentA[i] = new unsigned char[M];

        }
		this->A = weight;
		this->bits = bits;
		this->maxWeight = pow(2,bits) - 1;
		this->minWeight = 0;

		this->size = ceil((double) K / row_size) * ceil((double) M / (int)(col_size/(cells)));
        this->cells = cells;
		this->matrix_mapping(this->A);

	}

    void Energy(const float *input, uint64_t dims);
    void ReadLatency();




    double getWriteLatency(uint32_t row){
        return wLatency[row];
    }
    uint32_t getSize(){
        return size;
    }


    double getLatency()
    {
        return readLatency;
    }
    double getReadEnergy()
    {
        return readEnergy;
    }

	Memristor_CU(uint32_t M, uint32_t K, const double *weight, bool ot = WMI, short bits = 2, short row_size = 128, short col_size = 128, short cells = 12)
	{

	}

	void re_map(uint32_t MM, uint32_t KK, const float *weights, bool ot=WMI);


	void memristor_mm(const float *input, uint64_t dims, float *C, float alpha, float beta);


	void matrix_mapping(const float *weights);

	double getWriteEnergy(){
	    return writeEnergy;
	}


	~Memristor_CU(){
		delete tech;
		delete inputParameter;
		for(int i = 0 ; i < (int)(K/row_size)+1;i++)
        {
		    if(maxExponentA[i] != NULL)
		        delete[] maxExponentA[i];
		    maxExponentA[i] = NULL;
        }
		if(maxExponentA!=NULL)
		    delete[] maxExponentA; //TODO::
		if(cell!=NULL)
		    delete cell;
		if(signal_matrix!=NULL)
		    delete[] signal_matrix;
		if(currentSum !=NULL)
		    delete[] currentSum;
        if(operandsSumPositive !=NULL)
            delete[] operandsSumPositive;
        if(operandsSumNegative != NULL)
            delete[] operandsSumNegative;
		if(wLatency!=NULL)
		    delete[] wLatency;

		currentSum = NULL;
        operandsSumPositive = NULL;
        operandsSumNegative = NULL;
		tech = NULL;
		inputParameter = NULL;
		maxExponentA = NULL;
		cell = NULL;
		signal_matrix = NULL;




		for(int i = 0; i<(int)ceil((double)K / row_size) * (int)ceil((double)M / (col_size/cells)); i++)
		{
		    if(subArrays[i] != NULL){
		        delete subArrays[i];
		        subArrays[i] = NULL;
		    }
		    if(mp[i]!=NULL)
		    {
		        delete mp[i];
		        mp[i] = NULL;
		    }
		    if(mn[i]!=NULL)
		    {
		        delete mn[i];
		        mn[i] = NULL;
		    }
		}
		if(mn != NULL)
		    free(mn);
		if(mp!=NULL)
		    free(mp);
		mn = NULL;
		mp = NULL;


		if(subArrays!=NULL)
		    free(subArrays);
		if(adder!=NULL)
		    free(adder);
		if(mux!=NULL)
		    free(mux);
		if(muxDecoder!=NULL)
		    free(muxDecoder);
		if(dff!=NULL)
		    free(dff);
		if(subtractor!=NULL)
		    free(subtractor);

		subArrays = NULL;
		adder = NULL;
		mux = NULL;
		muxDecoder = NULL;
		dff = NULL;
		subtractor = NULL;

        wLatency = NULL;
		A = NULL;
	}


private:
	Array **mp= NULL,**mn= NULL;
	SubArray **subArrays= NULL;
	Adder *adder= NULL;
	Mux *mux= NULL;
	RowDecoder *muxDecoder= NULL;
	DFF *dff= NULL;
	Subtractor *subtractor= NULL;
	Technology *tech= NULL;
	MemCell *cell= NULL;
	InputParameter *inputParameter= NULL;
	short cells; //number of cells
	double readLatency = 0;
	double readEnergy = 0;
	double writeEnergy = 0;
	double writeLatency = 0;
	double *rLatency;
	double *wLatency;
	double heightNeuron;
	double widthNeuron;
    OperationType ot;
	const float *A = NULL;
	double *currentSum= NULL; //Used for calculating Energy

	double *operandsSumPositive= NULL; //Used for accelerating calculation
	double *operandsSumNegative= NULL; //Used for accelerating calculation

    //Matrix of M*K Multiply with Matrix of K*N
	uint32_t M; // M if ot=WMI or N if ot=IMW
	uint64_t K;
	unsigned char **maxExponentA= NULL;
	bool* signal_matrix = NULL;

    const int numColMuxed = 16;
    uint32_t writeDriverMuxed = 16;

	uint32_t size;


	double maxWeight;
	double minWeight;
	short row_size;
	short col_size;
	short bits; //bits of each cell
};



#endif /* NEURO_MEMRISTOR_CU_HPP_ */
