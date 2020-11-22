/*
 * Memristor_CU.cpp
 *
 *  Created on: 2019-8-26
 *      Author: jhxu
 */
#include "Memristor_CU.h"

#include "float_helper.h"
#include "../NeuroSim.h"
#include "../NeuroSim/SubArray.h"
#include "../Array.h"
#include "../Cell.h"
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <fstream>
extern "C" {
#include <cblas.h>
}
#include "omp.h"

template<class memoryType>
void Memristor_CU<memoryType>::matrix_mapping(const float *weights){
	const unsigned int *A1;


	currentSum = new double[M*K];

    operandsSumPositive = new double[M*K];
    operandsSumNegative = new double[M*K];
    for (int j = 0; j < M * K; ++j) {
        currentSum[j] = 0;
        operandsSumPositive[j] = 0;
        operandsSumNegative[j] = 0;
    }
    Array *tp = new Array(1,1<<bits,100);
    tp->Initialization<memoryType>();
    double *currentDictionary = new double[1<<bits];

    for(int i = 0; i < 1<<bits; i++){
        tp->WriteCell(0,i,-(maxWeight - minWeight) ,0,maxWeight,minWeight,true);
        tp->WriteCell(0,i,i,i,maxWeight,minWeight,true);
        currentDictionary[i] = tp->ReadCell(0,i);
    }

    delete tp;
    unsigned int *weight_mantissas = new unsigned int[M*K];

	if(ot == WMI)
	{
	    //dimensions of weight are M rows and K cols
		signal_matrix = new bool[K * M];
        wLatency = new double[K];
		int operands_per_row = (int)(col_size/(cells));
		unsigned short n = (int)ceil((double)K / row_size) * (int)ceil((double)M / operands_per_row);

		A1 = reinterpret_cast<const unsigned int*>(weights);
        mp = (Array**)malloc(sizeof(Array*)*n);
        mn = (Array**)malloc(sizeof(Array*)*n);

		subArrays = (SubArray**)malloc(sizeof(SubArray*)*n);
		adder = (Adder*)malloc(sizeof(Adder)*n);
		mux = (Mux*)malloc(sizeof(Mux)*n);
		muxDecoder = (RowDecoder*)malloc(sizeof(RowDecoder)*n);
		dff = (DFF*)malloc(sizeof(DFF)*n);
		subtractor = (Subtractor*)malloc(sizeof(Subtractor)*n);
		for (int i = 0; i < K * M; i++) {
			if(weights[i] == 0. || weights[i] == -0.)
			{
				weight_mantissas[i] = 0;
				continue;
			}
			if(get_exponent(A1[i]) == 255)
			{
				printf("weight = %f   Exponent = %d   Error: nan! \n", weights[i], get_exponent(A1[i]));
				exit(1);
			}
			weight_mantissas[i] = get_mantissas(A1[i]);
			signal_matrix[i] = 0;
		}
        for (int k = 0; k < K; ++k) {
            wLatency[k] = 0;
        }


//        Initialization
#pragma omp parallel for
        for (int i = 0; i < n; i++) {
            mp[i] = new Array(row_size, col_size, 100);
            mp[i]->Initialization<memoryType>();
            mn[i] = new Array(row_size, col_size, 100);
            mn[i]->Initialization<memoryType>();
            double heightNeuron, widthNeuron;

            new(adder + i) Adder(*inputParameter, *tech, *cell);
            new(mux + i) Mux(*inputParameter, *tech, *cell);
            new(dff + i) DFF(*inputParameter, *tech, *cell);
            new(subtractor + i) Subtractor(*inputParameter, *tech, *cell);
            new(muxDecoder + i) RowDecoder(*inputParameter, *tech, *cell);
//				subtractor[i] = Subtractor(inputParameter, tech, cell); //RealDevice use Subtractor?
            SubArray *subArray = NULL;
            NeuroSimSubArrayInitialize(subArray, mp[i], *inputParameter, *tech,
                                       *cell);
            NeuroSimSubArrayInitialize(subArray, mn[i], *inputParameter, *tech,
                                       *cell);

            NeuroSimSubArrayArea(subArray);
            NeuroSimSubArrayLeakagePower(subArray);
            NeuroSimNeuronInitialize(subArray, *inputParameter, *tech, *cell, adder[i],
                                     mux[i], muxDecoder[i], dff[i], subtractor[i]);
            NeuroSimNeuronArea(subArray, adder[i], mux[i], muxDecoder[i], dff[i], subtractor[i],
                               &heightNeuron, &widthNeuron);
            subArrays[i] = subArray;
        }


        double sumArrayWriteEnergy = 0;   // Use a temporary variable here since OpenMP does not support reduction on class member
        double sumNeuroSimWriteEnergy = 0;   // Use a temporary variable here since OpenMP does not support reduction on class member
        double sumWriteLatencyAnalogNVM = 0;   // Use a temporary variable here since OpenMP does not support reduction on class member

	   // maxExponentA = new unsigned char*[(int)(K/row_size)+1];
	   // for(int i = 0 ; i < (int)(K/row_size)+1; i++)
       // {
	   //     maxExponentA[i] = new unsigned char[M];

       // }
        unsigned char*** bbb = new unsigned char**[16];
        for(int i = 0 ; i < 16;i++){
            bbb[i] = new unsigned char*[(int)(K/row_size)+1];
            for(int j = 0 ; j < (int)(K/row_size+1);j++){
                bbb[i][j] = new unsigned char[M];
                for(int k = 0 ; k < M; k++)
                    bbb[i][j][k] = 0;
            }
        }

        int array_size = ceil((double)K / row_size); //Number of groups. Each group contains all the columns.
#pragma omp parallel for num_threads(16)
        for (int array_row = 0; array_row < array_size; array_row++) {
			for (int array_col = 0; array_col < ceil((double)M / operands_per_row); array_col++) { // array_col: Number of arrays that contains all the columns.
                for (int col = array_col * operands_per_row;
						col < M && col < operands_per_row * (array_col + 1); col++) { // col: The index of columns
					for (int row = array_row * row_size;
							row < K && row < row_size * (array_row + 1); row++) { // row: The index of columns
						//reset weight1
						//We modify the bound of weight from -1,1 to 0,1, since weight can only present positive value.
						for (int b = 0; b < (cells); b++) {
                            mp[(int)(array_row * ceil((double)M/operands_per_row))+array_col]->WriteCell(row - array_row * row_size, (col - array_col * operands_per_row) * (cells) + b,-(maxWeight - minWeight),
                                    0 /* delta_W=-(param->maxWeight-param->minWeight) will completely erase */,
                                    maxWeight, minWeight, true);
                            mn[(int)(array_row * ceil((double)M/operands_per_row))+array_col]->WriteCell(row - array_row * row_size, (col - array_col * operands_per_row) * (cells) + b,-(maxWeight - minWeight),
                                    0 /* delta_W=-(param->maxWeight-param->minWeight) will completely erase */,
                                    maxWeight, minWeight, true);

						}
					}
				}

				//Get the max exponent of weights that mapped in the same column.
				for (int row = 0; row < K; row++) {
					for (int col = 0; col < M; col++) {
						unsigned char exponentA = get_exponent(A1[col * K + row]); //Only map weight, the inputs are aligned in computation step, this may cause loss
                        int idx = omp_get_thread_num();
                        if (row - array_row * row_size == 0
						    		|| exponentA > bbb[idx][array_row][col])
                                bbb[idx][array_row][col] = exponentA;
							    //maxExponentA[array_row][col] = exponentA;
					}
				}

            }
        }

        for(int array_row = 0; array_row < array_size; array_row++){
            for(int col = 0; col < M; col++ ){
                unsigned char exponentA = 0;
                for(int i = 0; i < 16; i++){
                    if(bbb[i][array_row][col]>exponentA)
                        exponentA = bbb[i][array_row][col];
                }
                maxExponentA[array_row][col] = exponentA;
            }
        }

        //printf("end initial111\n");
#pragma omp parallel for
        for (int array_row = 0; array_row < array_size; array_row++) {
			for (int array_col = 0; array_col < ceil((double)M / operands_per_row); array_col++) {
                //Align the operands
				for (int row = array_row * row_size; row < K && row < row_size * (array_row + 1); row++) {
					for (int col = array_col * operands_per_row; col < M && col < operands_per_row * (array_col + 1); col++) {
						unsigned char exponentA = get_exponent(A1[col * K + row]);
						char diff = maxExponentA[array_row][col] - exponentA;
						//adjust the weight1
						weight_mantissas[col * K + row] = weight_mantissas[col * K + row]
								>> diff;
					}
				}
            }
        }

        for(int i = 0 ; i < 16; i++){
            for(int j = 0 ; j < (int)(K/row_size+1); j++)
                delete[] bbb[i][j];
            delete[] bbb[i];
        }

        delete[] bbb;
        //printf("end initial222\n");
#pragma omp parallel for reduction(+: sumArrayWriteEnergy, sumNeuroSimWriteEnergy, sumWriteLatencyAnalogNVM)
		for (int array_row = 0; array_row < array_size; array_row++) {
			for (int array_col = 0; array_col < ceil((double)M / operands_per_row); array_col++) { // array_col: Number of arrays that contains all the columns.



				for (int row = array_row * row_size;
                     row < K && row < row_size * (array_row + 1); row++) {

                    int numWriteOperationPerRow = 0; // # of write operations will be performed per row.
                    double numWriteOperation = 0;


                    double writeVoltageLTP = static_cast<eNVM*>(mp[0]->cell[0][0])->writeVoltageLTP;
                    double writeVoltageLTD = static_cast<eNVM*>(mp[0]->cell[0][0])->writeVoltageLTD;


					for (int col = array_col * operands_per_row; col < M && col < operands_per_row * (array_col + 1);
                         col++) {

					    bool weightChangeBatch = false;
                        double maxLatencyLTP = 0;
                        double maxLatencyLTD = 0;

                        double cSum = 0,opSumNegative = 0,opSumPositive = 0;


                        if (weights[col * K + row] < 0) //if weight1 is negative
                        {
                            signal_matrix[col * K + row] = 1; //Uselessa
                            for (int b = 0; b < (cells); b++) {  //Assume cells*colmux=colsize

                                //For WMI, we map the weight by transposed matrix
                                mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->WriteCell(
                                        row - array_row * row_size,
                                        (col - array_col * operands_per_row) * (cells) + b,
                                        (weight_mantissas[col * K + row]
                                                >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1),
                                        (weight_mantissas[col * K + row]
                                                >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1),
                                        maxWeight, minWeight,
                                        true);

                                //currentSum[col * K + row] means the sum of current when reading cells located in (col,row). (row,col) means the coordinate of cell in XB array.
                                //The operand at (col,row) is mapped into the cell with (row,col).
                                cSum += currentDictionary[(weight_mantissas[col * K + row]
                                                >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1)];
                                cSum += currentDictionary[0];

                                opSumNegative += currentDictionary[(weight_mantissas[col * K + row]
                                                >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1)] * (1<<(bits * ((cells-1) - b)));

                                opSumPositive +=  currentDictionary[0] * (1<<(bits * ((cells-1) - b)));



                                weightChangeBatch = weightChangeBatch || static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse;

                                if(static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTP > maxLatencyLTP)
                                    maxLatencyLTP = static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTP;
                                if(static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTD > maxLatencyLTD)
                                    maxLatencyLTD = static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTD;

                            }

#pragma omp critical
                            {
                                currentSum[col * K + row] += cSum;

                                operandsSumNegative[col * K + row] += opSumNegative;
                                operandsSumPositive[col * K + row] += opSumPositive;
                            }



                            numWriteOperationPerRow += weightChangeBatch;

                            for (int b = 0; b < cells; ++b) {
                                static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTP = maxLatencyLTP;
                                static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTD = maxLatencyLTD;
                                if(weightChangeBatch){
                                    if (static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->nonIdenticalPulse) {	// Non-identical write pulse scheme
                                        if (static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse > 0) {	// LTP
                                            static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTP = sqrt(static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageSquareSum / static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse);	// RMS value of LTP write voltage
                                            static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTD = static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VinitLTD + 0.5 * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VstepLTD * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->maxNumLevelLTD;	// Use average voltage of LTD write voltage
                                        } else if (static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse < 0) {	// LTD
                                            static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTP = static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VinitLTP + 0.5 * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VstepLTP * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->maxNumLevelLTP;    // Use average voltage of LTP write voltage
                                            static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTD = sqrt(static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageSquareSum / (-1*static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse));    // RMS value of LTD write voltage
                                        } else {	// Half-selected during LTP and LTD phases
                                            static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTP = static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VinitLTP + 0.5 * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VstepLTP * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->maxNumLevelLTP;    // Use average voltage of LTP write voltage
                                            static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTD = static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VinitLTD + 0.5 * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VstepLTD * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->maxNumLevelLTD;    // Use average voltage of LTD write voltage
                                        }
                                    }
                                }
                                static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->WriteEnergyCalculation(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapCol);
                                sumArrayWriteEnergy += static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeEnergy;
                            }

                            sumWriteLatencyAnalogNVM += maxLatencyLTP + maxLatencyLTD;
                            if(array_col == 0)
                            wLatency[row] += maxLatencyLTP + maxLatencyLTD;

                            if (static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->nonIdenticalPulse) { // Non-identical write pulse scheme
                                writeVoltageLTP = static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->VinitLTP + 0.5 * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->VstepLTP * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->maxNumLevelLTP;    // Use average voltage of LTP write voltage
                                writeVoltageLTD = static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->VinitLTD + 0.5 * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->VstepLTD * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->maxNumLevelLTD;    // Use average voltage of LTD write voltage
                            }
                            if (static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->cmosAccess) {  // 1T1R
                                // The energy on selected SLs is included in WriteCell()
                                sumArrayWriteEnergy += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireGateCapRow * tech->vdd * tech->vdd * 2;   // Selected WL (*2 means both LTP and LTD phases)
                                sumArrayWriteEnergy += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapRow * writeVoltageLTP * writeVoltageLTP;   // Selected BL (LTP phases)


                                //TODO: Verification
                                sumArrayWriteEnergy += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapCol * writeVoltageLTP * writeVoltageLTP * (col_size-cells);   // Unselected SLs (LTP phase)


                                // No LTD part because all unselected rows and columns are V=0
                            } else {
                                sumArrayWriteEnergy += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapRow * writeVoltageLTP * writeVoltageLTP;    // Selected WL (LTP phase)
                                sumArrayWriteEnergy += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapRow * writeVoltageLTP/2 * writeVoltageLTP/2 * (row_size - 1);  // Unselected WLs (LTP phase)
                                sumArrayWriteEnergy += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapCol * writeVoltageLTP/2 * writeVoltageLTP/2 * (col_size - cells);   // Unselected BLs (LTP phase)
                                sumArrayWriteEnergy += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapRow * writeVoltageLTD/2 * writeVoltageLTD/2 * (row_size - 1);    // Unselected WLs (LTD phase)
                                sumArrayWriteEnergy += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapCol * writeVoltageLTD/2 * writeVoltageLTD/2 * (col_size - cells); // Unselected BLs (LTD phase)
                            }


                            if (!static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->cmosAccess) { // Cross-point
                                for (int jj = 0; jj < col_size; jj++) { // Half-selected cells in the same row
                                    if (jj >= (col - array_col * operands_per_row) * (cells) && jj < (col - array_col * operands_per_row) * (cells+1)) { continue; } // Skip the selected cells
                                    sumArrayWriteEnergy += (writeVoltageLTP/2 * writeVoltageLTP/2 * static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][jj])->conductanceAtHalfVwLTP * maxLatencyLTP + writeVoltageLTD/2 * writeVoltageLTD/2 * static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][jj])->conductanceAtHalfVwLTD * maxLatencyLTD);
                                }
                                for (int kk = 0; kk < row_size; kk++) {    // Half-selected cells in other rows
                                    // Note that here is a bit inaccurate if using OpenMP, because the weight on other rows (threads) are also being updated
                                    if (kk == row - array_row * row_size) { continue; } // Skip the selected row
                                    for (int b = 0; b < cells; ++b) {
                                        sumArrayWriteEnergy += (writeVoltageLTP/2 * writeVoltageLTP/2 * static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->conductanceAtHalfVwLTP * maxLatencyLTP + writeVoltageLTD/2 * writeVoltageLTD/2 * static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->conductanceAtHalfVwLTD * maxLatencyLTD);
                                    }
                                }
                            }




                        } else {

                            for (int b = 0; b < (cells); b++) {
                                mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->WriteCell(
                                        row - array_row * row_size,
                                        (col - array_col * operands_per_row) * (cells) + b,
                                        (weight_mantissas[col * K + row]
                                                >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1), //Loop. Each loop get the values of bits that would be written into array. Finally, we get the highest cells*bits length values among 24 mantissas.
                                        (weight_mantissas[col * K + row]
                                                >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1),
                                        maxWeight, minWeight,
                                        true);
                                cSum += currentDictionary[(weight_mantissas[col * K + row]
                                                >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1)];
                                cSum += currentDictionary[0];

                                opSumPositive += currentDictionary[(weight_mantissas[col * K + row]
                                                >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1)] * (1<<(bits * ((cells-1) - b)));
                                opSumNegative +=  currentDictionary[0] * (1<<(bits * ((cells-1) - b)));

                                weightChangeBatch = weightChangeBatch || static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse;

                                if(static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTP > maxLatencyLTP)
                                    maxLatencyLTP = static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTP;
                                if(static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTD > maxLatencyLTD)
                                    maxLatencyLTD = static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTD;

                            }
#pragma omp critical
                            {
                                currentSum[col * K + row] += cSum;

                                operandsSumNegative[col * K + row] += opSumNegative;
                                operandsSumPositive[col * K + row] += opSumPositive;
                            }


                            numWriteOperationPerRow += weightChangeBatch;

                            for (int b = 0; b < cells; ++b) {
                                static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTP = maxLatencyLTP;
                                static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTD = maxLatencyLTD;
                                if(weightChangeBatch){
                                    if (static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->nonIdenticalPulse) {	// Non-identical write pulse scheme
                                        if (static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse > 0) {	// LTP
                                            static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTP = sqrt(static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageSquareSum / static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse);	// RMS value of LTP write voltage
                                            static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTD = static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VinitLTD + 0.5 * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VstepLTD * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->maxNumLevelLTD;	// Use average voltage of LTD write voltage
                                        } else if (static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse < 0) {	// LTD
                                            static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTP = static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VinitLTP + 0.5 * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VstepLTP * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->maxNumLevelLTP;    // Use average voltage of LTP write voltage
                                            static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTD = sqrt(static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageSquareSum / (-1*static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse));    // RMS value of LTD write voltage
                                        } else {	// Half-selected during LTP and LTD phases
                                            static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTP = static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VinitLTP + 0.5 * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VstepLTP * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->maxNumLevelLTP;    // Use average voltage of LTP write voltage
                                            static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTD = static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VinitLTD + 0.5 * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VstepLTD * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->maxNumLevelLTD;    // Use average voltage of LTD write voltage
                                        }
                                    }
                                }
                                static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->WriteEnergyCalculation(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapCol);
                                sumArrayWriteEnergy += static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeEnergy;
                            }

                            sumWriteLatencyAnalogNVM += maxLatencyLTP + maxLatencyLTD;

                            if(array_col == 0)
                            wLatency[row] += maxLatencyLTP + maxLatencyLTD;

                            if (static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->nonIdenticalPulse) { // Non-identical write pulse scheme
                                writeVoltageLTP = static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->VinitLTP + 0.5 * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->VstepLTP * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->maxNumLevelLTP;    // Use average voltage of LTP write voltage
                                writeVoltageLTD = static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->VinitLTD + 0.5 * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->VstepLTD * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->maxNumLevelLTD;    // Use average voltage of LTD write voltage
                            }
                            if (static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->cmosAccess) {  // 1T1R
                                // The energy on selected SLs is included in WriteCell()
                                sumArrayWriteEnergy += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireGateCapRow * tech->vdd * tech->vdd * 2;   // Selected WL (*2 means both LTP and LTD phases)
                                sumArrayWriteEnergy += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapRow * writeVoltageLTP * writeVoltageLTP;   // Selected BL (LTP phases)


                                //TODO: Verification
                                sumArrayWriteEnergy += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapCol * writeVoltageLTP * writeVoltageLTP * (col_size-cells);   // Unselected SLs (LTP phase)


                                // No LTD part because all unselected rows and columns are V=0
                            } else {
                                sumArrayWriteEnergy += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapRow * writeVoltageLTP * writeVoltageLTP;    // Selected WL (LTP phase)
                                sumArrayWriteEnergy += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapRow * writeVoltageLTP/2 * writeVoltageLTP/2 * (row_size - 1);  // Unselected WLs (LTP phase)
                                sumArrayWriteEnergy += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapCol * writeVoltageLTP/2 * writeVoltageLTP/2 * (col_size - cells);   // Unselected BLs (LTP phase)
                                sumArrayWriteEnergy += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapRow * writeVoltageLTD/2 * writeVoltageLTD/2 * (row_size - 1);    // Unselected WLs (LTD phase)
                                sumArrayWriteEnergy += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapCol * writeVoltageLTD/2 * writeVoltageLTD/2 * (col_size - cells); // Unselected BLs (LTD phase)
                            }


                            if (!static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->cmosAccess) { // Cross-point
                                for (int jj = 0; jj < col_size; jj++) { // Half-selected cells in the same row
                                    if (jj >= (col - array_col * operands_per_row) * (cells) && jj < (col - array_col * operands_per_row) * (cells+1)) { continue; } // Skip the selected cells
                                    sumArrayWriteEnergy += (writeVoltageLTP/2 * writeVoltageLTP/2 * static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][jj])->conductanceAtHalfVwLTP * maxLatencyLTP + writeVoltageLTD/2 * writeVoltageLTD/2 * static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][jj])->conductanceAtHalfVwLTD * maxLatencyLTD);
                                }
                                for (int kk = 0; kk < row_size; kk++) {    // Half-selected cells in other rows
                                    // Note that here is a bit inaccurate if using OpenMP, because the weight on other rows (threads) are also being updated
                                    if (kk == row - array_row * row_size) { continue; } // Skip the selected row
                                    for (int b = 0; b < cells; ++b) {
                                        sumArrayWriteEnergy += (writeVoltageLTP/2 * writeVoltageLTP/2 * static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->conductanceAtHalfVwLTP * maxLatencyLTP + writeVoltageLTD/2 * writeVoltageLTD/2 * static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->conductanceAtHalfVwLTD * maxLatencyLTD);
                                    }
                                }
                            }
                        }





					}

                    #pragma opm critical    // Use critical here since NeuroSim class functions may update its member variables
                    {
                        if (AnalogNVM *temp = dynamic_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])) {  // Analog eNVM
                            int sumNumWritePulse = 0;
                            for (int j = 0; j < col_size; j++) {
                                sumNumWritePulse += abs(static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][j])->numPulse);    // Note that LTD has negative pulse number
                                sumNumWritePulse += abs(static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][j])->numPulse);    // Note that LTD has negative pulse number
                            }
                            subArrays[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->numWritePulse = sumNumWritePulse / col_size;
                            double writeVoltageSquareSumRow = 0;
                            if (static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->nonIdenticalPulse) { // Non-identical write pulse scheme
                                for (int j = 0; j < col_size; j++) {
                                    writeVoltageSquareSumRow += static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][j])->writeVoltageSquareSum;
                                    writeVoltageSquareSumRow += static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][j])->writeVoltageSquareSum;
                                }
                                if (sumNumWritePulse > 0) {	// Prevent division by 0
                                    subArrays[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell.writeVoltage = sqrt(writeVoltageSquareSumRow / sumNumWritePulse);	// RMS value of write voltage in a row
                                } else {
                                    subArrays[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell.writeVoltage = 0;
                                }
                            }
                        }
                        sumNeuroSimWriteEnergy += NeuroSimSubArrayWriteEnergy(subArrays[(int)(array_row * ceil((double)M/operands_per_row)) + array_col], numWriteOperationPerRow, 0);
                    }

                    numWriteOperation += numWriteOperationPerRow;
                    sumNeuroSimWriteEnergy += NeuroSimSubArrayWriteEnergy(subArrays[(int)(array_row * ceil((double)M/operands_per_row)) + array_col], numWriteOperationPerRow, 0);


                }



			}
		}
        writeEnergy += sumArrayWriteEnergy + sumNeuroSimWriteEnergy;
	}
	else{
	    //IMW. The dimensions of weight are K rows and M(N) columns.

        wLatency = new double[K];
	    signal_matrix = new bool[K * M];
		int operands_per_row = (int)(col_size/(cells));
		unsigned short n = (int)ceil((double)K / row_size) * (int)ceil((double)M / operands_per_row);
		A1 = reinterpret_cast<const unsigned int*>(weights);
        mp = (Array**)malloc(sizeof(Array*)*n);
        mn = (Array**)malloc(sizeof(Array*)*n);
        subArrays = (SubArray**)malloc(sizeof(SubArray*)*n);
		adder = (Adder*)malloc(sizeof(Adder)*n);
		mux = (Mux*)malloc(sizeof(Mux)*n);
		muxDecoder = (RowDecoder*)malloc(sizeof(RowDecoder)*n);
		dff = (DFF*)malloc(sizeof(DFF)*n);
		subtractor = (Subtractor*)malloc(sizeof(Subtractor)*n);
		for (int i = 0; i < K * M; i++) {
			if(weights[i] == 0. || weights[i] == -0.)
			{
				weight_mantissas[i] = 0;
				continue;
			}
			if(get_exponent(A1[i]) == 255)
			{
				printf("weight = %f   Exponent = %d   Error: nan! \n", weights[i], get_exponent(A1[i]));
				exit(1);
			}
			weight_mantissas[i] = get_mantissas(A1[i]);
			signal_matrix[i] = 0;
		}

        for (int k = 0; k < K; ++k) {
            wLatency[k] = 0;
        }
#pragma omp parallel for
        for (int i = 0; i < n; i++) {
            mp[i] = new Array(row_size, col_size, 100);
            mp[i]->Initialization<memoryType>();
            mn[i] = new Array(row_size, col_size, 100);
            mn[i]->Initialization<memoryType>();
            double heightNeuron, widthNeuron;
            new(adder + i) Adder(*inputParameter, *tech, *cell);
            new(mux + i) Mux(*inputParameter, *tech, *cell);
            new(dff + i) DFF(*inputParameter, *tech, *cell);
            new(subtractor + i) Subtractor(*inputParameter, *tech, *cell);
            new(muxDecoder + i) RowDecoder(*inputParameter, *tech, *cell);
            //		subtractor[i] = Subtractor(inputParameter, tech, cell); //RealDevice use Subtractor?
            SubArray *subArray = NULL;
            NeuroSimSubArrayInitialize(subArray, mp[i], *inputParameter, *tech,
                                       *cell);
            NeuroSimSubArrayInitialize(subArray, mn[i], *inputParameter, *tech,
                                       *cell);
            NeuroSimSubArrayArea(subArray);
            NeuroSimSubArrayLeakagePower(subArray);
            NeuroSimNeuronInitialize(subArray, *inputParameter, *tech, *cell, adder[i],
                                     mux[i], muxDecoder[i], dff[i], subtractor[i]);
            NeuroSimNeuronArea(subArray, adder[i], mux[i], muxDecoder[i], dff[i], subtractor[i],
                               &heightNeuron, &widthNeuron);
            subArrays[i] = subArray;
        }


        double sumArrayWriteEnergy = 0;   // Use a temporary variable here since OpenMP does not support reduction on class member
        double sumNeuroSimWriteEnergy = 0;   // Use a temporary variable here since OpenMP does not support reduction on class member
        double sumWriteLatencyAnalogNVM = 0;   // Use a temporary variable here since OpenMP does not support reduction on class member

        int array_size = ceil((double)K / row_size);
#pragma omp parallel for
        for (int array_row = 0; array_row < array_size; array_row++) {
			for (int array_col = 0; array_col < ceil((double)M / operands_per_row); array_col++) { // array_col: Number of arrays that contains all the columns.
                for (int col = array_col * operands_per_row;
						col < M && col < operands_per_row * (array_col + 1); col++) { // col: The index of columns
					for (int row = array_row * row_size;
							row < K && row < row_size * (array_row + 1); row++) { // row: The index of columns
						//reset weight1
						//We modify the bound of weight from -1,1 to 0,1, since weight can only present positive value.
						for (int b = 0; b < (cells); b++) {
                            mp[(int)(array_row * ceil((double)M/operands_per_row))+array_col]->WriteCell(row - array_row * row_size, (col - array_col * operands_per_row) * (cells) + b,-(maxWeight - minWeight),
                                    0 /* delta_W=-(param->maxWeight-param->minWeight) will completely erase */,
                                    maxWeight, minWeight, true);
                            mn[(int)(array_row * ceil((double)M/operands_per_row))+array_col]->WriteCell(row - array_row * row_size, (col - array_col * operands_per_row) * (cells) + b,-(maxWeight - minWeight),
                                    0 /* delta_W=-(param->maxWeight-param->minWeight) will completely erase */,
                                    maxWeight, minWeight, true);

						}
					}
				}
				#pragma omp critical
                {
				for (int row = 0; row < K; row++) {
					for (int col = 0; col < M; col++) {
						unsigned char exponentA = get_exponent(A1[row * M + col]); //Only map weight1, the inputs are aligned in computation step, this may cause loss
						if (row - array_row * row_size == 0
								|| exponentA > maxExponentA[array_row][col])
							maxExponentA[array_row][col] = exponentA;
					}
				}
				for (int row = array_row * row_size; row < K && row < row_size * (array_row + 1); row++) {
					for (int col = array_col * operands_per_row; col < M && col < operands_per_row * (array_col + 1); col++) {
						unsigned char exponentA = get_exponent(A1[row * M + col]);
						char diff = maxExponentA[array_row][col] - exponentA;
						//adjust the weight1
						weight_mantissas[row * M + col] = weight_mantissas[row * M + col]
								>> diff;
                    }
                }
                }

            }
        }



#pragma omp parallel
		for (int array_row = 0; array_row < array_size; array_row++) {
			for (int array_col = 0; array_col < ceil((double)M / operands_per_row); array_col++) {

				for (int row = array_row * row_size;
                     row < K && row < row_size * (array_row + 1); row++) {


                    int numWriteOperationPerRow = 0; // # of write operations will be performed per row.
                    double numWriteOperation = 0;


                    double writeVoltageLTP = static_cast<eNVM*>(mp[0]->cell[0][0])->writeVoltageLTP;
                    double writeVoltageLTD = static_cast<eNVM*>(mp[0]->cell[0][0])->writeVoltageLTD;

					for (int col = array_col * operands_per_row; col < M && col < operands_per_row * (array_col + 1); col++) {

                        bool weightChangeBatch = false;
                        double maxLatencyLTP = 0;
                        double maxLatencyLTD = 0;
                        double cSum = 0,opSumNegative = 0,opSumPositive = 0;
                        if (weights[row * M + col] < 0) //if weight1 is negative
                        {
                            signal_matrix[row * M + col] = 1; //Useless
                            for (int b = 0; b < (cells); b++) {
                                mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->WriteCell(
                                        row - array_row * row_size,
                                        (col - array_col * operands_per_row) * (cells) + b,
                                        (weight_mantissas[row * M + col]
                                                >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1),
                                        (weight_mantissas[row * M + col]
                                                >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1),
                                        maxWeight, minWeight,
                                        true);

                                cSum += currentDictionary[(weight_mantissas[col * K + row]
                                                >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1)];
                                cSum += currentDictionary[0];

                                opSumNegative += currentDictionary[(weight_mantissas[col * K + row]
                                                >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1)] * (1<<(bits * ((cells-1) - b)));

                                opSumPositive +=  currentDictionary[0] * (1<<(bits * ((cells-1) - b)));





                                weightChangeBatch = weightChangeBatch || static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse;

                                if(static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTP > maxLatencyLTP)
                                    maxLatencyLTP = static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTP;
                                if(static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTD > maxLatencyLTD)
                                    maxLatencyLTD = static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTD;

                            }
#pragma omp critical
                            {
                                currentSum[col * K + row] += cSum;

                                operandsSumNegative[col * K + row] += opSumNegative;
                                operandsSumPositive[col * K + row] += opSumPositive;
                            }

                            numWriteOperationPerRow += weightChangeBatch;

                            for (int b = 0; b < cells; ++b) {
                                static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTP = maxLatencyLTP;
                                static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTD = maxLatencyLTD;
                                if(weightChangeBatch){
                                    if (static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->nonIdenticalPulse) {	// Non-identical write pulse scheme
                                        if (static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse > 0) {	// LTP
                                            static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTP = sqrt(static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageSquareSum / static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse);	// RMS value of LTP write voltage
                                            static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTD = static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VinitLTD + 0.5 * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VstepLTD * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->maxNumLevelLTD;	// Use average voltage of LTD write voltage
                                        } else if (static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse < 0) {	// LTD
                                            static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTP = static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VinitLTP + 0.5 * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VstepLTP * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->maxNumLevelLTP;    // Use average voltage of LTP write voltage
                                            static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTD = sqrt(static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageSquareSum / (-1*static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse));    // RMS value of LTD write voltage
                                        } else {	// Half-selected during LTP and LTD phases
                                            static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTP = static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VinitLTP + 0.5 * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VstepLTP * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->maxNumLevelLTP;    // Use average voltage of LTP write voltage
                                            static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTD = static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VinitLTD + 0.5 * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VstepLTD * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->maxNumLevelLTD;    // Use average voltage of LTD write voltage
                                        }
                                    }
                                }
                                static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->WriteEnergyCalculation(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapCol);
                                sumArrayWriteEnergy += static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeEnergy;
                            }

                            sumWriteLatencyAnalogNVM += maxLatencyLTP + maxLatencyLTD;

                            if(array_col == 0)
                            wLatency[row] += maxLatencyLTP + maxLatencyLTD;

                            if (static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->nonIdenticalPulse) { // Non-identical write pulse scheme
                                writeVoltageLTP = static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->VinitLTP + 0.5 * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->VstepLTP * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->maxNumLevelLTP;    // Use average voltage of LTP write voltage
                                writeVoltageLTD = static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->VinitLTD + 0.5 * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->VstepLTD * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->maxNumLevelLTD;    // Use average voltage of LTD write voltage
                            }
                            if (static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->cmosAccess) {  // 1T1R
                                // The energy on selected SLs is included in WriteCell()
                                sumArrayWriteEnergy += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireGateCapRow * tech->vdd * tech->vdd * 2;   // Selected WL (*2 means both LTP and LTD phases)
                                sumArrayWriteEnergy += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapRow * writeVoltageLTP * writeVoltageLTP;   // Selected BL (LTP phases)


                                //TODO: Verification
                                sumArrayWriteEnergy += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapCol * writeVoltageLTP * writeVoltageLTP * (col_size-cells);   // Unselected SLs (LTP phase)


                                // No LTD part because all unselected rows and columns are V=0
                            } else {
                                sumArrayWriteEnergy += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapRow * writeVoltageLTP * writeVoltageLTP;    // Selected WL (LTP phase)
                                sumArrayWriteEnergy += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapRow * writeVoltageLTP/2 * writeVoltageLTP/2 * (row_size - 1);  // Unselected WLs (LTP phase)
                                sumArrayWriteEnergy += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapCol * writeVoltageLTP/2 * writeVoltageLTP/2 * (col_size - cells);   // Unselected BLs (LTP phase)
                                sumArrayWriteEnergy += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapRow * writeVoltageLTD/2 * writeVoltageLTD/2 * (row_size - 1);    // Unselected WLs (LTD phase)
                                sumArrayWriteEnergy += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapCol * writeVoltageLTD/2 * writeVoltageLTD/2 * (col_size - cells); // Unselected BLs (LTD phase)
                            }


                            if (!static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->cmosAccess) { // Cross-point
                                for (int jj = 0; jj < col_size; jj++) { // Half-selected cells in the same row
                                    if (jj >= (col - array_col * operands_per_row) * (cells) && jj < (col - array_col * operands_per_row) * (cells+1)) { continue; } // Skip the selected cells
                                    sumArrayWriteEnergy += (writeVoltageLTP/2 * writeVoltageLTP/2 * static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][jj])->conductanceAtHalfVwLTP * maxLatencyLTP + writeVoltageLTD/2 * writeVoltageLTD/2 * static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][jj])->conductanceAtHalfVwLTD * maxLatencyLTD);
                                }
                                for (int kk = 0; kk < row_size; kk++) {    // Half-selected cells in other rows
                                    // Note that here is a bit inaccurate if using OpenMP, because the weight on other rows (threads) are also being updated
                                    if (kk == row - array_row * row_size) { continue; } // Skip the selected row
                                    for (int b = 0; b < cells; ++b) {
                                        sumArrayWriteEnergy += (writeVoltageLTP/2 * writeVoltageLTP/2 * static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->conductanceAtHalfVwLTP * maxLatencyLTP + writeVoltageLTD/2 * writeVoltageLTD/2 * static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->conductanceAtHalfVwLTD * maxLatencyLTD);
                                    }
                                }
                            }

                        } else {
                            for (int b = 0; b < (cells); b++) {
                                mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->WriteCell(
                                        row - array_row * row_size,
                                        (col - array_col * operands_per_row) * (cells) + b,
                                        (weight_mantissas[row * M + col]
                                                >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1),
                                        (weight_mantissas[row * M + col]
                                                >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1),
                                        maxWeight, minWeight,
                                        true);

                                cSum += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->ReadCell(row - array_row * row_size,(col - array_col * operands_per_row) * (cells) + b);
                                cSum += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->ReadCell(row - array_row * row_size,(col - array_col * operands_per_row) * (cells) + b);

                                opSumNegative += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->ReadCell(row - array_row * row_size,(col - array_col * operands_per_row) * (cells) + b) * (1<<(bits * ((cells-1) - b)));
                                opSumPositive += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->ReadCell(row - array_row * row_size,(col - array_col * operands_per_row) * (cells) + b) * (1<<(bits * ((cells-1) - b)));


                                cSum += currentDictionary[(weight_mantissas[col * K + row]
                                                >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1)];
                                cSum += currentDictionary[0];

                                opSumPositive += currentDictionary[(weight_mantissas[col * K + row]
                                                >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1)] * (1<<(bits * ((cells-1) - b)));

                                opSumNegative +=  currentDictionary[0] * (1<<(bits * ((cells-1) - b)));



                                weightChangeBatch = weightChangeBatch || static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse;

                                if(static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTP > maxLatencyLTP)
                                    maxLatencyLTP = static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTP;
                                if(static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTD > maxLatencyLTD)
                                    maxLatencyLTD = static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTD;

                            }
#pragma omp critical
                            {
                                currentSum[col * K + row] += cSum;

                                operandsSumNegative[col * K + row] += opSumNegative;
                                operandsSumPositive[col * K + row] += opSumPositive;
                            }


                            numWriteOperationPerRow += weightChangeBatch;

                            for (int b = 0; b < cells; ++b) {
                                static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTP = maxLatencyLTP;
                                static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTD = maxLatencyLTD;
                                if(weightChangeBatch){
                                    if (static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->nonIdenticalPulse) {	// Non-identical write pulse scheme
                                        if (static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse > 0) {	// LTP
                                            static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTP = sqrt(static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageSquareSum / static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse);	// RMS value of LTP write voltage
                                            static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTD = static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VinitLTD + 0.5 * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VstepLTD * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->maxNumLevelLTD;	// Use average voltage of LTD write voltage
                                        } else if (static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse < 0) {	// LTD
                                            static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTP = static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VinitLTP + 0.5 * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VstepLTP * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->maxNumLevelLTP;    // Use average voltage of LTP write voltage
                                            static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTD = sqrt(static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageSquareSum / (-1*static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse));    // RMS value of LTD write voltage
                                        } else {	// Half-selected during LTP and LTD phases
                                            static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTP = static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VinitLTP + 0.5 * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VstepLTP * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->maxNumLevelLTP;    // Use average voltage of LTP write voltage
                                            static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTD = static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VinitLTD + 0.5 * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VstepLTD * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->maxNumLevelLTD;    // Use average voltage of LTD write voltage
                                        }
                                    }
                                }
                                static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->WriteEnergyCalculation(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapCol);
                                sumArrayWriteEnergy += static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeEnergy;
                            }

                            sumWriteLatencyAnalogNVM += maxLatencyLTP + maxLatencyLTD;

                            if(array_col == 0)
                            wLatency[row] += maxLatencyLTP + maxLatencyLTD;

                            if (static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->nonIdenticalPulse) { // Non-identical write pulse scheme
                                writeVoltageLTP = static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->VinitLTP + 0.5 * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->VstepLTP * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->maxNumLevelLTP;    // Use average voltage of LTP write voltage
                                writeVoltageLTD = static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->VinitLTD + 0.5 * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->VstepLTD * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->maxNumLevelLTD;    // Use average voltage of LTD write voltage
                            }
                            if (static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->cmosAccess) {  // 1T1R
                                // The energy on selected SLs is included in WriteCell()
                                sumArrayWriteEnergy += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireGateCapRow * tech->vdd * tech->vdd * 2;   // Selected WL (*2 means both LTP and LTD phases)
                                sumArrayWriteEnergy += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapRow * writeVoltageLTP * writeVoltageLTP;   // Selected BL (LTP phases)


                                //TODO: Verification
                                sumArrayWriteEnergy += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapCol * writeVoltageLTP * writeVoltageLTP * (col_size-cells);   // Unselected SLs (LTP phase)


                                // No LTD part because all unselected rows and columns are V=0
                            } else {
                                sumArrayWriteEnergy += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapRow * writeVoltageLTP * writeVoltageLTP;    // Selected WL (LTP phase)
                                sumArrayWriteEnergy += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapRow * writeVoltageLTP/2 * writeVoltageLTP/2 * (row_size - 1);  // Unselected WLs (LTP phase)
                                sumArrayWriteEnergy += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapCol * writeVoltageLTP/2 * writeVoltageLTP/2 * (col_size - cells);   // Unselected BLs (LTP phase)
                                sumArrayWriteEnergy += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapRow * writeVoltageLTD/2 * writeVoltageLTD/2 * (row_size - 1);    // Unselected WLs (LTD phase)
                                sumArrayWriteEnergy += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapCol * writeVoltageLTD/2 * writeVoltageLTD/2 * (col_size - cells); // Unselected BLs (LTD phase)
                            }


                            if (!static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->cmosAccess) { // Cross-point
                                for (int jj = 0; jj < col_size; jj++) { // Half-selected cells in the same row
                                    if (jj >= (col - array_col * operands_per_row) * (cells) && jj < (col - array_col * operands_per_row) * (cells+1)) { continue; } // Skip the selected cells
                                    sumArrayWriteEnergy += (writeVoltageLTP/2 * writeVoltageLTP/2 * static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][jj])->conductanceAtHalfVwLTP * maxLatencyLTP + writeVoltageLTD/2 * writeVoltageLTD/2 * static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][jj])->conductanceAtHalfVwLTD * maxLatencyLTD);
                                }
                                for (int kk = 0; kk < row_size; kk++) {    // Half-selected cells in other rows
                                    // Note that here is a bit inaccurate if using OpenMP, because the weight on other rows (threads) are also being updated
                                    if (kk == row - array_row * row_size) { continue; } // Skip the selected row
                                    for (int b = 0; b < cells; ++b) {
                                        sumArrayWriteEnergy += (writeVoltageLTP/2 * writeVoltageLTP/2 * static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->conductanceAtHalfVwLTP * maxLatencyLTP + writeVoltageLTD/2 * writeVoltageLTD/2 * static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->conductanceAtHalfVwLTD * maxLatencyLTD);
                                    }
                                }
                            }
                        }


                        #pragma opm critical    // Use critical here since NeuroSim class functions may update its member variables
                        {
                            if (AnalogNVM *temp = dynamic_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])) {  // Analog eNVM
                                int sumNumWritePulse = 0;
                                for (int j = 0; j < col_size; j++) {
                                    sumNumWritePulse += abs(static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][j])->numPulse);    // Note that LTD has negative pulse number
                                    sumNumWritePulse += abs(static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][j])->numPulse);    // Note that LTD has negative pulse number
                                }
                                subArrays[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->numWritePulse = sumNumWritePulse / col_size;
                                double writeVoltageSquareSumRow = 0;
                                if (static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->nonIdenticalPulse) { // Non-identical write pulse scheme
                                    for (int j = 0; j < col_size; j++) {
                                        writeVoltageSquareSumRow += static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][j])->writeVoltageSquareSum;
                                        writeVoltageSquareSumRow += static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][j])->writeVoltageSquareSum;
                                    }
                                    if (sumNumWritePulse > 0) {	// Prevent division by 0
                                        subArrays[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell.writeVoltage = sqrt(writeVoltageSquareSumRow / sumNumWritePulse);	// RMS value of write voltage in a row
                                    } else {
                                        subArrays[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell.writeVoltage = 0;
                                    }
                                }
                            }
                            sumNeuroSimWriteEnergy += NeuroSimSubArrayWriteEnergy(subArrays[(int)(array_row * ceil((double)M/operands_per_row)) + array_col], numWriteOperationPerRow, 0);
                        }

                        numWriteOperation += numWriteOperationPerRow;
                        sumNeuroSimWriteEnergy += NeuroSimSubArrayWriteEnergy(subArrays[(int)(array_row * ceil((double)M/operands_per_row)) + array_col], numWriteOperationPerRow, 0);

                    }
				}
			}
		}
        writeEnergy += sumArrayWriteEnergy + sumNeuroSimWriteEnergy;
	}

    delete[] currentDictionary;
	A1 = NULL;
	delete[] weight_mantissas;
}

template<class memoryType>
void Memristor_CU<memoryType>::re_map(uint32_t MM, uint32_t KK, const float *weights, bool ot){

    const unsigned int *A1;

    if(currentSum != NULL)
        delete[] currentSum;
    Array *tp = new Array(1,1<<bits,100);
    tp->Initialization<memoryType>();
    double *currentDictionary = new double[1<<bits];

    for(int i = 0; i < 1<<bits; i++){
        tp->WriteCell(0,i,-(maxWeight - minWeight) ,0,maxWeight,minWeight,true);
        tp->WriteCell(0,i,i,i,maxWeight,minWeight,true);
        currentDictionary[i] = tp->ReadCell(0,i);
    }

    delete tp;
    currentSum = new double[M*K];
    unsigned int *weight_mantissas = new unsigned int[M*K];
    for (int j = 0; j < M * K; ++j) {
        currentSum[j] = 0;
        operandsSumPositive[j] = 0;
        operandsSumNegative[j] = 0;
    }
    delete[] maxExponentA;
    maxExponentA = new unsigned char*[(int)(K/row_size)+1];
    for(int i = 0 ; i < (int)(K/row_size)+1; i++)
    {
        maxExponentA[i] = new unsigned char[M];
    }
    int operands_per_row = (int)(col_size/(cells));
    int n1 = (int)ceil((double)K / row_size) * (int)ceil((double)M / operands_per_row);


    wLatency = new double[K];

    this->M = MM;
    this->K = KK;
    unsigned short n = (int)ceil((double)K / row_size) * (int)ceil((double)M / operands_per_row);
    A1 = reinterpret_cast<const unsigned int*>(weights);

    for (int i = 0; i < K * M; i++) {
        if (weights[i] == 0. || weights[i] == -0.) {
            weight_mantissas[i] = 0;
            continue;
        }
        if (get_exponent(A1[i]) == 255) {
            printf("weight = %f   Exponent = %d   Error: nan! \n", weights[i], get_exponent(A1[i]));
            exit(1);
        }
        weight_mantissas[i] = get_mantissas(A1[i]);
    }
    for (int k = 0; k < K; ++k) {
        wLatency[k] = 0;
    }




    if(n > n1) {
        mp = (Array **) realloc(mp, sizeof(Array *) * n);
        mn = (Array **) realloc(mn, sizeof(Array *) * n);

        subArrays = (SubArray **) realloc(subArrays, sizeof(SubArray *) * n);
        adder = (Adder *) realloc(adder, sizeof(Adder) * n);
        mux = (Mux *) realloc(mux, sizeof(Mux) * n);
        muxDecoder = (RowDecoder *) realloc(muxDecoder, sizeof(RowDecoder) * n);
        dff = (DFF *) realloc(dff, sizeof(DFF) * n);
        subtractor = (Subtractor *) realloc(subtractor, sizeof(Subtractor) * n);

#pragma omp parallel for
        for (int i = n1; i < n; i++) {
            mp[i] = new Array(row_size, col_size, 100);
            mp[i]->Initialization<memoryType>();
            mn[i] = new Array(row_size, col_size, 100);
            mn[i]->Initialization<memoryType>();
            double heightNeuron, widthNeuron;
            new(adder + i) Adder(*inputParameter, *tech, *cell);
            new(mux + i) Mux(*inputParameter, *tech, *cell);
            new(dff + i) DFF(*inputParameter, *tech, *cell);
            new(subtractor + i) Subtractor(*inputParameter, *tech, *cell);
            new(muxDecoder + i) RowDecoder(*inputParameter, *tech, *cell);
//			subtractor[i] = Subtractor(inputParameter, tech, cell); //RealDevice use Subtractor?
            SubArray *subArray = NULL;
            NeuroSimSubArrayInitialize(subArray, mp[i], *inputParameter, *tech,
                                       *cell);
            NeuroSimSubArrayInitialize(subArray, mn[i], *inputParameter, *tech,
                                       *cell);

            NeuroSimSubArrayArea(subArray);
            NeuroSimSubArrayLeakagePower(subArray);
            NeuroSimNeuronInitialize(subArray, *inputParameter, *tech, *cell, adder[i],
                                     mux[i], muxDecoder[i], dff[i], subtractor[i]);
            NeuroSimNeuronArea(subArray, adder[i], mux[i], muxDecoder[i], dff[i], subtractor[i],
                               &heightNeuron, &widthNeuron);
            subArrays[i] = subArray;
            
            for(int j = 0; j < row_size; j++){
                for (int k = 0; k < col_size; ++k) {
                    mp[i]->WriteCell(j,k,-(maxWeight-minWeight),0,maxWeight,minWeight,true);
                    mn[i]->WriteCell(j,k,-(maxWeight-minWeight),0,maxWeight,minWeight,true);
                }
            }
        }
    }

    if(ot == WMI){
        double sumArrayWriteEnergy = 0;   // Use a temporary variable here since OpenMP does not support reduction on class member
        double sumNeuroSimWriteEnergy = 0;   // Use a temporary variable here since OpenMP does not support reduction on class member
        double sumWriteLatencyAnalogNVM = 0;   // Use a temporary variable here since OpenMP does not support reduction on class member




        int array_size = ceil((double)K / row_size); //Number of groups. Each group contains all the columns.
#pragma omp parallel for reduction(+: sumArrayWriteEnergy, sumNeuroSimWriteEnergy, sumWriteLatencyAnalogNVM)
        for (int array_row = 0; array_row < array_size; array_row++) {
            for (int array_col = 0; array_col < ceil((double)M / operands_per_row); array_col++) { // array_col: Number of arrays that contains all the columns.

#pragma omp critical
                {
                //Get the max exponent of weights that mapped in the same column.
                for (int row = 0; row < K; row++) {
                    for (int col = 0; col < M; col++) {
                        unsigned char exponentA = get_exponent(A1[col * K + row]); //Only map weight, the inputs are aligned in computation step, this may cause loss
                        if (row - array_row * row_size == 0
                            || exponentA > maxExponentA[array_row][col])
                            maxExponentA[array_row][col] = exponentA;
                    }
                }

                //Align the operands
                for (int row = array_row * row_size; row < K && row < row_size * (array_row + 1); row++) {
                    for (int col = array_col * operands_per_row; col <M && col < operands_per_row * (array_col + 1); col++) {
                        unsigned char exponentA = get_exponent(A1[col * K + row]);
                        char diff = maxExponentA[array_row][col] - exponentA;
                        //adjust the weight1
                        weight_mantissas[col * K + row] = weight_mantissas[col * K + row]
                                >> diff;
                    }
                }
                }



                for (int row = array_row * row_size;
                     row < KK && row < row_size * (array_row + 1); row++) {

                    int numWriteOperationPerRow = 0; // # of write operations will be performed per row.
                    double numWriteOperation = 0;


                    double writeVoltageLTP = static_cast<eNVM*>(mp[0]->cell[0][0])->writeVoltageLTP;
                    double writeVoltageLTD = static_cast<eNVM*>(mp[0]->cell[0][0])->writeVoltageLTD;


                    for (int col = array_col * operands_per_row; col < M && col < operands_per_row * (array_col + 1);
                         col++) {

                        bool weightChangeBatch = false;
                        double maxLatencyLTP = 0;
                        double maxLatencyLTD = 0;

                        double cSum = 0,opSumNegative = 0,opSumPositive = 0;


                        if (weights[col * K + row] < 0) //if weight1 is negative
                        {
                            for (int b = 0; b < (cells); b++) {  //Assume cells*colmux=colsize

                                if(((int)(array_row * ceil((double)M/operands_per_row)) + array_col) < n1){
                                    int old_weight = mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->getWeight(row - array_row * row_size,(col - array_col * operands_per_row) * (cells) + b);
                                    int new_weight = (weight_mantissas[col * K + row]
                                            >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1);
//                                    if(col > 15)
//                                        printf("%d %d \n",old_weight,new_weight);
//                                    if(new_weight == old_weight)
//                                        ;
//                                    else
                                        mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->WriteCell(
                                                row - array_row * row_size,
                                                (col - array_col * operands_per_row) * (cells) + b,
                                                new_weight-old_weight,
                                                new_weight,
                                                maxWeight, minWeight,
                                                true);

                                }else{
                                    //For WMI, we map the weight by transposed matrix
                                    mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->WriteCell(
                                            row - array_row * row_size,
                                            (col - array_col * operands_per_row) * (cells) + b,
                                            (weight_mantissas[col * K + row]
                                                    >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1),
                                            (weight_mantissas[col * K + row]
                                                    >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1),
                                            maxWeight, minWeight,
                                            true);
                                }


                                cSum += currentDictionary[(weight_mantissas[col * K + row]
                                                >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1)];
                                cSum += currentDictionary[0];

                                opSumNegative += currentDictionary[(weight_mantissas[col * K + row]
                                                >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1)] * (1<<(bits * ((cells-1) - b)));

                                opSumPositive +=  currentDictionary[0] * (1<<(bits * ((cells-1) - b)));




                                weightChangeBatch = weightChangeBatch || static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse;

                                if(static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTP > maxLatencyLTP)
                                    maxLatencyLTP = static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTP;
                                if(static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTD > maxLatencyLTD)
                                    maxLatencyLTD = static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTD;

                            }
#pragma omp critical
                            {
                                currentSum[col * K + row] += cSum;

                                operandsSumNegative[col * K + row] += opSumNegative;
                                operandsSumPositive[col * K + row] += opSumPositive;
                            }


                            numWriteOperationPerRow += weightChangeBatch;

                            for (int b = 0; b < cells; ++b) {
                                static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTP = maxLatencyLTP;
                                static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTD = maxLatencyLTD;
                                if(weightChangeBatch){
                                    if (static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->nonIdenticalPulse) {	// Non-identical write pulse scheme
                                        if (static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse > 0) {	// LTP
                                            static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTP = sqrt(static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageSquareSum / static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse);	// RMS value of LTP write voltage
                                            static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTD = static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VinitLTD + 0.5 * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VstepLTD * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->maxNumLevelLTD;	// Use average voltage of LTD write voltage
                                        } else if (static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse < 0) {	// LTD
                                            static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTP = static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VinitLTP + 0.5 * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VstepLTP * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->maxNumLevelLTP;    // Use average voltage of LTP write voltage
                                            static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTD = sqrt(static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageSquareSum / (-1*static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse));    // RMS value of LTD write voltage
                                        } else {	// Half-selected during LTP and LTD phases
                                            static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTP = static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VinitLTP + 0.5 * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VstepLTP * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->maxNumLevelLTP;    // Use average voltage of LTP write voltage
                                            static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTD = static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VinitLTD + 0.5 * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VstepLTD * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->maxNumLevelLTD;    // Use average voltage of LTD write voltage
                                        }
                                    }
                                }
                                static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->WriteEnergyCalculation(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapCol);
                                sumArrayWriteEnergy += static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeEnergy;
                            }

                            sumWriteLatencyAnalogNVM += maxLatencyLTP + maxLatencyLTD;

                            if(array_col == 0)
                            wLatency[row] += maxLatencyLTP + maxLatencyLTD;

                            if (static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->nonIdenticalPulse) { // Non-identical write pulse scheme
                                writeVoltageLTP = static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->VinitLTP + 0.5 * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->VstepLTP * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->maxNumLevelLTP;    // Use average voltage of LTP write voltage
                                writeVoltageLTD = static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->VinitLTD + 0.5 * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->VstepLTD * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->maxNumLevelLTD;    // Use average voltage of LTD write voltage
                            }
                            if (static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->cmosAccess) {  // 1T1R
                                // The energy on selected SLs is included in WriteCell()
                                sumArrayWriteEnergy += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireGateCapRow * tech->vdd * tech->vdd * 2;   // Selected WL (*2 means both LTP and LTD phases)
                                sumArrayWriteEnergy += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapRow * writeVoltageLTP * writeVoltageLTP;   // Selected BL (LTP phases)


                                //TODO: Verification
                                sumArrayWriteEnergy += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapCol * writeVoltageLTP * writeVoltageLTP * (col_size-cells);   // Unselected SLs (LTP phase)


                                // No LTD part because all unselected rows and columns are V=0
                            } else {
                                sumArrayWriteEnergy += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapRow * writeVoltageLTP * writeVoltageLTP;    // Selected WL (LTP phase)
                                sumArrayWriteEnergy += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapRow * writeVoltageLTP/2 * writeVoltageLTP/2 * (row_size - 1);  // Unselected WLs (LTP phase)
                                sumArrayWriteEnergy += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapCol * writeVoltageLTP/2 * writeVoltageLTP/2 * (col_size - cells);   // Unselected BLs (LTP phase)
                                sumArrayWriteEnergy += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapRow * writeVoltageLTD/2 * writeVoltageLTD/2 * (row_size - 1);    // Unselected WLs (LTD phase)
                                sumArrayWriteEnergy += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapCol * writeVoltageLTD/2 * writeVoltageLTD/2 * (col_size - cells); // Unselected BLs (LTD phase)
                            }


                            if (!static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->cmosAccess) { // Cross-point
                                for (int jj = 0; jj < col_size; jj++) { // Half-selected cells in the same row
                                    if (jj >= (col - array_col * operands_per_row) * (cells) && jj < (col - array_col * operands_per_row) * (cells+1)) { continue; } // Skip the selected cells
                                    sumArrayWriteEnergy += (writeVoltageLTP/2 * writeVoltageLTP/2 * static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][jj])->conductanceAtHalfVwLTP * maxLatencyLTP + writeVoltageLTD/2 * writeVoltageLTD/2 * static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][jj])->conductanceAtHalfVwLTD * maxLatencyLTD);
                                }
                                for (int kk = 0; kk < row_size; K++) {    // Half-selected cells in other rows
                                    // Note that here is a bit inaccurate if using OpenMP, because the weight on other rows (threads) are also being updated
                                    if (kk == row - array_row * row_size) { continue; } // Skip the selected row
                                    for (int b = 0; b < cells; ++b) {
                                        sumArrayWriteEnergy += (writeVoltageLTP/2 * writeVoltageLTP/2 * static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->conductanceAtHalfVwLTP * maxLatencyLTP + writeVoltageLTD/2 * writeVoltageLTD/2 * static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->conductanceAtHalfVwLTD * maxLatencyLTD);
                                    }
                                }
                            }




                        } else {
                            for (int b = 0; b < (cells); b++) {

                                if(((int)(array_row * ceil((double)M/operands_per_row)) + array_col) < n1){
                                    int old_weight = mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->getWeight(row - array_row * row_size,(col - array_col * operands_per_row) * (cells) + b);
                                    int new_weight = (weight_mantissas[col * K + row]
                                            >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1);
//                                    if(new_weight == old_weight)
//                                        ;
//                                    else
                                        mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->WriteCell(
                                            row - array_row * row_size,
                                            (col - array_col * operands_per_row) * (cells) + b,
                                            new_weight-old_weight,
                                            new_weight,
                                            maxWeight, minWeight,
                                            true);
                                }else{
                                    mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->WriteCell(
                                            row - array_row * row_size,
                                            (col - array_col * operands_per_row) * (cells) + b,
                                            (weight_mantissas[col * K + row]
                                                    >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1), //Loop. Each loop get the values of bits that would be written into array. Finally, we get the highest cells*bits length values among 24 mantissas.
                                            (weight_mantissas[col * K + row]
                                                    >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1),
                                            maxWeight, minWeight,
                                            true);
                                }

                                cSum += currentDictionary[(weight_mantissas[col * K + row]
                                                >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1)];
                                cSum += currentDictionary[0];

                                opSumPositive += currentDictionary[(weight_mantissas[col * K + row]
                                                >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1)] * (1<<(bits * ((cells-1) - b)));

                                opSumNegative +=  currentDictionary[0] * (1<<(bits * ((cells-1) - b)));




                                weightChangeBatch = weightChangeBatch || static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse;

                                if(static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTP > maxLatencyLTP)
                                    maxLatencyLTP = static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTP;
                                if(static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTD > maxLatencyLTD)
                                    maxLatencyLTD = static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTD;

                            }
#pragma omp critical
                            {
                                currentSum[col * K + row] += cSum;

                                operandsSumNegative[col * K + row] += opSumNegative;
                                operandsSumPositive[col * K + row] += opSumPositive;
                            }


                            numWriteOperationPerRow += weightChangeBatch;

                            for (int b = 0; b < cells; ++b) {
                                static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTP = maxLatencyLTP;
                                static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTD = maxLatencyLTD;
                                if(weightChangeBatch){
                                    if (static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->nonIdenticalPulse) {	// Non-identical write pulse scheme
                                        if (static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse > 0) {	// LTP
                                            static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTP = sqrt(static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageSquareSum / static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse);	// RMS value of LTP write voltage
                                            static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTD = static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VinitLTD + 0.5 * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VstepLTD * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->maxNumLevelLTD;	// Use average voltage of LTD write voltage
                                        } else if (static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse < 0) {	// LTD
                                            static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTP = static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VinitLTP + 0.5 * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VstepLTP * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->maxNumLevelLTP;    // Use average voltage of LTP write voltage
                                            static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTD = sqrt(static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageSquareSum / (-1*static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse));    // RMS value of LTD write voltage
                                        } else {	// Half-selected during LTP and LTD phases
                                            static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTP = static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VinitLTP + 0.5 * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VstepLTP * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->maxNumLevelLTP;    // Use average voltage of LTP write voltage
                                            static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTD = static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VinitLTD + 0.5 * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VstepLTD * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->maxNumLevelLTD;    // Use average voltage of LTD write voltage
                                        }
                                    }
                                }
                                static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->WriteEnergyCalculation(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapCol);
                                sumArrayWriteEnergy += static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeEnergy;
                            }

                            sumWriteLatencyAnalogNVM += maxLatencyLTP + maxLatencyLTD;

                            if(array_col == 0)
                            wLatency[row] += maxLatencyLTP + maxLatencyLTD;

                            if (static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->nonIdenticalPulse) { // Non-identical write pulse scheme
                                writeVoltageLTP = static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->VinitLTP + 0.5 * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->VstepLTP * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->maxNumLevelLTP;    // Use average voltage of LTP write voltage
                                writeVoltageLTD = static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->VinitLTD + 0.5 * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->VstepLTD * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->maxNumLevelLTD;    // Use average voltage of LTD write voltage
                            }
                            if (static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->cmosAccess) {  // 1T1R
                                // The energy on selected SLs is included in WriteCell()
                                sumArrayWriteEnergy += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireGateCapRow * tech->vdd * tech->vdd * 2;   // Selected WL (*2 means both LTP and LTD phases)
                                sumArrayWriteEnergy += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapRow * writeVoltageLTP * writeVoltageLTP;   // Selected BL (LTP phases)


                                //TODO: Verification
                                sumArrayWriteEnergy += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapCol * writeVoltageLTP * writeVoltageLTP * (col_size-cells);   // Unselected SLs (LTP phase)


                                // No LTD part because all unselected rows and columns are V=0
                            } else {
                                sumArrayWriteEnergy += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapRow * writeVoltageLTP * writeVoltageLTP;    // Selected WL (LTP phase)
                                sumArrayWriteEnergy += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapRow * writeVoltageLTP/2 * writeVoltageLTP/2 * (row_size - 1);  // Unselected WLs (LTP phase)
                                sumArrayWriteEnergy += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapCol * writeVoltageLTP/2 * writeVoltageLTP/2 * (col_size - cells);   // Unselected BLs (LTP phase)
                                sumArrayWriteEnergy += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapRow * writeVoltageLTD/2 * writeVoltageLTD/2 * (row_size - 1);    // Unselected WLs (LTD phase)
                                sumArrayWriteEnergy += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapCol * writeVoltageLTD/2 * writeVoltageLTD/2 * (col_size - cells); // Unselected BLs (LTD phase)
                            }


                            if (!static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->cmosAccess) { // Cross-point
                                for (int jj = 0; jj < col_size; jj++) { // Half-selected cells in the same row
                                    if (jj >= (col - array_col * operands_per_row) * (cells) && jj < (col - array_col * operands_per_row) * (cells+1)) { continue; } // Skip the selected cells
                                    sumArrayWriteEnergy += (writeVoltageLTP/2 * writeVoltageLTP/2 * static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][jj])->conductanceAtHalfVwLTP * maxLatencyLTP + writeVoltageLTD/2 * writeVoltageLTD/2 * static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][jj])->conductanceAtHalfVwLTD * maxLatencyLTD);
                                }
                                for (int kk = 0; kk < row_size; kk++) {    // Half-selected cells in other rows
                                    // Note that here is a bit inaccurate if using OpenMP, because the weight on other rows (threads) are also being updated
                                    if (kk == row - array_row * row_size) { continue; } // Skip the selected row
                                    for (int b = 0; b < cells; ++b) {
                                        sumArrayWriteEnergy += (writeVoltageLTP/2 * writeVoltageLTP/2 * static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->conductanceAtHalfVwLTP * maxLatencyLTP + writeVoltageLTD/2 * writeVoltageLTD/2 * static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->conductanceAtHalfVwLTD * maxLatencyLTD);
                                    }
                                }
                            }
                        }





                    }

#pragma opm critical    // Use critical here since NeuroSim class functions may update its member variables
                    {
                        if (AnalogNVM *temp = dynamic_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])) {  // Analog eNVM
                            int sumNumWritePulse = 0;
                            for (int j = 0; j < col_size; j++) {
                                sumNumWritePulse += abs(static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][j])->numPulse);    // Note that LTD has negative pulse number
                                sumNumWritePulse += abs(static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][j])->numPulse);    // Note that LTD has negative pulse number
                            }
                            subArrays[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->numWritePulse = sumNumWritePulse / col_size;
                            double writeVoltageSquareSumRow = 0;
                            if (static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->nonIdenticalPulse) { // Non-identical write pulse scheme
                                for (int j = 0; j < col_size; j++) {
                                    writeVoltageSquareSumRow += static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][j])->writeVoltageSquareSum;
                                    writeVoltageSquareSumRow += static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][j])->writeVoltageSquareSum;
                                }
                                if (sumNumWritePulse > 0) {	// Prevent division by 0
                                    subArrays[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell.writeVoltage = sqrt(writeVoltageSquareSumRow / sumNumWritePulse);	// RMS value of write voltage in a row
                                } else {
                                    subArrays[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell.writeVoltage = 0;
                                }
                            }
                        }
                        sumNeuroSimWriteEnergy += NeuroSimSubArrayWriteEnergy(subArrays[(int)(array_row * ceil((double)M/operands_per_row)) + array_col], numWriteOperationPerRow, 0);
                    }

                    numWriteOperation += numWriteOperationPerRow;
                    sumNeuroSimWriteEnergy += NeuroSimSubArrayWriteEnergy(subArrays[(int)(array_row * ceil((double)M/operands_per_row)) + array_col], numWriteOperationPerRow, 0);


                }



            }
        }
        writeEnergy += sumArrayWriteEnergy + sumNeuroSimWriteEnergy;
    }
    else{
        double sumArrayWriteEnergy = 0;   // Use a temporary variable here since OpenMP does not support reduction on class member
        double sumNeuroSimWriteEnergy = 0;   // Use a temporary variable here since OpenMP does not support reduction on class member
        double sumWriteLatencyAnalogNVM = 0;   // Use a temporary variable here since OpenMP does not support reduction on class member



        int array_size = ceil((double)K / row_size);
#pragma omp parallel for
        for (int array_row = 0; array_row < array_size; array_row++) {
            for (int array_col = 0; array_col < ceil((double)M / operands_per_row); array_col++) {

#pragma omp critical
                {
                for (int row = 0; row < K; row++) {
                    for (int col = 0; col < M; col++) {
                        unsigned char exponentA = get_exponent(A1[row * M + col]); //Only map weight1, the inputs are aligned in computation step, this may cause loss
                        if (row - array_row * row_size == 0
                            || exponentA > maxExponentA[array_row][col])
                            maxExponentA[array_row][col] = exponentA;
                    }
                }
                for (int row = array_row * row_size; row < K && row < row_size * (array_row + 1); row++) {
                    for (int col = array_col * operands_per_row; col < M && col < operands_per_row * (array_col + 1); col++) {
                        unsigned char exponentA = get_exponent(A1[row * M + col]);
                        char diff = maxExponentA[array_row][col] - exponentA;
                        //adjust the weight1
                        weight_mantissas[row * M + col] = weight_mantissas[row * M + col]
                                >> diff;
                    }
                }
                }

                for (int row = array_row * row_size;
                     row < K && row < row_size * (array_row + 1); row++) {


                    int numWriteOperationPerRow = 0; // # of write operations will be performed per row.
                    double numWriteOperation = 0;


                    double writeVoltageLTP = static_cast<eNVM*>(mp[0]->cell[0][0])->writeVoltageLTP;
                    double writeVoltageLTD = static_cast<eNVM*>(mp[0]->cell[0][0])->writeVoltageLTD;

                    for (int col = array_col * operands_per_row; col < M && col < operands_per_row * (array_col + 1);
                         col++) {

                        bool weightChangeBatch = false;
                        double maxLatencyLTP = 0;
                        double maxLatencyLTD = 0;
                        double cSum = 0,opSumNegative = 0,opSumPositive = 0;
                        if (weights[row * M + col] < 0) //if weight1 is negative
                        {
                            signal_matrix[row * M + col] = 1; //Useless
                            for (int b = 0; b < (cells); b++) {
                                if(((int)(array_row * ceil((double)M/operands_per_row)) + array_col) < n1) {
                                    int old_weight = mn[(int) (array_row * ceil((double) M / operands_per_row)) +
                                                        array_col]->getWeight(row - array_row * row_size,(col - array_col * operands_per_row) * (cells) + b);
                                    int new_weight = (weight_mantissas[col * K + row]
                                            >> (24 - (b + 1) * bits)) & ((int) pow(2, bits) - 1);
//                                    if (new_weight == old_weight);
//                                    else
                                        mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->WriteCell(
                                                row - array_row * row_size,
                                                (col - array_col * operands_per_row) * (cells) + b,
                                                new_weight-old_weight,
                                                new_weight,
                                                maxWeight, minWeight,
                                                true);
                                }
                                else{
                                    mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->WriteCell(
                                            row - array_row * row_size,
                                            (col - array_col * operands_per_row) * (cells) + b,
                                            (weight_mantissas[row * M + col]
                                                    >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1),
                                            (weight_mantissas[row * M + col]
                                                    >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1),
                                            maxWeight, minWeight,
                                            true);
                                }


                                cSum += currentDictionary[(weight_mantissas[col * K + row]
                                                >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1)];
                                cSum += currentDictionary[0];

                                opSumNegative += currentDictionary[(weight_mantissas[col * K + row]
                                                >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1)] * (1<<(bits * ((cells-1) - b)));

                                opSumPositive +=  currentDictionary[0] * (1<<(bits * ((cells-1) - b)));




                                weightChangeBatch = weightChangeBatch || static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse;

                                if(static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTP > maxLatencyLTP)
                                    maxLatencyLTP = static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTP;
                                if(static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTD > maxLatencyLTD)
                                    maxLatencyLTD = static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTD;

                            }
#pragma omp critical
                            {
                                currentSum[col * K + row] += cSum;

                                operandsSumNegative[col * K + row] += opSumNegative;
                                operandsSumPositive[col * K + row] += opSumPositive;
                            }

                            numWriteOperationPerRow += weightChangeBatch;

                            for (int b = 0; b < cells; ++b) {
                                static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTP = maxLatencyLTP;
                                static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTD = maxLatencyLTD;
                                if(weightChangeBatch){
                                    if (static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->nonIdenticalPulse) {	// Non-identical write pulse scheme
                                        if (static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse > 0) {	// LTP
                                            static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTP = sqrt(static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageSquareSum / static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse);	// RMS value of LTP write voltage
                                            static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTD = static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VinitLTD + 0.5 * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VstepLTD * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->maxNumLevelLTD;	// Use average voltage of LTD write voltage
                                        } else if (static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse < 0) {	// LTD
                                            static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTP = static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VinitLTP + 0.5 * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VstepLTP * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->maxNumLevelLTP;    // Use average voltage of LTP write voltage
                                            static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTD = sqrt(static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageSquareSum / (-1*static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse));    // RMS value of LTD write voltage
                                        } else {	// Half-selected during LTP and LTD phases
                                            static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTP = static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VinitLTP + 0.5 * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VstepLTP * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->maxNumLevelLTP;    // Use average voltage of LTP write voltage
                                            static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTD = static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VinitLTD + 0.5 * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VstepLTD * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->maxNumLevelLTD;    // Use average voltage of LTD write voltage
                                        }
                                    }
                                }
                                static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->WriteEnergyCalculation(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapCol);
                                sumArrayWriteEnergy += static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeEnergy;
                            }

                            sumWriteLatencyAnalogNVM += maxLatencyLTP + maxLatencyLTD;

                            if(array_col == 0)
                            wLatency[row] += maxLatencyLTP + maxLatencyLTD;

                            if (static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->nonIdenticalPulse) { // Non-identical write pulse scheme
                                writeVoltageLTP = static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->VinitLTP + 0.5 * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->VstepLTP * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->maxNumLevelLTP;    // Use average voltage of LTP write voltage
                                writeVoltageLTD = static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->VinitLTD + 0.5 * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->VstepLTD * static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->maxNumLevelLTD;    // Use average voltage of LTD write voltage
                            }
                            if (static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->cmosAccess) {  // 1T1R
                                // The energy on selected SLs is included in WriteCell()
                                sumArrayWriteEnergy += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireGateCapRow * tech->vdd * tech->vdd * 2;   // Selected WL (*2 means both LTP and LTD phases)
                                sumArrayWriteEnergy += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapRow * writeVoltageLTP * writeVoltageLTP;   // Selected BL (LTP phases)


                                //TODO: Verification
                                sumArrayWriteEnergy += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapCol * writeVoltageLTP * writeVoltageLTP * (col_size-cells);   // Unselected SLs (LTP phase)


                                // No LTD part because all unselected rows and columns are V=0
                            } else {
                                sumArrayWriteEnergy += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapRow * writeVoltageLTP * writeVoltageLTP;    // Selected WL (LTP phase)
                                sumArrayWriteEnergy += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapRow * writeVoltageLTP/2 * writeVoltageLTP/2 * (row_size - 1);  // Unselected WLs (LTP phase)
                                sumArrayWriteEnergy += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapCol * writeVoltageLTP/2 * writeVoltageLTP/2 * (col_size - cells);   // Unselected BLs (LTP phase)
                                sumArrayWriteEnergy += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapRow * writeVoltageLTD/2 * writeVoltageLTD/2 * (row_size - 1);    // Unselected WLs (LTD phase)
                                sumArrayWriteEnergy += mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapCol * writeVoltageLTD/2 * writeVoltageLTD/2 * (col_size - cells); // Unselected BLs (LTD phase)
                            }


                            if (!static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->cmosAccess) { // Cross-point
                                for (int jj = 0; jj < col_size; jj++) { // Half-selected cells in the same row
                                    if (jj >= (col - array_col * operands_per_row) * (cells) && jj < (col - array_col * operands_per_row) * (cells+1)) { continue; } // Skip the selected cells
                                    sumArrayWriteEnergy += (writeVoltageLTP/2 * writeVoltageLTP/2 * static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][jj])->conductanceAtHalfVwLTP * maxLatencyLTP + writeVoltageLTD/2 * writeVoltageLTD/2 * static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][jj])->conductanceAtHalfVwLTD * maxLatencyLTD);
                                }
                                for (int kk = 0; kk < row_size; kk++) {    // Half-selected cells in other rows
                                    // Note that here is a bit inaccurate if using OpenMP, because the weight on other rows (threads) are also being updated
                                    if (kk == row - array_row * row_size) { continue; } // Skip the selected row
                                    for (int b = 0; b < cells; ++b) {
                                        sumArrayWriteEnergy += (writeVoltageLTP/2 * writeVoltageLTP/2 * static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->conductanceAtHalfVwLTP * maxLatencyLTP + writeVoltageLTD/2 * writeVoltageLTD/2 * static_cast<eNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->conductanceAtHalfVwLTD * maxLatencyLTD);
                                    }
                                }
                            }

                        } else {
                            for (int b = 0; b < (cells); b++) {
                                if(((int)(array_row * ceil((double)M/operands_per_row)) + array_col) < n1) {
                                    int old_weight = mp[(int) (array_row * ceil((double) M / operands_per_row)) +
                                                        array_col]->getWeight(row - array_row * row_size,(col - array_col * operands_per_row) * (cells) + b);
                                    int new_weight = (weight_mantissas[col * K + row]
                                            >> (24 - (b + 1) * bits)) & ((int) pow(2, bits) - 1);
//                                    if (new_weight == old_weight);
//                                    else
                                        mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->WriteCell(
                                                row - array_row * row_size,
                                                (col - array_col * operands_per_row) * (cells) + b,
                                                new_weight - old_weight,
                                                new_weight,
                                                maxWeight, minWeight,
                                                true);
                                }
                                else{
                                    mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->WriteCell(
                                            row - array_row * row_size,
                                            (col - array_col * operands_per_row) * (cells) + b,
                                            (weight_mantissas[row * M + col]
                                                    >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1),
                                            (weight_mantissas[row * M + col]
                                                    >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1),
                                            maxWeight, minWeight,
                                            true);
                                }

                                cSum += currentDictionary[(weight_mantissas[col * K + row]
                                                >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1)];
                                cSum += currentDictionary[0];

                                opSumPositive += currentDictionary[(weight_mantissas[col * K + row]
                                                >> (24-(b+1)*bits)) & ((int)pow(2,bits)-1)] * (1<<(bits * ((cells-1) - b)));

                                opSumNegative +=  currentDictionary[0] * (1<<(bits * ((cells-1) - b)));






                                weightChangeBatch = weightChangeBatch || static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse;

                                if(static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTP > maxLatencyLTP)
                                    maxLatencyLTP = static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTP;
                                if(static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTD > maxLatencyLTD)
                                    maxLatencyLTD = static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTD;

                            }
#pragma omp critical
                            {
                                currentSum[col * K + row] += cSum;

                                operandsSumNegative[col * K + row] += opSumNegative;
                                operandsSumPositive[col * K + row] += opSumPositive;
                            }


                            numWriteOperationPerRow += weightChangeBatch;

                            for (int b = 0; b < cells; ++b) {
                                static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTP = maxLatencyLTP;
                                static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeLatencyLTD = maxLatencyLTD;
                                if(weightChangeBatch){
                                    if (static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->nonIdenticalPulse) {	// Non-identical write pulse scheme
                                        if (static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse > 0) {	// LTP
                                            static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTP = sqrt(static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageSquareSum / static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse);	// RMS value of LTP write voltage
                                            static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTD = static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VinitLTD + 0.5 * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VstepLTD * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->maxNumLevelLTD;	// Use average voltage of LTD write voltage
                                        } else if (static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse < 0) {	// LTD
                                            static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTP = static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VinitLTP + 0.5 * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VstepLTP * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->maxNumLevelLTP;    // Use average voltage of LTP write voltage
                                            static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTD = sqrt(static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageSquareSum / (-1*static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->numPulse));    // RMS value of LTD write voltage
                                        } else {	// Half-selected during LTP and LTD phases
                                            static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTP = static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VinitLTP + 0.5 * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VstepLTP * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->maxNumLevelLTP;    // Use average voltage of LTP write voltage
                                            static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeVoltageLTD = static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VinitLTD + 0.5 * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->VstepLTD * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->maxNumLevelLTD;    // Use average voltage of LTD write voltage
                                        }
                                    }
                                }
                                static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->WriteEnergyCalculation(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapCol);
                                sumArrayWriteEnergy += static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->writeEnergy;
                            }

                            sumWriteLatencyAnalogNVM += maxLatencyLTP + maxLatencyLTD;

                            if(array_col == 0)
                            wLatency[row] += maxLatencyLTP + maxLatencyLTD;

                            if (static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->nonIdenticalPulse) { // Non-identical write pulse scheme
                                writeVoltageLTP = static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->VinitLTP + 0.5 * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->VstepLTP * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->maxNumLevelLTP;    // Use average voltage of LTP write voltage
                                writeVoltageLTD = static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->VinitLTD + 0.5 * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->VstepLTD * static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->maxNumLevelLTD;    // Use average voltage of LTD write voltage
                            }
                            if (static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->cmosAccess) {  // 1T1R
                                // The energy on selected SLs is included in WriteCell()
                                sumArrayWriteEnergy += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireGateCapRow * tech->vdd * tech->vdd * 2;   // Selected WL (*2 means both LTP and LTD phases)
                                sumArrayWriteEnergy += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapRow * writeVoltageLTP * writeVoltageLTP;   // Selected BL (LTP phases)


                                //TODO: Verification
                                sumArrayWriteEnergy += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapCol * writeVoltageLTP * writeVoltageLTP * (col_size-cells);   // Unselected SLs (LTP phase)


                                // No LTD part because all unselected rows and columns are V=0
                            } else {
                                sumArrayWriteEnergy += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapRow * writeVoltageLTP * writeVoltageLTP;    // Selected WL (LTP phase)
                                sumArrayWriteEnergy += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapRow * writeVoltageLTP/2 * writeVoltageLTP/2 * (row_size - 1);  // Unselected WLs (LTP phase)
                                sumArrayWriteEnergy += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapCol * writeVoltageLTP/2 * writeVoltageLTP/2 * (col_size - cells);   // Unselected BLs (LTP phase)
                                sumArrayWriteEnergy += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapRow * writeVoltageLTD/2 * writeVoltageLTD/2 * (row_size - 1);    // Unselected WLs (LTD phase)
                                sumArrayWriteEnergy += mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->wireCapCol * writeVoltageLTD/2 * writeVoltageLTD/2 * (col_size - cells); // Unselected BLs (LTD phase)
                            }


                            if (!static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->cmosAccess) { // Cross-point
                                for (int jj = 0; jj < col_size; jj++) { // Half-selected cells in the same row
                                    if (jj >= (col - array_col * operands_per_row) * (cells) && jj < (col - array_col * operands_per_row) * (cells+1)) { continue; } // Skip the selected cells
                                    sumArrayWriteEnergy += (writeVoltageLTP/2 * writeVoltageLTP/2 * static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][jj])->conductanceAtHalfVwLTP * maxLatencyLTP + writeVoltageLTD/2 * writeVoltageLTD/2 * static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][jj])->conductanceAtHalfVwLTD * maxLatencyLTD);
                                }
                                for (int kk = 0; kk < row_size; kk++) {    // Half-selected cells in other rows
                                    // Note that here is a bit inaccurate if using OpenMP, because the weight on other rows (threads) are also being updated
                                    if (kk == row - array_row * row_size) { continue; } // Skip the selected row
                                    for (int b = 0; b < cells; ++b) {
                                        sumArrayWriteEnergy += (writeVoltageLTP/2 * writeVoltageLTP/2 * static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->conductanceAtHalfVwLTP * maxLatencyLTP + writeVoltageLTD/2 * writeVoltageLTD/2 * static_cast<eNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][(col - array_col * operands_per_row) * (cells) + b])->conductanceAtHalfVwLTD * maxLatencyLTD);
                                    }
                                }
                            }
                        }


#pragma opm critical    // Use critical here since NeuroSim class functions may update its member variables
                        {
                            if (AnalogNVM *temp = dynamic_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])) {  // Analog eNVM
                                int sumNumWritePulse = 0;
                                for (int j = 0; j < col_size; j++) {
                                    sumNumWritePulse += abs(static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][j])->numPulse);    // Note that LTD has negative pulse number
                                    sumNumWritePulse += abs(static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][j])->numPulse);    // Note that LTD has negative pulse number
                                }
                                subArrays[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->numWritePulse = sumNumWritePulse / col_size;
                                double writeVoltageSquareSumRow = 0;
                                if (static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[0][0])->nonIdenticalPulse) { // Non-identical write pulse scheme
                                    for (int j = 0; j < col_size; j++) {
                                        writeVoltageSquareSumRow += static_cast<AnalogNVM*>(mn[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][j])->writeVoltageSquareSum;
                                        writeVoltageSquareSumRow += static_cast<AnalogNVM*>(mp[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell[row - array_row * row_size][j])->writeVoltageSquareSum;
                                    }
                                    if (sumNumWritePulse > 0) {	// Prevent division by 0
                                        subArrays[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell.writeVoltage = sqrt(writeVoltageSquareSumRow / sumNumWritePulse);	// RMS value of write voltage in a row
                                    } else {
                                        subArrays[(int)(array_row * ceil((double)M/operands_per_row)) + array_col]->cell.writeVoltage = 0;
                                    }
                                }
                            }
                            sumNeuroSimWriteEnergy += NeuroSimSubArrayWriteEnergy(subArrays[(int)(array_row * ceil((double)M/operands_per_row)) + array_col], numWriteOperationPerRow, 0);
                        }

                        numWriteOperation += numWriteOperationPerRow;
                        sumNeuroSimWriteEnergy += NeuroSimSubArrayWriteEnergy(subArrays[(int)(array_row * ceil((double)M/operands_per_row)) + array_col], numWriteOperationPerRow, 0);

                    }
                }
            }
        }
        writeEnergy += sumArrayWriteEnergy + sumNeuroSimWriteEnergy;
    }



    A1 = NULL;
    delete[] currentDictionary;
    delete[] weight_mantissas;
}



template<class memoryType>
void Memristor_CU<memoryType>::memristor_mm(const float *input, uint64_t dims, float *C, float alpha, float beta)
{
    //printf("start mvm\n");
	unsigned int *input_mantissas;
	unsigned char *maxExponent;


	if(ot == WMI)
	{
		int operand_per_row = (int)(col_size/(cells));
		const unsigned int *B1 = reinterpret_cast<const unsigned int*>(input);
		unsigned short ak = ceil((double) K / row_size), am = ceil((double) M / operand_per_row);
		double sumArrayReadEnergy = 0;
		input_mantissas =	new unsigned int[dims * K];
		for(uint64_t i = 0; i < dims*K; i++)
		{
			if(input[i] == 0. || input[i] == -0.)
			{
				input_mantissas[i] = 0;
				continue;
			}
			if(get_exponent(B1[i]) == 255)
			{

                printf("B[i] = %f  dims = %ld   K = %d   Error: nan! \n", B1[i], dims, K);
                //for(int ii = 0 ; ii < K; ii++)
                //{
                //    for(int jj = 0 ; jj < dims; jj++)
                //        printf("%f  ",B1[i]);
                //    printf("\n");
                //}
				//exit(1);
                input_mantissas[i] = 0;
                continue;
			}
			input_mantissas[i] = get_mantissas(B1[i]);

		}
        //aligning the exponent
		maxExponent = new unsigned char[dims];
		for (uint64_t col = 0; col < dims; col++) {
			for (uint64_t row = 0; row < K; row++) {
				unsigned char exponentB = get_exponent(B1[row * dims + col]);
				if (row == 0|| exponentB > maxExponent[col])
					maxExponent[col] = exponentB;
			}
		}
		for (uint64_t row = 0; row < K; row++) {
			for (uint64_t col = 0; col < dims; col++) {
				unsigned char exponentB = get_exponent(B1[row * dims + col]);
				char diff = maxExponent[col] - exponentB;
				//adjust the input mantissas
				input_mantissas[row * dims + col] =	input_mantissas[row * dims + col] >> diff;
			}
		}

        for (uint64_t ii = 0; ii < am; ii++) {
            for (uint64_t jj = 0; jj < ak; jj++) {
                double readVoltage = static_cast<eNVM*>(mp[jj]->cell[0][0])->readVoltage;
                double readPulseWidth =	static_cast<eNVM*>(mp[jj]->cell[0][0])->readPulseWidth;

		#pragma omp parallel for reduction(+: sumArrayReadEnergy)
                for (uint64_t i_col = 0; i_col < dims; i_col++) {
                    for (uint64_t w_col = ii * operand_per_row; w_col < M && w_col < (ii + 1) * operand_per_row; w_col++) {
                        long long dot_product_sum = 0, dot_product_sum_negative = 0;
                        //for (int b = 0; b < (cells); b++) {
                            unsigned long long immediate_product = 0, immediate_product_negative = 0;
//								for (int n = 0; n < 24; n++) {
                            double Isum = 0; // weighted sum current
                            double IsumMax = 0; // Max weighted sum current
                            double IsumMin = 0;
                            double inputSum = 0;

                            double Isum_negative = 0;    // weighted sum current
                            double IsumMax_negative = 0; // Max weighted sum current
                            double IsumMin_negative = 0;
                            double inputSum_negative = 0;

                            for (uint64_t i_row = jj * row_size; i_row < K && i_row < (jj + 1) * row_size; i_row++) {
                                if (!(B1[i_row * dims + i_col] >> 31)) {
//											if ((input_mantissas[i_row * dims + i_col] >> n) & 1) {

                                    //Read out the current.  i_row,w_col
                                    Isum += input_mantissas[i_row * dims + i_col]*operandsSumPositive[w_col * K + i_row];//mp[jj * am + ii]->ReadCell(i_row - jj * row_size, (w_col - ii * operand_per_row) * (cells) + b); //row:i_col
                                    inputSum += input_mantissas[i_row * dims + i_col]*mp[jj * am + ii]->GetMinCellReadCurrent(i_row - jj * row_size, (w_col - ii * operand_per_row));
                                    Isum_negative += input_mantissas[i_row * dims + i_col]* operandsSumNegative[w_col * K + i_row] ;//mn[jj * am + ii]->ReadCell(i_row - jj * row_size, (w_col - ii * operand_per_row) * (cells) + b); //row:i_col
                                    inputSum_negative += input_mantissas[i_row * dims + i_col]*mn[jj * am + ii]->GetMinCellReadCurrent(i_row - jj * row_size, (w_col - ii * operand_per_row));


//											}
                                }

                                IsumMax += mp[jj * am + ii]->GetMaxCellReadCurrent(i_row - jj * row_size, (w_col - ii * operand_per_row));
                                IsumMin += mp[jj * am + ii]->GetMinCellReadCurrent(i_row - jj * row_size, (w_col - ii * operand_per_row));
                            }

                            long long outputDigits, outputDigits_negative;

                            if ((jj + 1) * row_size > K)
                            {
                                outputDigits = CurrentToDigits(Isum - inputSum, (IsumMax - IsumMin) * row_size / (K - jj * row_size),bits,row_size); //need to be further checked
                                outputDigits_negative = CurrentToDigits(Isum_negative - inputSum_negative, (IsumMax - IsumMin) * row_size / (K - jj * row_size),bits,row_size); //need to be further checked

                            }
                            else
                            {
                                outputDigits = CurrentToDigits(Isum - inputSum, (IsumMax - IsumMin),bits,row_size);
                                outputDigits_negative = CurrentToDigits(Isum_negative - inputSum_negative, (IsumMax - IsumMin),bits,row_size); //need to be further checked
                            }

                            immediate_product += outputDigits;// << n;
                            immediate_product_negative += outputDigits_negative;// << n;


                            Isum = 0;
                            inputSum = 0;
                            Isum_negative = 0;
                            inputSum_negative = 0;
                            for (uint64_t i_row = jj * row_size; i_row < K && i_row < (jj + 1) * row_size; i_row++) {
                                if ((B1[i_row * dims + i_col] >> 31)) {
//											if ((input_mantissas[i_row * dims + i_col] >> n) & 1) {
                                    Isum += input_mantissas[i_row * dims + i_col]*operandsSumPositive[w_col * K + i_row];//mp[jj * am + ii]->ReadCell(i_row - jj * row_size, (w_col - ii * operand_per_row) * (cells) + b); //row:i_col
                                    inputSum += input_mantissas[i_row * dims + i_col]*mp[jj * am + ii]->GetMinCellReadCurrent(i_row - jj * row_size, (w_col - ii * operand_per_row));
                                    Isum_negative += input_mantissas[i_row * dims + i_col]* operandsSumNegative[w_col * K + i_row] ;//mn[jj * am + ii]->ReadCell(i_row - jj * row_size, (w_col - ii * operand_per_row) * (cells) + b); //row:i_col
                                    inputSum_negative += input_mantissas[i_row * dims + i_col]*mn[jj * am + ii]->GetMinCellReadCurrent(i_row - jj * row_size, (w_col - ii * operand_per_row));



//											}
                                }
                            }

                            if ((jj + 1) * row_size > K)
                            {
                                outputDigits = CurrentToDigits(Isum - inputSum, (IsumMax - IsumMin) * row_size / (K - jj * row_size),bits,row_size); //need to be further checked
                                outputDigits_negative = CurrentToDigits(Isum_negative - inputSum_negative, (IsumMax - IsumMin) * row_size / (K - jj * row_size),bits,row_size); //need to be further checked

                            }
                            else
                            {
                                outputDigits = CurrentToDigits(Isum - inputSum, (IsumMax - IsumMin),bits,row_size);
                                outputDigits_negative = CurrentToDigits(Isum_negative - inputSum_negative, (IsumMax - IsumMin),bits,row_size); //need to be further checked
                            }

                            immediate_product += outputDigits_negative;// << n;
                            immediate_product_negative += outputDigits;// << n;


//								}
//                                if(i_col == 0 && w_col == 0 && ii == 0 && jj == 0 && b == 0)
//                                    printf("%X", immediate_product);

                            //dot_product_sum += immediate_product<<(bits * ((cells-1) - b));


                            //dot_product_sum_negative += immediate_product_negative<<(bits * ((cells-1) - b));


                            dot_product_sum += immediate_product;


                            dot_product_sum_negative += immediate_product_negative;
//                        }

//                        printf("%X %X %d %d %d ",dot_product_sum,dot_product_sum_negative,maxExponent[i_col] , maxExponentA[jj][w_col],(254+(cells*bits-1)+23));

                        dot_product_sum -= dot_product_sum_negative;



                        if (jj == 0) {
                            //Convert to floating-point format.
                            C[w_col * dims + i_col] = beta * C[w_col * dims + i_col] + alpha * dot_product_sum
                                                                                       * pow(2, maxExponent[i_col] + maxExponentA[jj][w_col] - (254+(cells*bits-1)+23)/*300*/);
                        } else {
                            C[w_col * dims + i_col] += alpha * dot_product_sum * pow(2, maxExponent[i_col] + maxExponentA[jj][w_col] - (254+(cells*bits-1)+23));
                        }
                        unsigned int *CC = reinterpret_cast<unsigned int*>(C);
			            if(get_exponent(CC[w_col *dims + i_col]) == 255){
                            printf("nan\n");
                            float *EE = new float[dims*M];
                        	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, dims, K, 1., A, K, input,
	                            dims, 0., EE, dims);
                            printf("maxExponent = %d maxExponentA=%d  dot_product_sum = %llu  M = %ld K = %ld N = %ld   result = %.4e  hex = %x  cb = %d  exponent = %d  result=%llu ", maxExponent[i_col], maxExponentA[jj][w_col], dot_product_sum, M, K, dims, 
                                   EE[w_col*dims+i_col], *reinterpret_cast<unsigned int*>(&EE[w_col*dims+i_col]),cells*bits-1,maxExponent[i_col] + maxExponentA[jj][w_col] - (254+(cells*bits-1)+23), *reinterpret_cast<uint32_t*>(&(C[w_col * dims + i_col]) ));
                            std::ofstream out1("weight"),out2("input");
                            ////out1.open();
                            ////out2.open();
                            for(uint64_t i = 0; i < M*K; i++){
                                out1<<A[i]<<",";
                            }
                            for(uint64_t i = 0; i < dims*K; i++){
                                out2<<input[i]<<",";
                            }
                            out1.close();
                            out2.close();
                            exit(3);
                        }
                    }
                }

            }
        }
	}
	else{
		int operand_per_row = (int)(col_size/(cells));
		const unsigned int *B1 = reinterpret_cast<const unsigned int*>(input);
		unsigned short ak = ceil((double) K / row_size), am = ceil((double) M / operand_per_row);
		double sumArrayReadEnergy = 0;
		input_mantissas =	new unsigned int[dims * K];
		for(uint64_t i = 0; i < dims*K; i++)
		{
			if(input[i] == 0. || input[i] == -0.)
			{
				input_mantissas[i] = 0;
				continue;
			}
			if(get_exponent(B1[i]) == 255)
			{
				printf("Error: nan! \n");
				exit(1);
			}
			input_mantissas[i] = get_mantissas(B1[i]);
		}
		maxExponent = new unsigned char[dims];
		for (uint64_t col = 0; col < K; col++) {
			for (uint64_t row = 0; row < dims; row++) {
				unsigned char exponentB = get_exponent(B1[row * K + col]);
				if (col == 0|| exponentB > maxExponent[row])
					maxExponent[row] = exponentB;
			}
		}
		for (uint64_t row = 0; row < dims; row++) {
			for (uint64_t col = 0; col < K; col++) {
				unsigned char exponentB = get_exponent(B1[row * K + col]);
				char diff = maxExponent[row] - exponentB;
				//adjust the input mantissas
				input_mantissas[row * K + col] = input_mantissas[row * K + col] >> diff;


			}
		}
        for (uint64_t ii = 0; ii < am; ii++) {
            for (uint64_t jj = 0; jj < ak; jj++) {
                double readVoltage = static_cast<eNVM*>(mp[jj]->cell[0][0])->readVoltage;
                double readPulseWidth =	static_cast<eNVM*>(mp[jj]->cell[0][0])->readPulseWidth;
                //aligning the exponent
		#pragma omp parallel for reduction(+: sumArrayReadEnergy)
                for (uint64_t i_row = 0; i_row < dims; i_row++) {
                    for (uint64_t w_col = ii * operand_per_row; w_col < M && w_col < (ii + 1) * operand_per_row; w_col++) {
                        long long dot_product_sum = 0, dot_product_sum_negative = 0;
//                        for (int n = 0; n < 24; n++) {

//                            for (int b = 0; b < (cells); b++) {
                                unsigned long long immediate_product = 0, immediate_product_negative = 0;

                                double Isum = 0; // weighted sum current
                                double IsumMax = 0; // Max weighted sum current
                                double IsumMin = 0;
                                double inputSum = 0;

                                double Isum_negative = 0;    // weighted sum current
                                double IsumMax_negative = 0; // Max weighted sum current
                                double IsumMin_negative = 0;
                                double inputSum_negative = 0;

                                for (uint64_t i_col = jj * row_size; i_col < K && i_col < (jj + 1) * row_size; i_col++) {
                                    if (!(B1[i_row * K + i_col] >> 31)) {
//                                        if ((input_mantissas[i_row * K + i_col] >> n) & 1) { i_col,w_col
                                            Isum += input_mantissas[i_row * K + i_col]*operandsSumPositive[w_col * K + i_col]; //row:i_col
                                            inputSum += input_mantissas[i_row * K + i_col]*mp[jj * am + ii]->GetMinCellReadCurrent(i_col - jj * row_size, (w_col - ii * operand_per_row) * (cells));

                                            Isum_negative += input_mantissas[i_row * K + i_col]*operandsSumNegative[w_col * K + i_col]; //row:i_col
                                            inputSum_negative += input_mantissas[i_row * K + i_col]*mn[jj * am + ii]->GetMinCellReadCurrent(i_col - jj * row_size, (w_col - ii * operand_per_row) * (cells));



//                                        }
                                    }

                                    IsumMax += mp[jj * am + ii]->GetMaxCellReadCurrent(i_col - jj * row_size, (w_col - ii * operand_per_row));
                                    IsumMin += mp[jj * am + ii]->GetMinCellReadCurrent(i_col - jj * row_size, (w_col - ii * operand_per_row));
                                }
                                sumArrayReadEnergy += Isum * readVoltage * readPulseWidth;
                                long long outputDigits, outputDigits_negative;

                                if ((jj + 1) * row_size > K)
                                {
                                    outputDigits = CurrentToDigits(Isum - inputSum, (IsumMax - IsumMin) * row_size / (K - jj * row_size),bits,row_size); //need to be further checked
                                    outputDigits_negative = CurrentToDigits(Isum_negative - inputSum_negative, (IsumMax - IsumMin) * row_size / (K - jj * row_size),bits,row_size); //need to be further checked

                                }
                                else
                                {
                                    outputDigits = CurrentToDigits(Isum - inputSum, (IsumMax - IsumMin),bits,row_size);
                                    outputDigits_negative = CurrentToDigits(Isum_negative - inputSum_negative, (IsumMax - IsumMin),bits,row_size); //need to be further checked
                                }

                                immediate_product += outputDigits;// << (bits * ((cells-1) - b));
                                immediate_product_negative += outputDigits_negative;// << (bits * ((cells-1) - b));

//
                                Isum = 0;
                                inputSum = 0;
                                Isum_negative = 0;
                                inputSum_negative = 0;
                                for (uint64_t i_col = jj * row_size; i_col < K && i_col < (jj + 1) * row_size; i_col++) {
                                    if ((B1[i_row * K + i_col] >> 31)) {
//                                        if ((input_mantissas[i_row * K + i_col] >> n) & 1) {
                                            Isum += input_mantissas[i_row * K + i_col]*operandsSumPositive[w_col * K + i_col]; //row:i_col
                                            inputSum += input_mantissas[i_row * K + i_col]*mp[jj * am + ii]->GetMinCellReadCurrent(i_col - jj * row_size, (w_col - ii * operand_per_row) * (cells));

                                            Isum_negative += input_mantissas[i_row * K + i_col]*operandsSumNegative[w_col * K + i_col]; //row:i_col
                                            inputSum_negative += input_mantissas[i_row * K + i_col]*mn[jj * am + ii]->GetMinCellReadCurrent(i_col - jj * row_size, (w_col - ii * operand_per_row) * (cells));



//                                        }
                                    }
                                }

                                if ((jj + 1) * row_size > K)
                                {
                                    outputDigits = CurrentToDigits(Isum - inputSum, (IsumMax - IsumMin) * row_size / (K - jj * row_size),bits,row_size); //need to be further checked
                                    outputDigits_negative = CurrentToDigits(Isum_negative - inputSum_negative, (IsumMax - IsumMin) * row_size / (K - jj * row_size),bits,row_size); //need to be further checked

                                }
                                else
                                {
                                    outputDigits = CurrentToDigits(Isum - inputSum, (IsumMax - IsumMin),bits,row_size);
                                    outputDigits_negative = CurrentToDigits(Isum_negative - inputSum_negative, (IsumMax - IsumMin),bits,row_size); //need to be further checked
                                }

                                immediate_product += outputDigits_negative;// << (bits * ((cells-1) - b));
                                immediate_product_negative += outputDigits;// << (bits * ((cells-1) - b));



                                dot_product_sum += immediate_product;
                                dot_product_sum_negative += immediate_product_negative;
//                            }

//                        }


                        dot_product_sum -= dot_product_sum_negative;
                        if (jj == 0) {
                            C[i_row * M + w_col] = beta * C[i_row * M + w_col] + alpha * dot_product_sum
                                                                                 * pow(2, maxExponent[i_row] + maxExponentA[jj][w_col] - (254+(cells*bits-1)+23));
                        } else {
                            C[i_row * M + w_col] += alpha * dot_product_sum * pow(2, maxExponent[i_row] + maxExponentA[jj][w_col] - (254+(cells*bits-1)+23));
                        }

                    }
                }

            }
        }
	}




	delete[] input_mantissas;
//	delete[] input_signal;
	delete[] maxExponent;
	input_mantissas = NULL;
//	input_signal = NULL;
	maxExponent = NULL;
}


template<class memoryType>
void Memristor_CU<memoryType>::Energy(const float *input, uint64_t dims){
    int operand_per_row = (int)(col_size/(cells));
    const unsigned int *B1 = reinterpret_cast<const unsigned int*>(input);
    unsigned short ak = ceil((double) K / row_size), am = ceil((double) M / operand_per_row);
    unsigned char *maxExponent;
    unsigned int *input_mantissas;
    input_mantissas =	new unsigned int[dims * K];
    unsigned short n = ak * am;
    unsigned char *count = new unsigned char[dims * K];

    readEnergy = 0;
    int total_num = 0;

    if(ot == WMI){
        for(uint64_t i = 0; i < dims*K; i++)
        {
            count[i] = 0;
            if(input[i] == 0. || input[i] == -0.)
            {
                input_mantissas[i] = 0;
                continue;
            }
            if(get_exponent(B1[i]) == 255)
            {
                printf("Error: nan! \n");
                exit(1);
            }
            input_mantissas[i] = get_mantissas(B1[i]);
        }
        //	input_signal[i] = (input[i] >= 0);
        maxExponent = new unsigned char[dims];
        for (uint64_t col = 0; col < dims; col++) {
            for (uint64_t row = 0; row < K; row++) {
                unsigned char exponentB = get_exponent(B1[row * dims + col]);
                if (row == 0|| exponentB > maxExponent[col])
                    maxExponent[col] = exponentB;
            }
        }
        for (uint64_t row = 0; row < K; row++) {
            for (uint64_t col = 0; col < dims; col++) {
                unsigned char exponentB = get_exponent(B1[row * dims + col]);
                char diff = maxExponent[col] - exponentB;
                //adjust the input mantissas
                input_mantissas[row * dims + col] =	input_mantissas[row * dims + col] >> diff;
            }
        }
    }else{
        for(uint64_t i = 0; i < dims*K; i++)
        {
            count[i] = 0;
            if(input[i] == 0. || input[i] == -0.)
            {
                input_mantissas[i] = 0;
                continue;
            }
            if(get_exponent(B1[i]) == 255)
            {
                printf("Error: nan! \n");
                exit(1);
            }
            input_mantissas[i] = get_mantissas(B1[i]);
        }
        maxExponent = new unsigned char[dims];
        for (uint64_t col = 0; col < K; col++) {
            for (uint64_t row = 0; row < dims; row++) {
                unsigned char exponentB = get_exponent(B1[row * K + col]);
                if (col == 0|| exponentB > maxExponent[row])
                    maxExponent[row] = exponentB;
            }
        }
        for (uint64_t row = 0; row < dims; row++) {
            for (uint64_t col = 0; col < K; col++) {
                unsigned char exponentB = get_exponent(B1[row * K + col]);
                char diff = maxExponent[row] - exponentB;
                //adjust the input mantissas
                input_mantissas[row * K + col] = input_mantissas[row * K + col] >> diff;


            }
        }
    }

    for (uint64_t i = 0; i < dims*K; ++i) {
        for(uint64_t j = 0; j < 24; j++)
        {
            if(input_mantissas[i]>>j&1)
                count[i]++;
        }
    }



    for (uint64_t ii = 0; ii < am; ii++) {
        for (uint64_t jj = 0; jj < ak; jj++) {
            double sumArrayReadEnergy = 0;
            double sumArrayReadEnergyN = 0;
            int i = jj * am + ii;
            double readVoltage = 0;
            double readPulseWidth = 0;
            readVoltage = static_cast<eNVM*>(mp[jj * am + ii]->cell[0][0])->readVoltage;
            readPulseWidth = static_cast<eNVM*>(mp[jj * am + ii]->cell[0][0])->readPulseWidth;


            if(ot == WMI){
//#pragma omp parallel for reduction(+: sumArrayReadEnergy, sumArrayReadEnergyN)
                for (uint64_t i_col = 0; i_col < dims; i_col++) {
                    for (uint64_t w_col = ii * operand_per_row; w_col < M && w_col < (ii + 1) * operand_per_row; w_col++) {
                        double Isum = 0; // weighted sum current
                        double Isum_negative = 0;    // weighted sum current

                        for (uint64_t i_row = jj * row_size; i_row < K && i_row < (jj + 1) * row_size; i_row++) {

                            Isum += count[i_row * dims + i_col]*currentSum[w_col*K+i_row];//mp[jj * am + ii]->ReadCell(i_row - jj * row_size, (w_col - ii * operand_per_row) * (cells) + b); //row:i_col
                            sumArrayReadEnergy += (cells)*count[i_row * dims + i_col]*mp[jj * am + ii]->wireCapRow * readVoltage * readVoltage; // Selected BLs (1T1R) or Selected WLs (cross-point)
                            sumArrayReadEnergyN += (cells)*count[i_row * dims + i_col]*mn[jj * am + ii]->wireCapRow * readVoltage * readVoltage; // Selected BLs (1T1R) or Selected WLs (cross-point)

                        }

                        //printf("Isum = %.4e  I_ReadCell = %.4e sumEnergy = %.4e \n",Isum,mp[0]->ReadCell(0,0),Isum*readVoltage*readPulseWidth);
                        sumArrayReadEnergy += Isum * readVoltage * readPulseWidth;
                        sumArrayReadEnergyN += Isum_negative * readVoltage * readPulseWidth;
                    }
                }
            }else{
//#pragma omp parallel for reduction(+: sumArrayReadEnergy, sumArrayReadEnergyN)
                for (uint64_t i_row = 0; i_row < dims; i_row++) {
                    for (uint64_t w_col = ii * operand_per_row; w_col < M && w_col < (ii + 1) * operand_per_row; w_col++) {
                        double Isum = 0; // weighted sum current
                        double Isum_negative = 0;    // weighted sum current

                        for (uint64_t i_col = jj * row_size; i_col < K && i_col < (jj + 1) * row_size; i_col++) {

                            Isum += count[i_row * K + i_col]*currentSum[i_col * M + w_col];//mp[jj * am + ii]->ReadCell(i_row - jj * row_size, (w_col - ii * operand_per_row) * (cells) + b); //row:i_col
                            sumArrayReadEnergy += (cells)*count[i_row * K + i_col]*mp[jj * am + ii]->wireCapRow * readVoltage * readVoltage; // Selected BLs (1T1R) or Selected WLs (cross-point)

                            sumArrayReadEnergyN += (cells)*count[i_row * K + i_col]*mn[jj * am + ii]->wireCapRow * readVoltage * readVoltage; // Selected BLs (1T1R) or Selected WLs (cross-point)

                        }
                        sumArrayReadEnergy += Isum * readVoltage * readPulseWidth;
                        sumArrayReadEnergyN += Isum_negative * readVoltage * readPulseWidth;
                    }
                }

            }
            if(ii == 0 && jj == 0)
            printf("supply:%f    WLs open = %f  \n",sumArrayReadEnergy,24*dims*numColMuxed*mp[jj*am+ii]->wireGateCapRow*tech->vdd*tech->vdd*(K-jj*row_size));
            if (static_cast<eNVM*>(mp[jj * am + ii]->cell[0][0])->cmosAccess) {  // 1T1R
                if ((jj + 1) * row_size > K)
                {
                    sumArrayReadEnergy += 24*dims * numColMuxed * mp[jj * am + ii]->wireGateCapRow * tech->vdd * tech->vdd * (K - jj * row_size); // All WLs open
                    sumArrayReadEnergyN += 24*dims * numColMuxed * mn[jj * am + ii]->wireGateCapRow * tech->vdd * tech->vdd * (K - jj * row_size); // All WLs open
                }
                else
                {
                    sumArrayReadEnergyN += 24*dims * numColMuxed * mn[jj * am + ii]->wireGateCapRow * tech->vdd * tech->vdd * row_size; // All WLs open
                    sumArrayReadEnergy += 24*dims * numColMuxed * mp[jj * am + ii]->wireGateCapRow * tech->vdd * tech->vdd * row_size; // All WLs open
                }

            }
            mp[jj * am + ii]->readEnergy += sumArrayReadEnergy;
            mn[jj * am + ii]->readEnergy += sumArrayReadEnergyN;
            double SubarrayDynamicEnergy = 0;
//                for (int b = 0; b < (cells); b++) {
            int numBatchReadSynapse = (int)ceil((double)col_size/numColMuxed);

            for (uint64_t j = 0; j < dims; ++j) {
                for (uint64_t w_col = 0; w_col < operand_per_row*cells; w_col+=numBatchReadSynapse) {
                    SubarrayDynamicEnergy += 4*NeuroSimSubArrayReadEnergy(subArrays[jj * am + ii]); //2 times input  2 times array
                    SubarrayDynamicEnergy += 4*NeuroSimNeuronReadEnergy(subArrays[jj * am + ii], adder[jj * am + ii], mux[jj * am + ii], muxDecoder[jj * am + ii], dff[jj * am + ii], subtractor[jj * am + ii]);
                    //printf("readEnergy = %.4e %.4e \n",NeuroSimSubArrayReadEnergy(subArrays[jj*am+ii])+NeuroSimNeuronReadEnergy(subArrays[jj*am+ii],adder[jj*am+ii],mux[jj*am+ii],muxDecoder[jj*am+ii],dff[jj*am+ii],subtractor[jj*am+ii]),SubarrayDynamicEnergy);
//                        if(ii == 0 && jj == 0)
//                        printf("%.4e %.4e \n",NeuroSimSubArrayReadEnergy(subArrays[jj * am + ii]),NeuroSimNeuronReadEnergy(subArrays[jj * am + ii], adder[jj * am + ii], mux[jj * am + ii], muxDecoder[jj * am + ii], dff[jj * am + ii], subtractor[jj * am + ii]));

                }
            }
            printf("SubarrayDynamicEnergy = %.4e, readEnergy = %.4e  numBatchReadSynapse = %d \n",SubarrayDynamicEnergy,mp[jj*am+ii]->readEnergy,numBatchReadSynapse);
//                }
            readEnergy += SubarrayDynamicEnergy + mp[jj * am + ii]->readEnergy + mn[jj * am + ii]->readEnergy;

        }
    }

    delete[] input_mantissas;
    delete[] maxExponent;
    delete[] count;
    input_mantissas = NULL;
    maxExponent = NULL;
    count = NULL;


//    printf("read latency=%.4e s\n", readLatency);
//    printf("read energy=%.4e s\n", readEnergy);

}

template<class memoryType>
void Memristor_CU<memoryType>::ReadLatency()
{
    //The latency of vector multiply with Matrix

    int operand_per_row = (int)(col_size/(cells));
    unsigned short ak = ceil((double) K / row_size), am = ceil((double) M / operand_per_row);
    unsigned short n = ak * am;


    readLatency = 0;
    double SubarrayReadLatency = 0;
    //numColMuxed = 2;
    int numBatchReadSynapse = (int)ceil((double)operand_per_row*cells/numColMuxed);

    for (int w_col = 0; w_col < operand_per_row*cells; w_col+=numBatchReadSynapse) {
        SubarrayReadLatency += 2*NeuroSimSubArrayReadLatency(subArrays[0]);
        SubarrayReadLatency += 2*NeuroSimNeuronReadLatency(subArrays[0], adder[0], mux[0], muxDecoder[0], dff[0], subtractor[0]);
        //printf("%.7e %.7e \n",2*NeuroSimSubArrayReadLatency(subArrays[0]),2*NeuroSimNeuronReadLatency(subArrays[0],adder[0],mux[0],muxDecoder[0],dff[0],subtractor[0]));
    }
    readLatency += SubarrayReadLatency;
}




