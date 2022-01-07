/*
 * @Author: jhxu
 * @LastEditors: jhxu
 * @LastEditTime: 2022-01-05 17:23:59
 * @FilePath: /src/MHSim/alu.cpp
 */

#include "alu.h"

void ALU::preprocess(Address src, uint64_t &memCycle){
    //Subtract, shift, ReLU, and max. We assume it is done in 1 cycles.
    memCycle += 1;
}