/*
 * @Author: jhxu
 * @LastEditors: jhxu
 * @LastEditTime: 2022-01-05 17:22:20
 * @FilePath: /src/MHSim/selector.cpp
 */

#include "selector.h"
#include <stdio.h>
bool Selector::trade_off(uint32_t M, uint32_t K, uint32_t N, uint64_t frequency, uint32_t bits){
    float cpu_time = 0, MBA_time = 0;
    uint32_t mul_time = 0, accu_time = 0;
    mul_time = N * K * M;
//    accu_time = (K-1) * N * M;
// FMUL and FADD are both 5 cycles.
    mul_time *= 5;
    accu_time *= 5;
    cpu_time = 1.0*(mul_time + accu_time)/frequency;
    MBA_time = N * 0.0000120811 /120 * bits / 24;
    //printf("mul_time = %lu accu_time = %lu frequency = %llu  cpu_time = %.4e  MBA_time = %.4e \n", mul_time, accu_time, frequency, cpu_time, MBA_time);
    if(MBA_time < cpu_time && K != 1)
        return true;
    else
        return false;
}
