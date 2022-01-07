/*
 * @Author: jhxu
 * @LastEditors: jhxu
 * @LastEditTime: 2022-01-05 17:22:28
 * @FilePath: /src/MHSim/selector.h
 */
#ifndef ZSIM_SELECTOR_H
#define ZSIM_SELECTOR_H
#include <stdint.h>

class Selector{
public:
    Selector(){}
    bool trade_off(uint32_t M, uint32_t K, uint32_t N, uint64_t frequency, uint32_t input_bits = 24);
};


#endif
