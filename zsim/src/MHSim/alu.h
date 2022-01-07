/*
 * @Author: jhxu
 * @LastEditors: jhxu
 * @LastEditTime: 2022-01-05 17:23:53
 * @FilePath: /src/MHSim/alu.h
 */


#ifndef ZSIM_ALU_H
#define ZSIM_ALU_H
#include "buffer.h"
#include "../memory_hierarchy.h"
class ALU{
public:
    ALU()
    {}
    ALU(Buffer *_buffer):buffer(_buffer){}
    ~ALU(){}
    void preprocess(Address src, uint64_t &memCycle);
    void process(){}
private:
    Buffer *buffer;
};

#endif //ZSIM_ALU_H
