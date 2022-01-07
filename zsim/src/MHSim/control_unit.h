/*
 * @Author: jhxu
 * @LastEditTime: 2022-01-05 17:22:52
 * @LastEditors: jhxu
 * @FilePath: /src/MHSim/control_unit.h
 */


#ifndef ZSIM_CONTROLL_UNIT_H
#define ZSIM_CONTROLL_UNIT_H

#include "alu.h"
enum OPID{
    map = 0,
    mvm = 1,
    alu = 2,
    alui = 3,
    copy = 4,
    load = 5,
    store = 6
};

struct Instruction{
    OPID opid;
    Address src;
    Address dst;
    uint8_t reserve;
};

class Control_Unit{
public:
    Control_Unit()
    {
        futex_init(&filterLock);
    }

    ~Control_Unit(){}

    Control_Unit(Buffer *_buffer, ALU *_alu):buffer(_buffer),Alu(_alu)
    {
        futex_init(&filterLock);
    }

    uint64_t execute(Instruction ins, uint64_t& cycles);
private:
    Buffer *buffer;
    ALU *Alu;

    lock_t filterLock;
};

#endif //ZSIM_CONTROLL_UNIT_H
