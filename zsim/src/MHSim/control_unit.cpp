/*
 * @Author: jhxu
 * @LastEditTime: 2022-01-05 17:22:48
 * @LastEditors: jhxu
 * @FilePath: /src/MHSim/control_unit.cpp
 */

#include "control_unit.h"

uint64_t Control_Unit::execute(Instruction ins, uint64_t& cycles){

    MESIState dummyState;
    switch(ins.opid){
    case load:
        {
        dummyState = MESIState::E;
        MemReq req = {ins.src, GETS, 0, &dummyState, cycles, &filterLock, dummyState, 0, 0};
        cycles = buffer->access(req);
        }
        break;
    case alu:{}
        break;
    case store:
        {
        dummyState = MESIState::M;
        MemReq req = {ins.dst, PUTX, 0, &dummyState, cycles, &filterLock, dummyState, 0, 0};
        cycles = buffer->access(req);
        }
        break;
    case map:
        {}
        break;
    case mvm:
        {}
        break;
    case OPID::copy:
        {}
        break;
    case alui:
        {}
        break;
    }
    cycles += 2;
    return 0;
}
