/*
 * @Author: jhxu
 * @LastEditors: jhxu
 * @LastEditTime: 2022-01-05 17:23:44
 * @FilePath: /src/MHSim/buffer.cpp
 */


#include "buffer.h"

uint64_t Buffer::access(MemReq& req){

    if(req.type == PUTX)
        return mems->access(req);

    for (int i = 0; i < numLines; ++i) {
        if((bufferEntry[i].address <= req.lineAddr && bufferEntry[i].address + lineSize > req.lineAddr)){
            return req.cycle;
        }
    }
    int x = 0;
    uint64_t minCycle = wait(x);
    if(minCycle > req.cycle){
        req.cycle = minCycle;
    }

    uint64_t respCycle = mems->access(req);

    setFlags(req.lineAddr,uselessIdx,bufferEntry[uselessIdx].type,respCycle);
    assert(respCycle >= req.cycle);
/*
if ((req.type == GETS || req.type == GETX) && eventRecorders[req.srcId]) {
    Address addr = req.lineAddr<<lineBits;
    MemAccReqEvent* memEv = new (eventRecorders[req.srcId]->alloc<MemAccReqEvent>()) MemAccReqEvent(nullptr, false, addr);
    TimingRecord tr = {addr, req.cycle, respCycle, req.type, memEv, memEv};
    eventRecorders[req.srcId]->pushRecord(tr);
}
*/
    return respCycle;
}

uint64_t Buffer::wait(int &idx){
    uint64_t minCycle = 0;
    for (int i = 0; i < numLines; ++i) {
        if(bufferEntry[i].isUseless){
            if(minCycle == 0){
                uselessIdx = i;
                minCycle = bufferEntry[i].memCycle;
            }
            else{
                if(bufferEntry[i].memCycle<minCycle){
                    uselessIdx = i;
                    minCycle = bufferEntry[i].memCycle;
                }

            }
        }
    }
    idx = uselessIdx;
    return minCycle;
}