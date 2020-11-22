//
// Created by jiahong on 20-1-25.
//

#ifndef ZSIM_BUFFER_H
#define ZSIM_BUFFER_H

#include "../galloc.h"
#include "coherence.h"
#include "neuro/Memristor_CU.h"
#include "iostream"
#include "../memory_hierarchy.h"
using namespace std;

class Buffer{
private:
    struct BufferEntry{
        bool isUseless;
        Address address;
        OperationTypes type;
        uint64_t memCycle;

        void clear(){isUseless = true; address = 0;memCycle = 0;}
    };
    BufferEntry *bufferEntry;
    uint32_t numLines;
    uint32_t cursor;
    uint64_t latency;
    uint32_t lineSize;
    uint32_t uselessIdx;
    MemObject *mems;
public:
    Buffer(uint32_t numLines, uint32_t lineSize, MemObject *_mems/*Depends to the amount of XB arrays*/):mems(_mems){
        numLines = numLines;
        cursor = 0;
        latency = 100;
        lineSize = lineSize/4;
        bufferEntry = new BufferEntry[numLines];
        for (int i = 0; i < numLines; ++i) {
            bufferEntry[i].clear();
        }
    }
    Buffer(){
        numLines = 1024;
        cursor = 0;
        latency = 100;
        lineSize = 1024/4;
        bufferEntry = new BufferEntry[numLines];
        for (int i = 0; i < numLines; ++i) {
            bufferEntry[i].clear();
        }
    }

    void setBufferLine(uint32_t BUFFER_LINE_BYTES){
        lineSize = BUFFER_LINE_BYTES;
    }

    void setFlags(Address address, uint32_t idx, OperationTypes type, uint64_t cycle) {
        bufferEntry[idx].address = address;
        bufferEntry[idx].isUseless = false;
        bufferEntry[idx].type = type;
        bufferEntry[idx].memCycle = cycle;
        cursor = idx;
    }

    void setUseless(uint32_t idx, uint64_t cycle){
        bufferEntry[idx].isUseless = true;
        bufferEntry[idx].memCycle = cycle;
        bufferEntry[idx].address = NULL;
    }

    uint32_t findOffset(){
        int i = 0;
        for (i = 0; i < numLines; ++i) {
            if(bufferEntry[(cursor-i)%numLines].address != 0 && bufferEntry[(cursor-i)%numLines].type == Map)
                continue;
            else break;
        }
        return (cursor-i+1)%numLines;

    }

    int findAvailable(uint64_t memCycle){
        for (int i = 0; i < numLines; ++i) {
            if(bufferEntry[i].isUseless && bufferEntry[i].memCycle <= memCycle){
                uselessIdx = i;
                return i;
            }
        }

        return -1;
    }

    uint64_t wait(){
        uint64_t minCycle = 0;
        for (int i = 0; i < numLines; ++i) {
            if(bufferEntry[i].isUseless){
                if(minCycle == 0)
                    minCycle = bufferEntry[i].memCycle;
                else{
                    if(bufferEntry[i].memCycle<minCycle){
                        minCycle = bufferEntry[i].memCycle;
                    }

                }
            }
        }
        return minCycle;
    }

    bool inline isUseless(int i){
        return bufferEntry[i].isUseless;
    }
    uint64_t inline get(int i){
        return bufferEntry[i].memCycle;
    }

    uint64_t access(MemReq& req){
       // switch (req.type) {
       //     case PUTS:
       //     case PUTX:
       //         *req.state = I;
       //         break;
       //     case GETS:
       //         *req.state = req.is(MemReq::NOEXCL)? S : E;
       //         break;
       //     case GETX:
       //         *req.state = M;
       //         break;

       //     default: panic("!?");
       // }
        if(req.type == PUTX)
            return mems->access(req);

        for (int i = 0; i < numLines; ++i) {
            if(bufferEntry[i].address==req.lineAddr && !bufferEntry[i].isUseless)
                return req.cycle;
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

    ~Buffer(){
        delete[] bufferEntry;
    }

};


#endif //ZSIM_BUFFER_H
