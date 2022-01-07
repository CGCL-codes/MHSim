/*
 * @Author: jhxu
 * @LastEditors: jhxu
 * @LastEditTime: 2022-01-05 21:38:52
 * @FilePath: /src/MHSim/buffer.h
 */


#ifndef ZSIM_BUFFER_H
#define ZSIM_BUFFER_H

#include "../galloc.h"
#include "coherence.h"
#include "GEMM/Memristor_CU.h"
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
        this->numLines = numLines;
        cursor = 0;
        latency = 100;
        this->lineSize = lineSize;
        bufferEntry = new BufferEntry[numLines];
        for (int i = 0; i < numLines; ++i) {
            bufferEntry[i].clear();
        }
    }
    Buffer(){
        numLines = 1024;
        cursor = 0;
        latency = 100;
        lineSize = 1024;
        bufferEntry = new BufferEntry[numLines];
        for (int i = 0; i < numLines; ++i) {
            bufferEntry[i].clear();
        }
    }



    void setFlags(Address address, uint32_t idx, OperationTypes type, uint64_t cycle) {
        bufferEntry[idx].address = address;
        bufferEntry[idx].isUseless = true;

        //bufferEntry[idx].isUseless = false; // We use FIFO strategy.
        bufferEntry[idx].type = type;
        bufferEntry[idx].memCycle = cycle;
        cursor = idx;
    }

    void setUseless(uint32_t idx, uint64_t cycle){
        bufferEntry[idx].isUseless = true;
        bufferEntry[idx].memCycle = cycle;
        bufferEntry[idx].address = NULL;
    }

    uint32_t getLineSize(){
        return lineSize;
    }
    int findAvailable(uint64_t memCycle){
        for (int i = 0; i < numLines; ++i) {
            if(bufferEntry[i].isUseless && bufferEntry[i].memCycle <= memCycle){
                //uselessIdx = i;
                return i;
            }
        }
        return -1;
    }

    uint64_t wait(int &idx);

    bool inline isUseless(int i){
        return bufferEntry[i].isUseless;
    }
    uint64_t inline get(int i){
        return bufferEntry[i].memCycle;
    }

    uint64_t access(MemReq& req);

    uint32_t getUselessIdx(){
        return uselessIdx;
    }

    ~Buffer(){
        delete[] bufferEntry;
    }

};


#endif //ZSIM_BUFFER_H
