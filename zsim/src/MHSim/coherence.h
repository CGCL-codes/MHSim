/*
 * @Author: jhxu
 * @LastEditors: jhxu
 * @LastEditTime: 2022-01-05 21:38:56
 * @FilePath: /src/MHSim/coherence.h
 */



#ifndef ZSIM_COHERENCE_H
#define ZSIM_COHERENCE_H

#include "../memory_hierarchy.h"
#include "../locks.h"
#include "../mutex.h"
#include "GEMM/Memristor_CU.h"


typedef enum {
    Map = 1,
    Compute = 2
} OperationTypes;

struct AccReq{
    Address lineAddr;
    OperationTypes type;
    uint32_t K;
    uint32_t M;
    uint32_t N;
    OperationType mode;
    Address otherAddr;
    Address resultAddr;
    uint32_t replications;
};

class MemristorMESICC {
private:
    MESIState* array;
    uint32_t numLines;
    lock_t mlock;
    rwmutex rwm;
public:
    MemristorMESICC(uint32_t _numLines) : numLines(_numLines){
        futex_init(&mlock);
    }

    MemristorMESICC(){
        numLines = 128;
        futex_init(&mlock);
    }

    int startAccess(AccReq& req);//initial locking, address races; returns true if access should be skipped; may change req!
    bool shouldAllocate(AccReq& req); //called when we don't find req's lineAddr in the array             look_up
    void processEviction();//called if shouldAllocate returns true
    void processAccess(AccReq& req);
    void endAccess();
    void processWriteBack();

};
#endif //ZSIM_COHERENCE_H
