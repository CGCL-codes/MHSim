//
// Created by jiahong on 20-2-18.
//

#ifndef ZSIM_COHERENCE_H
#define ZSIM_COHERENCE_H

#include "../memory_hierarchy.h"
#include "../locks.h"
#include "../mutex.h"
#include "neuro/Memristor_CU.h"


typedef enum {
    Map,
    Compute
} OperationTypes;

struct AccReq{
    Address lineAddr;
    OperationTypes type;
    uint32_t K;
    uint32_t M;
    uint32_t size;
    OperationType mode;
    Address otherAddr;
    Address resultAddr;
};
//
//
//
class MemristorMESICC {
private:
    MESIState* array;
    uint32_t numLines;
    lock_t mlock;
    rwmutex rwm;
public:
    MemristorMESICC(uint32_t _numLines) : numLines(_numLines){
//        array = gm_calloc<MESIState>(numLines);
//        for (unsigned int i = 0; i < numLines; ++i) {
//            array[i] = I;
//        }
        futex_init(&mlock);
    }

    MemristorMESICC(){
        numLines = 128;
//        array = gm_calloc<MESIState>(numLines);
//        for (unsigned int i = 0; i < numLines; ++i) {
//            array[i] = I;
//        }
        futex_init(&mlock);
    }

    int startAccess(AccReq& req);//initial locking, address races; returns true if access should be skipped; may change req!
    bool shouldAllocate(AccReq& req); //called when we don't find req's lineAddr in the array             look_up
    void processEviction();//called if shouldAllocate returns true
    void processAccess(AccReq& req);
    void endAccess();
    void processWriteBack();

};
//
//
//struct coherence_table{
//    uint8_t offset;
//
//};
//
#endif //ZSIM_COHERENCE_H
