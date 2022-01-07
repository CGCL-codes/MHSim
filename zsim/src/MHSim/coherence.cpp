/*
 * @Author: jhxu
 * @LastEditors: jhxu
 * @LastEditTime: 2022-01-05 17:23:07
 * @FilePath: /src/MHSim/coherence.cpp
 */

#include "coherence.h"

int MemristorMESICC::startAccess(AccReq& req) {
//    futex_lock(&mlock);
    return -1;//return the index of tile that contains the resource.
}


bool MemristorMESICC::shouldAllocate(AccReq& req){
    return true; // true when the req's resources are specified
}
void MemristorMESICC::processEviction() {

}
void MemristorMESICC::processAccess(AccReq& req) {
    if(req.type==Compute){
        //rwm.rdLock();
    }
    else if(req.type == Map){
        //rwm.wrLock();

        //rwm.wrUnlock();
    }

}
void MemristorMESICC::endAccess() {

    //rwm.rdUnlock();
//    futex_unlock(&mlock)
}
