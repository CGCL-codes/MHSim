//
// Created by jiahong on 19-11-23.
//
#include "chip.h"

//Chip<RealDevice> *chip = new Chip<RealDevice>();

void* notifyProcess(void  *args){
    Tile<RealDevice>* tile = (Tile<RealDevice>*)args;
    while(!Chip<RealDevice>::getIsEnd())
        tile->prefetch();
    return NULL;
}


template <class memoryType>
bool Chip<memoryType>::isEnd = false;


template <>
bool Chip<RealDevice>::memristor_mm(Address weight, Address operands, Address output, uint32_t M, uint32_t N, uint32_t K){
    mapinfo* minfo = look_up(weight);
    if(minfo == NULL)
    {
        uint16_t size = ceil((double) K / row_size) * ceil((double) M / (int)(col_size/(cells)));
        //printf("Size = %d \n",size);
        AccReq accReq = {weight,Map,K,M,size,WMI};//TODO::
        //std::cout<<std::hex<<weight<<std::endl;
        mcc->startAccess(accReq);

        if(mcc->shouldAllocate(accReq)){ // Decide whether we should allocate resources to map the request.
            int idx = findAvailableTile(size);
            if(idx!=-1){
                //std::cout<<idx<<endl;
                addMapInfo(weight,idx,size);
                printf("Map \n");
                (*tiles)[idx]->addQueue(accReq);
                (*tiles)[idx]->addQueue({operands,Compute,K,N,0,WMI,weight,output});
                //Tile<RealDevice>* t = (*tiles)[idx];
                //assert(t);
                //std::thread t(notifyProcess,iii);
                //printf("before create \n");
                pthread_create(&pt[idx],NULL,notifyProcess,(*tiles)[idx]);
                //printf("created \n");
                //t.join();
                //tile_threads[idx] = thread(notifyProcess,t);
                //(*tiles)[idx]->notifyProcess();


                //for(N){ add queue (compute)};


            }
            else{
                //Perform evicting
                //...
                //
                return false; //comment it when completed function shouldAllocate()
                //tiles

            }
            //  map(weight,offset1);



        }
        else{
            //TODO:use CPU
            return false;
        }


    }
    else{
        AccReq accReq = {operands,Compute,K,N,0,WMI,weight,output};
        (*tiles)[minfo->idx]->addQueue(accReq);
         // (tiles)[minfo->idx]->compute(weight,operands,N);

    }
    return true;
}
