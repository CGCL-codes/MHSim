/*
 * @Author: jhxu
 * @LastEditors: jhxu
 * @LastEditTime: 2022-01-05 19:13:33
 * @FilePath: /src/MHSim/chip.cpp
 */
#include "chip.h"

//Chip<RealDevice> *chip = new Chip<RealDevice>();

void* notifyProcess(void  *args){
    Tile<RealDevice>* tile = (Tile<RealDevice>*)args;
    while(!Chip<RealDevice>::getIsEnd()){
        tile->prefetch();
        usleep(3000);
    }
    return NULL;
}

template <class memoryType>
bool Chip<memoryType>::isEnd = false;

template <>
bool Chip<RealDevice>::premap(Matrix m, std::vector<XB> *xbs, uint16_t depth){
    bool is_fully_map = false;
    std::vector<XB> temp;
    if(depth == 0){
        temp = std::vector<XB>(*xbs);
    }


    uint64_t ac = 0;
    for(int i = 0; i < xbs->size(); ++i){
        ac += (*xbs)[i].available_capacity();
    }
    //printf("M.size = %d * %d \n", m.K, m.M);
    if(m.size() > ac)
        return false;
    else{
        for(int i = 0; i < xbs->size(); ++i){
            if((*xbs)[i].compare(m)){
                (*xbs)[i].map(m);
                sort(xbs->begin(), xbs->end());
                is_fully_map = true;
                break;
            }
        }
        if(!is_fully_map){
            Matrix mm = m - (*xbs)[xbs->size()-1];
            //printf("M.size = %d * %d  %d * %d \n", m.K, m.M, mm.K, mm.M);
            sort(xbs->begin(), xbs->end());
            bool b = false;
            b = premap(m, xbs, depth+1);
            if(!mm.isNull())
                b &= premap(mm, xbs, depth+1);
            is_fully_map = b;
        }
    }



    if(depth == 0){
        if(!is_fully_map){
            xbs->assign(temp.begin(), temp.end());
        }
    }
    return is_fully_map;
}
/**
 * @description: Try to map the matrix into multiple tiles
 * @param {*}
 * @return {*} Whether the matrix is entirely mapped
 */
template <>
bool Chip<RealDevice>::mapMultiTiles(Matrix m, int &idx){
    std::vector<XB> *xbs = new std::vector<XB>, *xbs_temp = new std::vector<XB>, *xb_ti, **xb_tiles;
    xb_ti = (*tiles)[0]->getXBs();
    if((*tiles)[0]->get_cap() != 0)
    xbs_temp->insert(xbs_temp->end(), xb_ti->begin(), xb_ti->end());
    bool isMap = false;
    uint64_t cap = 0;
    cap += (*tiles)[0]->get_cap();
    for(int i = 1; i < tiles->size(); ++i){
        if((*tiles)[i]->get_cap() == 0){
            continue;
        }
        cap += (*tiles)[i]->get_cap();
        xb_ti = (*tiles)[i]->getXBs();
        xbs_temp->insert(xbs_temp->end(), xb_ti->begin(), xb_ti->end());
        if(cap < m.size()){
            continue;
        }
        isMap = premap(m, xbs_temp, 0);
        if(isMap){
            idx = i;
            // printf("  Used tiles = %d\n", i+1);
            xb_tiles = new std::vector<XB>*[i+1];
            for(int j = 0; j <= i; ++j){
                if((*tiles)[j]->get_cap() != 0){
                    xb_tiles[j] = (*tiles)[j]->getXBs();
                    std::vector<XB> x;
                    xb_tiles[j]->swap(x);
                }

            }

            for(int j = 0; j < xbs_temp->size(); ++j){
                XB xb = (*xbs_temp)[j];
                xb_tiles[xb.num_tiles]->push_back(xb);
            }
            for(int j = 0; j <= i; ++j){
                (*tiles)[j]->syn();
            }
            delete[] xb_tiles;
            break;
        }
    }

    delete xbs;
    delete xbs_temp;
    return isMap;
}

template <>
bool Chip<RealDevice>::memristor_mm(Address weight, Address operands, Address output, uint32_t M, uint32_t N, uint32_t K){
    mapinfo* minfo = look_up(weight);
    if(minfo == NULL)
    {
        uint32_t size = M * K * cells;
        //printf("M = %d N = %d K = %d\n", M, N, K);
        AccReq accReq = {weight,Map,K,M,N,IMW,operands,output};
        mcc->startAccess(accReq);
        if(mcc->shouldAllocate(accReq)){ // Decide whether we should allocate resources to map the request.
            uint32_t replications = repl[repl_index++];
            if(N == 100)
                replications = 1;
            printf("Matrix %d*%d  replications = %d   %d  N = %d \n", M, K, replications, replications * M * K * 8/128/128, N);

            for(int iii = 0; iii < replications; ++iii){
                int idx = -1;
                Matrix m(K, M * cells);
                bool bb = false;
                for(int i = 0; i < tiles->size(); ++i){
                    if((*tiles)[i]->premap(m,0)){
                        idx = i;
                        bb = true;
                        break;
                    }
                }
                if(bb){

                    if(iii == 0){
                        // printf("map into a tile \n");
                        addMapInfo(weight,idx,size,0,replications);
                        (*tiles)[idx]->addQueue(accReq);
                        pthread_create(&pt[idx],NULL,notifyProcess,(*tiles)[idx]);
                        uint64_t mapping_time = 0;
                        while(!(*tiles)[idx]->isMapped(weight));

                        mapping_time = (*tiles)[idx]->time(weight);
                        uint64_t core0Cycle = zinfo->cores[0]->getCycles();
                        (*mmap)[weight]->completeCycle = core0Cycle + mapping_time;
                    }

                    if(iii == replications - 1){
                        (*tiles)[(*mmap)[weight]->idx]->addQueue({operands,Compute,K,M,N,IMW,weight,output,replications});
                        return true;
                    }
                }
                else{
                    if(mapMultiTiles(m, idx)){

                        if(iii == 0){
                            addMapInfo(weight,idx,size,0,replications);                            
                            (*tiles)[idx]->addQueue(accReq);
                            pthread_create(&pt[idx],NULL,notifyProcess,(*tiles)[idx]);
                            uint64_t mapping_time = 0;
                            while(!(*tiles)[idx]->isMapped(weight));

                            mapping_time = (*tiles)[idx]->time(weight);
                            uint64_t maxCycle = zinfo->cores[0]->getCycles();
                            (*mmap)[weight]->completeCycle = maxCycle + mapping_time;
                        }

                        if(iii == replications - 1){
                            (*tiles)[(*mmap)[weight]->idx]->addQueue({operands,Compute,K,M,N,IMW,weight,output,replications});
                            return true;
                        }
                    }
                }
            }
            printf("Matrix(%dx%d)   Address = %llX is too large to be mapped into tile. Use CPU to process it. \n", K, M, weight);
                    //Deliver to processors.
                    //...
                    //
            return false; 
            


        }
        else{
            //use CPU
            return false;
        }

    }
    else{
        if(minfo->completeCycle < zinfo->cores[0]->getCycles() + zinfo->phaseLength){
            // printf("replications = %d \n", (*mmap)[weight]->replications);
            // printf("dims = %d \n", M);
            AccReq accReq = {operands,Compute,K,M,N,IMW,weight,output,(*mmap)[weight]->replications};
            (*tiles)[minfo->idx]->addQueue(accReq);

        }else{
            return false;
        }
    }
    return true;
}
template <>
void Chip<RealDevice>::init(){
    for(int i = 0; i < 100; i++) repl[i] = 1;
    // repl[0] = 128;
}
template <>
void Chip<RealDevice>::initStats(AggregateStat* parentStat) {
    AggregateStat* mbaStat = new AggregateStat();
    mbaStat->init("MBA", "Memristor-based accelerator chip");

    auto x = [this]() { return getTotalCycle(); };
    LambdaStat<decltype(x)>* cyclesStat = new LambdaStat<decltype(x)>(x);
    cyclesStat->init("mvm cycles", "The total latency of processing all MVMs");

    auto y = [this]() { return getMemAccessCycle(); };
    LambdaStat<decltype(y)>* mCyclesStat = new LambdaStat<decltype(y)>(y);
    mCyclesStat->init("mem access", "Memory access latency during MVMs");


    mbaStat->append(cyclesStat);
    mbaStat->append(mCyclesStat);
    parentStat->append(mbaStat);
}