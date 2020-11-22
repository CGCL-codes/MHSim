//
// Created by jiahong on 19-11-23.
//

#ifndef ZSIM_TILE_H
#define ZSIM_TILE_H
#include "../memory_hierarchy.h"
#include "neuro/Memristor_CU.hpp"
#include "Cell.h"
#include "coherence.h"
#include <map>
#include "alu.h"
#include "buffer.h"
#include <queue>
#include <iostream>
using namespace std;

template<class memoryType>
class Tile{
public:
    Tile(){
        PEs = new std::map<Address,Memristor_CU<memoryType>*>();
        req_queue = new queue<AccReq>();
        memCycle = 0;
        arrayCycle = 0;
        mapCycle = 0;
        access_count = 0;
        readCycle = 0;
        computeCycle = 0;
        count = 0;
        buffer = new Buffer();
        futex_init(&filterLock);
        capacity = 1000;
        available_capacity = 1000;
        col_size = 128;
        row_size = 128;
        bits = 2;
        cells = 6;
        isEnd = false;
    }


    Tile(Buffer *_buffer, uint16_t _bits, uint16_t _row_size, uint16_t _col_size,uint16_t _capacity, uint16_t _cells):buffer(_buffer),bits(_bits),row_size(_row_size),col_size(_col_size),capacity(_capacity),cells(_cells){
        PEs = new std::map<Address,Memristor_CU<memoryType>*>();
        req_queue = new queue<AccReq>();
        memCycle = 0;
        arrayCycle = 0;
        mapCycle = 0;
        access_count = 0;
        readCycle = 0;
        computeCycle = 0;
        count = 0;
        futex_init(&filterLock);
        futex_init(&mtx);
        futex_init(&process_mtx);
        //capacity = 1000;
        available_capacity = _capacity;
        isEnd = false;
        cnt = 0;
    }

    void  setXBs(uint16_t bits,uint16_t row_size,uint16_t col_size){
        this->bits = bits;
        this->row_size = row_size;
        this->col_size = col_size;
    }

    Memristor_CU<memoryType>* find(Address addr)
    {
        if(PEs->count(addr))
            return (*PEs)[addr];
        else
            return NULL;
    }
    void add(Address addr,Memristor_CU<memoryType>* mcu)
    {
        PEs->insert(std::pair<Address,Memristor_CU<memoryType>*>(addr,mcu));
    }
    void map(Address weight,int offset){

    }
//    void compute(Address weight, Address operands, int offset){
//        Memristor_CU<memoryType>* mcu = find(weight);
//        mcu->Latency();
//    }
    void addQueue(AccReq accReq){
        futex_lock(&mtx);
        req_queue->push(accReq);

        //std::cout<<"1111111     reqSize = "<<req_queue->size()<<std::endl;

        futex_unlock(&mtx);

        //req_queue->size()
    }
    void test(){}
    void notifyProcess(){
        //thread(&Tile::test,this);
        //work_thread = std::thread(&Tile::test,this);
        //pthread_create(&Tile::wt,NULL,prefetch);
        //pthread_join(wt,NULL);
        //if(work_thread.joinable())
        //    work_thread.join();
    }



    void printState(){
        std::cout<<"MemCycle = "<<memCycle<<"\t arrayCycle = "<<arrayCycle<<std::endl;
        std::cout<<"Waiting count = "<<count<<"\t access count = "<<access_count<<std::endl;
        std::cout<<"ReadCycle = "<<readCycle<<"\t MapCycle = "<<mapCycle<<"\t ComputeCycle = "<<computeCycle<<std::endl;
    }

    void prefetch(){

        futex_lock(&mtx);
        bool a = req_queue->empty();
        futex_unlock(&mtx);
        if(a){
            usleep(10);
            return;
        }
        futex_lock(&process_mtx);

        AccReq  accReq = req_queue->front();
        futex_lock(&mtx);
        req_queue->pop();
        futex_unlock(&mtx);


        if(accReq.mode == WMI){
            if(accReq.type == Map){

                //First lock (Writer)
                //futex_lock(&filterLock);
                //Fetch Weight
                for(int i = 0; i < accReq.K; i++){
                    int idx = buffer->findAvailable(memCycle);

                    if(idx == -1){
                        memCycle = buffer->wait();
                        count ++;
                        idx = buffer->findAvailable(memCycle);

                        buffer->setFlags(accReq.lineAddr,idx,Map,memCycle);
                    }
                    uint64_t temp = memCycle;
                    buffer->setBufferLine(accReq.M);
                    for (int j = 0; j < accReq.M; ++j) {
                        //TODO: Access Memory
                        //Access all the columns at the i row
                        access_count ++;
                        MESIState dummyState = MESIState::I;
                        //MemReq req = {accReq.lineAddr+(/*j*accReq.K+*/i), GETS, 0, &dummyState, memCycle, &filterLock, dummyState, 0, 0};
                        MemReq req = {accReq.lineAddr+j*accReq.K, GETS, 0, &dummyState, memCycle, &filterLock, dummyState, 0, 0};
                        memCycle = buffer->access(req);
                    }
                    readCycle += memCycle - temp;
                    //if(readCycle > memCycle){
                    //    printf("rdCycles = %llu   memCycle = %llu  temp = %llu   delta = %llu  \n",readCycle,memCycle,temp,memCycle-temp);
                    //}
                    //Total rows as one batch
                    //Trigger write drivers to write cells
//                    std::cout<<"idx    "<<idx<<std::endl;
                    process(MAX(memCycle,arrayCycle),accReq);
                    buffer->setUseless(idx,arrayCycle);

                }

                //futex_unlock(&filterLock);
                //Unlock
            }
            else{
                //Fetch Input

                for(int i = 0; i < accReq.M /*It means the # of columns of Input*/; i++){
                    int idx = buffer->findAvailable(memCycle);

                    if(idx == -1){
                        memCycle = buffer->wait();
                        count ++;
                        idx = buffer->findAvailable(memCycle);
                        buffer->setFlags(accReq.lineAddr,idx,Compute,memCycle);
                    }
                    uint64_t temp = memCycle;
                    //printf("5 %u ",(unsigned int)pthread_self());

                    //futex_lock(&filterLock);
                    buffer->setBufferLine(accReq.K);
                    for (int j = 0; j < accReq.K; ++j) {
                        //TODO: Access Memory
                        //Access all the rows at the j column
                        //memCycle += buffer->access(accReq.lineAddr+(j*accReq.M+i));
                        MESIState dummyState = MESIState::I;
                        access_count ++;
                        //MemReq req = {accReq.lineAddr+(/*j*accReq.M+*/i), GETS, 0, &dummyState, memCycle, &filterLock, dummyState, 0, 0};
                        MemReq req = {accReq.lineAddr+j*accReq.M, GETS, 0, &dummyState, memCycle, &filterLock, dummyState, 0, 0};
                        memCycle = buffer->access(req);

                    }
                    readCycle += memCycle - temp;
                    //Total rows as an input vector to multiply with matrix
                    //Trigger read drivers to read cells and other components to compose the result
//                    std::cout<<"idx    "<<idx<<std::endl;
                    process(MAX(memCycle,arrayCycle),accReq);
                    buffer->setUseless(idx,arrayCycle);

                    // write back
                    MESIState dummyState = MESIState::I;
                    MemReq wbReq = {accReq.resultAddr+accReq.M, PUTX, 0, &dummyState, memCycle, &filterLock, dummyState,0,0};
                    memCycle = buffer->access(wbReq);
                    //lock to perform read (Reader)
                    //process(MAX(memCycle,arrayCycle));
                    //unlock
                }

                //futex_unlock(&filterLock);
            }
            //printf("6 %u\n",(unsigned int)pthread_self());

        }
        else{
            if(accReq.type == Map){

                //First lock (Writer)

                for (int i = 0; i < accReq.K /*It means the # of columns of Weight*/; ++i) {
                    buffer->setBufferLine(accReq.M);
                    for (int j = 0; j < accReq.M; ++j) {
                        //Access all the rows at the i column
                        //memCycle += buffer->access(accReq.lineAddr+(i*accReq.M+j));
                        access_count ++;
                    }
                    //Total rows as one batch
                    //buffer->setFlags(idx,Map);
                    //process(MAX(memCycle,arrayCycle));
                }

                //unlock
            }
            else{
                for (int i = 0; i < accReq.M; ++i) {
                    buffer->setBufferLine(accReq.K);
                    for (int j = 0; j < accReq.K; ++j) {
                        //Access all the columns at the i rows
                        //memCycle += buffer->access(accReq.lineAddr+(i*accReq.K+j));
                        access_count ++;
                    }
                    //Trigger read drivers to read cells and other components to compose the result
                    //buffer->setFlags(idx,Compute);

                    //lock to perform read (Reader)
                    //process(MAX(memCycle,arrayCycle));
                    //unlock
                }

            }

        }



        futex_unlock(&process_mtx);

        //TODO: start Map or Compute





    }


    void process(uint64_t curCycle, AccReq arq){
        //buffer->cursor
        //futex_lock(&mtx);
        //AccReq arq = req_queue->front();
        //futex_unlock(&mtx);

        if(arq.type == Map){
            Memristor_CU<RealDevice>* mcu = find(arq.lineAddr);
            if(mcu==NULL)
            {
                printf("before\n");
                mcu = new Memristor_CU<RealDevice>(arq.M, arq.K, reinterpret_cast<float *>(arq.lineAddr), arq.mode, bits, row_size, col_size, cells);
                printf("initialized\n");
                add(arq.lineAddr,mcu);
                available_capacity -= mcu->getSize();
            }
            arrayCycle = curCycle + mcu->getWriteLatency(buffer->findOffset()) * 2e9;
            mapCycle += mcu->getWriteLatency(buffer->findOffset()) * 2e9;
        }else{
            //printf("compute\n");
            while(find(arq.otherAddr) == NULL)
                usleep(10);
            Memristor_CU<RealDevice>* mcu = find(arq.otherAddr);
            if(mcu->getLatency() == 0)
                mcu->ReadLatency();
            arrayCycle = curCycle + mcu->getLatency() * 2e9;
            computeCycle += mcu->getLatency() * 2e9; //2e9 is the clock cycle in NeuroSim
        }

    }


    uint16_t getAvailableCapacity(){
        return available_capacity;
    }

    uint16_t caculateSize(uint32_t M, uint32_t K, uint32_t cells, uint32_t colSize, uint32_t rowSize){
        return (int)ceil((double)K/rowSize)*(int)ceil((double)M*cells/colSize);
    }

    void end(){
        isEnd = true;
    }

    void wait_for_end(){
        while(!req_queue->empty()){
            usleep(100);
        }

        isEnd = true;

    }

    ~Tile(){
        delete[] PEs;
        delete buffer;
        std::cout<<"req remain:"<<req_queue->size()<<std::endl;
        //work_thread.join();
    }

    uint16_t getBits(){
        return bits;
    }
    uint16_t getRowSize(){return row_size;}
    uint16_t getColSize(){return col_size;}
    uint16_t getCapacity(){return capacity;}
private:
    //TODO::
    Buffer *buffer;
    std::map<Address,Memristor_CU<memoryType>*> *PEs;
    queue<AccReq> *req_queue; //control_unit just sends the req to each tiles
    //TODO::
    ALU alu;
    //uint16_t capacity,available_capacity;  //number of crossbar
    uint64_t memCycle,readCycle;
    uint64_t arrayCycle,mapCycle,computeCycle;
    uint64_t count,access_count;
    uint16_t bits,row_size,col_size,capacity,available_capacity,cells;
    //std::thread work_thread;
    //pthread_t wt;
    bool isEnd;
    lock_t mtx,process_mtx;
    lock_t filterLock;
    uint64_t cnt;
};

#endif //ZSIM_TILE_H
