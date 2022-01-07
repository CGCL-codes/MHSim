/*
 * @Author: jhxu
 * @LastEditors: jhxu
 * @LastEditTime: 2022-01-05 19:29:02
 * @FilePath: /src/MHSim/chip.h
 */

#ifndef ZSIM_CHIP_H
#define ZSIM_CHIP_H

#include <map>
#include "tile.h"
#include "Cell.h"
#include "coherence.h"
#include <thread>
#include <pthread.h>
#include "../core.h"
#include "../zsim.h"
#include "../stats.h"
using namespace std;

struct mapinfo{
    Address address;
    uint16_t idx;
    uint16_t size; //
    uint64_t completeCycle;
    uint32_t replications;
};

template <class memoryType>
class Chip{
public:
    Chip(){
        tiles = new std::vector<Tile<memoryType>*>(10);
        this->energy = 0.;
        this->latency = 0.;
        mcc = new MemristorMESICC();
        mmap = new std::map<Address ,mapinfo*>();
        row_size = 128;
        col_size = 128;
        tidx = 0;
        for (int i = 0; i < 10; ++i) {
            Tile<memoryType> *tile = new Tile<memoryType>();
            add(i,tile);
        }
        cells = (*tiles)[0]->getCells();
        pt = new pthread_t[10];
        repl = new uint32_t[50];
        for(int i = 0; i < 100; i++) repl[i] = 1;
        repl_index = 0;
    }

    void initStats(AggregateStat* parentStat);

    Chip(std::vector<Tile<memoryType>*> *_tiles):tiles(_tiles){
        energy = 0.;
        latency = 0.;
        cells = (*_tiles)[0]->getCells();
        row_size = (*_tiles)[0]->getRowSize();
        col_size = (*_tiles)[0]->getColSize();
        mcc = new MemristorMESICC();
        mmap = new std::map<Address,mapinfo*>();
        tidx = 0;
        pt = new pthread_t[(*_tiles).size()];
        repl = new uint32_t[50];

        repl_index = 0;
        init();
    }


    void init();
    Tile<memoryType>* find(int addr){
        if((*tiles)[addr]!=NULL)
            return (*tiles)[addr];
        else
            return NULL;
    }

    // void printAll(){
    //     uint64_t memcycle = 0, xbcycle = 0;
    //     for (int i = 0; i < tiles->size(); ++i) {
    //         memcycle += (*tiles)[i]->getMemCycle();
    //         xbcycle += (*tiles)[i]->getXBCycle();
    //     }
    //     printf("MemCycle = %llu  , XBcycle = %llu \n", memcycle, xbcycle);
    // }

    uint64_t getTotalCycle(){
        uint64_t xbcycle = 0;
        for (int i = 0; i < tiles->size(); ++i) {
            xbcycle += (*tiles)[i]->getXBCycle();
        }
        return xbcycle;
    }

    uint64_t getMemAccessCycle(){
        uint64_t memcycle = 0;
        for (int i = 0; i < tiles->size(); ++i) {
            memcycle += (*tiles)[i]->getMemCycle();
        }
        return memcycle;
    }
    
    bool mapMultiTiles(Matrix m, int &idx);

    bool memristor_mm(Address weight, Address operands, Address output, uint32_t M, uint32_t N, uint32_t K);


    void add(int addr,Tile<memoryType> *tile)
    {
       (*tiles)[addr] = tile;
    }
    void addLatency(double latency){
        this->latency += latency;
    }
    void addEnergy(double energy){
        this->energy += energy;
    }
    double getLatency(){
        return this->latency;
    }
    double getEnergy(){
        return this->energy;
    }
    static inline bool getIsEnd(){
        return isEnd;
    }


    mapinfo* look_up(Address addr){
        if(mmap->count(addr))
            return (*mmap)[addr];
        else return NULL;
    }

    void addMapInfo(Address weight, uint16_t idx, uint16_t size, uint64_t completeCycle, uint32_t replications = 1){
        mapinfo *minfo = new mapinfo{weight,idx,size, completeCycle, replications};
        mmap->insert(std::pair<Address,mapinfo*>(weight,minfo));
    }

    void map(Address addr,int offset){

    }

    void compute(Address addr,int offset){

    }
    void wait_for_end(){

        for (int i = 0; i < (*tiles).size(); ++i) {
            // printf("wait for tile %d \n",i);
            (*tiles)[i]->wait_for_end();
        }
        isEnd = true;
    }

    bool getMask(){return Masked;}
    void setMask(bool m){Masked = m;}

    bool premap(Matrix m, std::vector<XB> *xbs, uint16_t depth = 0);

    ~Chip(){
        void* rtval[tiles->size()];
        for(int i = 0 ; i < tiles->size(); ++i) {pthread_join(pt[i],&rtval[i]);}//delete tiles[i];}
        delete tiles;
        delete[] mcc;
        delete[] mmap;
        delete[] pt;
    }
private:
    std::vector<Tile<memoryType>*> *tiles;
    bool Masked;
    double latency,energy;
    uint16_t cells, col_size, row_size,tidx;
    MemristorMESICC *mcc;
    std::map<Address ,mapinfo*> *mmap;
    static bool isEnd;
    pthread_t *pt;
    uint32_t *repl;
    int repl_index;
};

#endif //ZSIM_CHIP_H
