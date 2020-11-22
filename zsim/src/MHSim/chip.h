//
// Created by jiahong on 19-11-23.
//

#ifndef ZSIM_CHIP_H
#define ZSIM_CHIP_H

#include <map>
#include "tile.h"
#include "Cell.h"
#include "coherence.h"
#include <thread>
#include <pthread.h>
using namespace std;

struct mapinfo{
    Address address;
    uint16_t idx;
    uint16_t size; //
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
        cells = 6;
        row_size = 128;
        col_size = 128;
        tidx = 0;
        for (int i = 0; i < 10; ++i) {
            Tile<memoryType> *tile = new Tile<memoryType>();
            add(i,tile);
        }
        pt = new pthread_t[10];
    }

    Chip(std::vector<Tile<memoryType>*> *_tiles):tiles(_tiles){
        energy = 0.;
        latency = 0.;
        cells = 6;
        row_size = (*_tiles)[0]->getRowSize();
        col_size = (*_tiles)[0]->getColSize();
        mcc = new MemristorMESICC();
        mmap = new std::map<Address,mapinfo*>();
        tidx = 0;
        pt = new pthread_t[(*_tiles).size()];
    }


    Tile<memoryType>* find(int addr){
        if((*tiles)[addr]!=NULL)
            return (*tiles)[addr];
        else
            return NULL;
    }

    void printAll(){
        for (int i = 0; i < tiles->size(); ++i) {
            std::cout<<i<<" tile  "<<(((*tiles)[i])==NULL)<<std::endl;
            (*tiles)[i]->printState();
        }
    }

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

    void addMapInfo(Address weight, uint16_t idx, uint16_t size){
        mapinfo *minfo = new mapinfo{weight,idx,size};
        mmap->insert(std::pair<Address,mapinfo*>(weight,minfo));
    }

    void map(Address addr,int offset){

    }

    void compute(Address addr,int offset){

    }
    void wait_for_end(){
        
        for (int i = 0; i < (*tiles).size(); ++i) {
            printf("wait for tile %d \n",i);
            (*tiles)[i]->wait_for_end();
        }
        isEnd = true;
    }


    int findAvailableTile(int size){
        //typename std::map<int,Tile<memoryType>*>::iterator iter;
        //iter = tiles.begin();

        for (int i = 0; i < (*tiles).size(); i++){
            if((*tiles)[(i+tidx)%(*tiles).size()]->getAvailableCapacity()>size){
                tidx = (i+tidx+1)%(*tiles).size();
                return (tidx-1)%(*tiles).size();
            }
        }
        return -1;
    }


    ~Chip(){
        void* rtval[tiles->size()];
        for(int i = 0 ; i < tiles->size(); ++i) {pthread_join(pt[i],&rtval[i]);}//delete tiles[i];}
        delete tiles;
        delete[] mcc;
        delete[] mmap;
        delete[] pt;
    }
private:
    //TODO::
    //Control control;
    //TODO:
    //Coherence_table ct;
    std::vector<Tile<memoryType>*> *tiles;
    double latency,energy;
    uint16_t cells, col_size, row_size,tidx;
    MemristorMESICC *mcc;
    std::map<Address ,mapinfo*> *mmap;
    static bool isEnd;
    pthread_t *pt;
};

#endif //ZSIM_CHIP_H
