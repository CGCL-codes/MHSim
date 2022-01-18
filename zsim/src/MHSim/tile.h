/*
 * @Author: jhxu
 * @LastEditTime: 2022-01-05 21:39:00
 * @LastEditors: jhxu
 * @FilePath: /src/MHSim/tile.h
 */

#ifndef ZSIM_TILE_H
#define ZSIM_TILE_H
#include "../memory_hierarchy.h"
#include "GEMM/Memristor_CU.hpp"
#include "Cell.h"
#include "coherence.h"
#include <map>
#include "alu.h"
#include "buffer.h"
#include "control_unit.h"
#include <queue>
#include <iostream>
#include <algorithm>

using namespace std;

struct Item{
    uint16_t index;
    uint16_t cols;
    uint16_t rows;
    uint8_t start_col;
    uint8_t start_row;
    uint8_t size;
};

class Matrix;
class XB;

class XB{
public:
    XB(uint64_t rows, uint64_t cols, uint32_t _num_tiles = 0){remain_rows = rows; remain_cols = cols; num_tiles = _num_tiles;}
    XB(){remain_rows = 0; remain_cols = 0; num_tiles = 0;}
    XB(const XB &xb){remain_rows = xb.remain_rows; remain_cols = xb.remain_cols; num_tiles = xb.num_tiles;}
    uint64_t available_capacity(){
        return remain_rows*remain_cols;
    }
    bool operator<(XB &xb){
        if(available_capacity()<xb.available_capacity()){
            return true;
        }else{
            if(available_capacity() == xb.available_capacity()){
                return num_tiles > xb.num_tiles;
            }else{
                return false;
            }
        }
    }
    bool operator>(XB &xb){
        if(available_capacity()>xb.available_capacity()){
            return true;
        }else{
            if(available_capacity() == xb.available_capacity()){
                return num_tiles < xb.num_tiles;
            }else{
                return false;
            }
        }
    }
    bool compare(Matrix m);
    void map(Matrix m);
    void set(XB xb){
        this->remain_cols = xb.remain_cols;
        this->remain_rows = xb.remain_rows;
    }
    void print_size(){
        printf("%d * %d \n", remain_rows, remain_cols);
    }
    uint64_t remain_rows;
    uint64_t remain_cols;
    uint32_t num_tiles;
};

class Matrix{
public:
    Matrix(uint64_t K, uint64_t M){
        this->K = K;
        this->M = M;
    }
    Matrix(const Matrix &m){K = m.K; M = m.M;}
    bool operator>=(XB xb){
        return M >= xb.remain_cols && K >= xb.remain_rows;
    }
    bool operator<=(XB xb){
        return M <= xb.remain_cols && K <= xb.remain_rows;
    }
    Matrix operator-(XB &xb){
        if(xb.compare(*this)){
            printf("error: xb is larger than this matrix, Matrix = %d*%d   XB = %d*%d \n", K, M, xb.remain_rows, xb.remain_cols);
        }
        Matrix m(0,0);
        if(xb.remain_rows >= K && xb.remain_cols < M){
            M -= xb.remain_cols;
            xb.remain_rows -= K;
        }else if(xb.remain_rows < K && xb.remain_cols >= M){
            K -= xb.remain_rows;
            xb.remain_cols -= M;
        }else{ // only K >= remain_rows and M >= remain_cols
            m.M = M - xb.remain_cols;
            m.K = xb.remain_rows;
            K = K - xb.remain_rows;
            
            xb.remain_rows = 0;
            xb.remain_cols = 0;
        }
        return m;
    }
    uint64_t size(){return M*K;}
    bool isNull(){return (M == 0) && (K == 0);}
    uint64_t M;
    uint64_t K;
};


template<class memoryType>
class Tile{
public:
    Tile();


    Tile(Buffer *_buffer, uint16_t _bits, uint16_t _row_size, uint16_t _col_size,uint16_t _capacity, uint16_t _cells, uint32_t index);

    void  setXBs(uint16_t bits, uint16_t row_size, uint16_t col_size, uint16_t writeMux = 16){
        this->bits = bits;
        this->row_size = row_size;
        this->col_size = col_size;
        this->writeMux = writeMux;
        writeDrivers = row_size/writeMux;
    }

    Memristor_CU<memoryType>* find(Address addr)
    {
        if(PEs->count(addr))
            return (*PEs)[addr];
        else
            return NULL;
    }


    float* getmLat(Address addr){
        if(mLat->count(addr))
            return (*mLat)[addr];
        else return NULL;
    }

    void addmLat(Address addr, float* ml){
        mLat->insert(std::pair<Address, float*>(addr, ml));
    }

    float getcLat(Address addr){
        if(cLat->count(addr))
            return (*cLat)[addr];
        else return 0;
    }

    float addcLat(Address addr, float lat)
    {
        cLat->insert(std::pair<Address, float>(addr,lat));
    }

    bool isMapped(Address addr){
        bool b = false;
        futex_lock(&query_mtx);
        b = mtTable->count(addr);
        futex_unlock(&query_mtx);
        return b;
    }

    uint64_t time(Address addr)
    {
        uint64_t rt = 0;

        futex_lock(&query_mtx);
        if(mtTable->count(addr))
            rt = (*mtTable)[addr];

        futex_unlock(&query_mtx);
        return rt;
    }

    void add(Address addr,Memristor_CU<memoryType>* mcu)
    {
        PEs->insert(std::pair<Address,Memristor_CU<memoryType>*>(addr,mcu));
    }


    /**
     * @description: Pre-allocate the XBs
     * @param {Matrix} m
     * The weight matrix
     * @return {*}  Whether the matrix can be fully mapped into XBs
     */    
    bool premap(Matrix m, uint16_t depth = 0);
    
    void syn(){
        std::vector<XB> x;
        XB_vector_premap->swap(x);
        XB_vector_premap->insert(XB_vector_premap->end(), XB_vector->begin(), XB_vector->end());
    }

    uint64_t get_cap(){
        uint64_t cap = 0;
        for(int i = 0; i < XB_vector->size(); ++i)
            cap += (*XB_vector)[i].available_capacity();
        return cap;
    }

    void addQueue(AccReq accReq){
        futex_lock(&mtx);
        req_queue->push(accReq);
        futex_unlock(&mtx);
    }

    void printState(){
        std::cout<<"MemCycle = "<<memCycle<<"\t arrayCycle = "<<arrayCycle<<std::endl;
        std::cout<<"Waiting count = "<<count<<"\t access count = "<<access_count<<std::endl;
        std::cout<<"ReadCycle = "<<readCycle<<"\t MapCycle = "<<mapCycle<<"\t ComputeCycle = "<<computeCycle<<std::endl;
    }

    void prefetch();

    uint16_t getCells(){return cells;}


    void process(uint64_t curCycle, AccReq arq, uint16_t offset = 0);


    int64_t getAvailableCapacity(){
        int64_t available = 0;
        for(int i = 0; i < XB_capacity->size(); ++i)
            available += (*XB_capacity)[i];
        if(available - occupation < 0)
            return 0;
        return available - occupation;
    }

    uint16_t caculateSize(uint32_t M, uint32_t K, uint32_t cells, uint32_t colSize, uint32_t rowSize){
        return (int)ceil((double)K/rowSize)*(int)ceil((double)M*cells/colSize);
    }

    void end(){
        isEnd = true;
    }

    void setOccupation(uint64_t ocp){
        occupation += ocp;
    }

    void wait_for_end(){
        while(!req_queue->empty()){
            printf("req_queue size = %d \n", req_queue->size());
            usleep(100);
        }

        isEnd = true;

    }

    std::vector<XB> *getXBs(){
        return XB_vector;
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
    uint64_t getMemCycle(){return readCycle;}
    uint64_t getXBCycle(){return computeCycle;}
    uint16_t getRowSize(){return row_size;}
    uint16_t getColSize(){return col_size;}
    uint16_t getCapacity(){return capacity;}
private:
    //TODO::
    Buffer *buffer;
    std::map<Address,Memristor_CU<memoryType>*> *PEs;
    std::map<Address, float> *cLat;
    std::map<Address, float*> *mLat;
    std::map<Address, uint64_t> *mtTable;
    std::map<Address,Item> *mapping_table;
    std::vector<XB> *XB_vector, *XB_vector_premap;
    std::vector<uint32_t> *XB_capacity;
    queue<AccReq> *req_queue; //control_unit just sends the req to each tiles
    //TODO::
    Control_Unit *control_unit;
    ALU *alu;
    //uint16_t capacity,available_capacity;  //number of crossbar
    uint64_t memCycle,readCycle;
    uint64_t tileCycles;
    uint64_t arrayCycle,mapCycle,computeCycle, mpCycle;
    uint64_t count,access_count;
    uint64_t bits,row_size,col_size,cells;
    int64_t  capacity,available_capacity,occupation;
    //std::thread work_thread;
    //pthread_t wt;
    bool isEnd;
    lock_t mtx, process_mtx, query_mtx;
    lock_t filterLock;
    uint64_t cnt;
    uint16_t writeMux;
    uint16_t writeDrivers;
};

#endif //ZSIM_TILE_H
