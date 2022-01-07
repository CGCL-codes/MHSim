/*
 * @Author: jhxu
 * @LastEditors: jhxu
 * @LastEditTime: 2022-01-05 17:24:50
 * @FilePath: /src/MHSim/tile.cpp
 */
#include "tile.h"

bool XB::compare(Matrix m){return (m.K<=remain_rows)&&(m.M<=remain_cols);}
void XB::map(Matrix m){
    if(compare(m)){
       uint64_t rs = remain_rows - m.K;
        uint64_t cs = remain_cols - m.M;
        if(remain_rows * cs >= remain_cols * rs){
            remain_cols = cs;
        }else{
            remain_rows = rs;
        }
    }
}


template<>
void Tile<RealDevice>::process(uint64_t curCycle, AccReq arq, uint16_t offset){
    if(arq.type == Map){
        float* ml = getmLat(arq.lineAddr);
        if(ml==NULL)
        {
            Memristor_CU<RealDevice> *mcu;
            // printf("M = %d  K = %d  cells = %d   Address = %llX   mode == WMI? %d \n", arq.K, arq.M, cells, arq.lineAddr, arq.mode == WMI);
            mcu = new Memristor_CU<RealDevice>(arq.M, arq.K, reinterpret_cast<float *>(arq.lineAddr), arq.mode, bits, row_size, col_size, cells);
            add(arq.lineAddr,mcu);
            mcu->ReadLatency(cells * bits);
            addcLat(arq.lineAddr, mcu->getLatency());
            // printf("latency = %.4e", mcu->getLatency());
            ml = new float[arq.K];
            for(int i = 0; i < arq.K; ++i)
                ml[i] = mcu->getWriteLatency(i);
            addmLat(arq.lineAddr, ml);
            delete mcu;
        }
        arrayCycle = curCycle + ml[offset] * 3e9;
        mapCycle += ml[offset] * 3e9;

    }else{
        while(getcLat(arq.otherAddr) == 0.)
            usleep(10);
        arrayCycle = curCycle + getcLat(arq.otherAddr) * 3e9;
    }

}

template <>
Tile<RealDevice>::Tile(){
    PEs = new std::map<Address,Memristor_CU<RealDevice>*>();
    cLat = new std::map<Address, float>();
    mLat = new  std::map<Address, float*>();
    mtTable = new std::map<Address,uint64_t>();
    req_queue = new queue<AccReq>();
    XB_capacity = new std::vector<uint32_t>(16);
    XB_vector = new std::vector<XB>(16);
    XB_vector_premap = new std::vector<XB>(16);
    for(int i = 0; i < XB_capacity->size(); i++){
        (*XB_capacity)[i] = row_size * col_size;
        XB xb(row_size, col_size, 0);
        (*XB_vector)[i] = xb;
        (*XB_vector_premap)[i] = xb;
    }
    memCycle = 0;
    arrayCycle = 0;
    mapCycle = 0;
    tileCycles = 0;
    mpCycle = 0;
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

    futex_init(&query_mtx);
    alu = new ALU(buffer);
    control_unit = new Control_Unit(buffer,alu);
    occupation = 0;
}
template <>
void Tile<RealDevice>::prefetch(){
    futex_lock(&mtx);
    bool a = req_queue->empty();
    futex_unlock(&mtx);
    if(a){
        usleep(100);
        return;
    }
    futex_lock(&process_mtx);

    futex_lock(&mtx);
    if(req_queue->empty())
    {
        futex_unlock(&mtx);
        futex_unlock(&process_mtx);
        return;
    }
    AccReq  accReq = req_queue->front();
    req_queue->pop();
    futex_unlock(&mtx);
    if(accReq.M == 0)
        printf("error M = %d  %d \n", accReq.M, req_queue->size());
    uint32_t lineSize = buffer->getLineSize();
    int16_t cnt = 0;
    if(accReq.mode == WMI){
        if(accReq.type == Map){

            //First lock (Writer)
            //futex_lock(&filterLock);
            //Fetch Weight
            uint64_t mapTemp = arrayCycle;
            for(int i = 0; i < accReq.K; i++){
                int idx = buffer->findAvailable(memCycle);

                if(idx == -1){
                    memCycle = buffer->wait(idx);
                    count ++;
                    // idx = buffer->getUselessIdx();

                    buffer->setFlags(accReq.lineAddr,idx,Map,memCycle);
                }
                uint64_t temp = memCycle;
                for (int j = 0; j < accReq.M; j++) {
                    //TODO: Access Memory
                    //Access all the columns at the i row
                    access_count ++;
                    control_unit->execute({load,accReq.lineAddr + j * accReq.K + i,0,0},memCycle);
                    alu->preprocess(accReq.lineAddr + j * accReq.K + i, memCycle);

                    if(++cnt % lineSize == 0){
                        cnt = 0;
                        buffer->setUseless(idx, memCycle);
                        buffer->wait(idx);
                        // idx = buffer->getUselessIdx();
                    }
                }
                readCycle += memCycle - temp;
                //if(readCycle > memCycle){
                //    printf("rdCycles = %llu   memCycle = %llu  temp = %llu   delta = %llu  \n",readCycle,memCycle,temp,memCycle-temp);
                //}
                //Total rows as one batch
                //Trigger write drivers to write cells
//                    std::cout<<"idx    "<<idx<<std::endl;
                process(MAX(memCycle,arrayCycle),accReq);
                buffer->setUseless(idx,memCycle);
                cnt = 0;

            }

            tileCycles = arrayCycle;

            futex_lock(&query_mtx);


            mtTable->insert(std::pair<Address,uint64_t>(accReq.lineAddr, arrayCycle - mapTemp));
            futex_unlock(&query_mtx);
            //futex_unlock(&filterLock);
            //Unlock
        }
        else{
            //Fetch Input
            uint64_t computeTemp = arrayCycle;
            printf("replications = %d\n", accReq.replications);
            for(int i = 0; i < accReq.M /*It means the # of columns of Input*/; i+= accReq.replications){
                int idx = buffer->findAvailable(memCycle);

                if(idx == -1){
                    memCycle = buffer->wait(idx);
                    count ++;
                    // idx = buffer->getUselessIdx();
                    buffer->setFlags(accReq.lineAddr,idx,Compute,memCycle);
                }
                uint64_t temp = memCycle;
                //printf("5 %u ",(unsigned int)pthread_self());

                //futex_lock(&filterLock);
                //buffer->setBufferLine(accReq.K);
                for (int j = 0; j < accReq.K; ++j) {
                    //TODO: Access Memory
                    //Access all the rows at the j column
                    //memCycle += buffer->access(accReq.lineAddr+(j*accReq.M+i));


                    control_unit->execute({load,accReq.lineAddr + i * accReq.K + j,0,0},memCycle);
                    alu->preprocess(accReq.lineAddr + i * accReq.K + j, memCycle);

                    if(++cnt % lineSize == 0){
                        cnt = 0;

                        buffer->setUseless(idx, memCycle);
                        buffer->wait(idx);
                    }

                }
                readCycle += memCycle - temp;
                //Total rows as an input vector to multiply with matrix
                //Trigger read drivers to read cells and other components to compose the result
//                    std::cout<<"idx    "<<idx<<std::endl;
                process(MAX(memCycle,arrayCycle),accReq);
                buffer->setUseless(idx,memCycle);
                cnt = 0;
                // write back
                tileCycles = MAX(memCycle, arrayCycle);


                for(int ii = 0; ii*lineSize < accReq.N; ii++)
                    control_unit->execute({store, 0, accReq.resultAddr + accReq.N * i + lineSize * ii, lineSize}, tileCycles);

                //lock to perform read (Reader)
                //process(MAX(memCycle,arrayCycle));
                //unlock
            }

            computeCycle += arrayCycle - computeTemp;

            //futex_unlock(&filterLock);
        }
        //printf("6 %u\n",(unsigned int)pthread_self());
    }
    else{
        if(accReq.type == Map){

            // printf("   Tile map \n");
            //First lock (Writer)

            uint64_t mapTemp = arrayCycle;

            for(int i = 0; i < accReq.K; i++){
                int idx = buffer->findAvailable(memCycle);
                if(idx == -1){
                    memCycle = buffer->wait(idx);
                    count ++;
                    // idx = buffer->findAvailable(memCycle);
                    buffer->setFlags(accReq.lineAddr,idx,Map,memCycle);
                }

                uint64_t temp = memCycle;
                ///buffer->setBufferLine(accReq.M);
                for (int j = 0; j < accReq.M; j ++) {
                    //TODO: Access Memory
                    //Access all the columns at the i row
                    access_count ++;

                    //MemReq req = {accReq.lineAddr+(/*j*accReq.K+*/i), GETS, 0, &dummyState, memCycle, &filterLock, dummyState, 0, 0};


                    control_unit->execute({load,accReq.lineAddr + i * accReq.M + j,0,0},memCycle);
                    alu->preprocess(accReq.lineAddr + i * accReq.M + j, memCycle);

                    if(++cnt % lineSize == 0){
                        
                        // printf("idx = %d \n", idx);
                        cnt = 0;
                        buffer->setUseless(idx, memCycle);
                        buffer->wait(idx);
                    }



                }
                process(MAX(memCycle,arrayCycle),accReq, i);
                buffer->setUseless(idx, memCycle);
                cnt = 0;
                readCycle += memCycle - temp;
            }
            tileCycles = arrayCycle;
            arrayCycle - mapTemp; // map overhead;
            futex_lock(&query_mtx);

            mtTable->insert(std::pair<Address,uint64_t>(accReq.lineAddr, arrayCycle - mapTemp));
            futex_unlock(&query_mtx);
        }
        else{

            //First lock (Writer)

            uint64_t computeTemp = arrayCycle;
            // printf("aaaaaa\n");
            for(int i = 0; i < accReq.N; i += accReq.replications){
                int idx = buffer->findAvailable(memCycle);

                if(idx == -1){
                    memCycle = buffer->wait(idx);
                    count ++;

                    buffer->setFlags(accReq.lineAddr,idx, Compute, memCycle);
                }
                uint64_t temp = memCycle;
                //buffer->setBufferLine(accReq.M);
                for (int j = 0; j < accReq.K; j ++) {
                    //TODO: Access Memory
                    //Access all the columns at the i row
                    access_count ++;
                    //MemReq req = {accReq.lineAddr+(/*j*accReq.K+*/i), GETS, 0, &dummyState, memCycle, &filterLock, dummyState, 0, 0};
                    control_unit->execute({load,accReq.lineAddr + i * accReq.M + j,0,0},memCycle);

                    alu->preprocess(accReq.lineAddr + i * accReq.K + j, memCycle);
                    if(++cnt % lineSize == 0){
                        cnt = 0;
                        buffer->setUseless(idx, memCycle);
                        buffer->wait(idx);
                    }
                }
                buffer->setUseless(idx, memCycle);
                cnt = 0;
                process(MAX(memCycle,arrayCycle),accReq);
                tileCycles = MAX(memCycle, arrayCycle);

                readCycle += memCycle - temp;
                // printf("memcycle = %d  temp = %d  arraycycle = %d  ctemp = %d \n", memCycle, temp, arrayCycle, computeTemp);
                for(int ii = 0; ii*lineSize < accReq.N; ii++)
                    control_unit->execute({store, 0, accReq.resultAddr + accReq.N * i + ii * lineSize, lineSize}, tileCycles);

                if(accReq.N % 10 != 0)
                {
                    //max pooling
                    for(int ii  = 0; ii < accReq.replications; ii++){
                        alu->preprocess(0, arrayCycle);
                        alu->preprocess(0, arrayCycle);
                    }
                }

                //MemReq req = {accReq.resultAddr + , }
                // printf("array cycle = %d  %d\n", arrayCycle, i);
            }

            computeCycle += arrayCycle - computeTemp;
            // printf("computeCycle = %d \n", computeCycle);
            // printf("memCycle = %d   arrayCycle = %d \n", memCycle, arrayCycle);
                //lock to perform read (Reader)
                //process(MAX(memCycle,arrayCycle));
                //unlock

        }
            //unlock



    }
    //delete[] (float*)accReq.lineAddr;



    futex_unlock(&process_mtx);

    //TODO: start Map or Compute

}

template <>
Tile<RealDevice>::Tile(Buffer *_buffer, uint16_t _bits, uint16_t _row_size, uint16_t _col_size,uint16_t _capacity, uint16_t _cells, uint32_t index):buffer(_buffer),bits(_bits),row_size(_row_size),col_size(_col_size),capacity(_capacity),cells(_cells){
    PEs = new std::map<Address,Memristor_CU<RealDevice>*>();
    cLat = new std::map<Address, float>();
    mLat = new std::map<Address, float*>();
    mtTable = new std::map<Address,uint64_t>();
    req_queue = new queue<AccReq>();
    memCycle = 0;
    arrayCycle = 0;
    mapCycle = 0;
    access_count = 0;
    readCycle = 0;
    tileCycles = 0;
    mpCycle = 0;
    computeCycle = 0;
    count = 0;
    futex_init(&filterLock);
    futex_init(&mtx);
    futex_init(&process_mtx);
    futex_init(&query_mtx);
    //capacity = 1000;
    available_capacity = _capacity;
    isEnd = false;
    alu = new ALU(buffer);
    XB_capacity = new std::vector<uint32_t>(_capacity);
    XB_vector = new std::vector<XB>(_capacity);
    XB_vector_premap = new std::vector<XB>(_capacity);
    for(int i = 0; i < XB_capacity->size(); i++){
        (*XB_capacity)[i] = row_size * col_size;
        XB xb(row_size, col_size, index);
        (*XB_vector)[i] = xb;
        (*XB_vector_premap)[i] = xb;
    }
    control_unit = new Control_Unit(buffer,alu);
    cnt = 0;
    occupation = 0;
    // printf("lineSize = %d \n", _buffer->getLineSize());
}

template <>
bool Tile<RealDevice>::premap(Matrix m, uint16_t depth){
    uint64_t ac = 0;
    for(int i = 0; i < XB_capacity->size(); ++i){
        ac += (*XB_vector_premap)[i].available_capacity(); //TODO:
//            ac += xbs[i].available_capacity;
    }
    bool is_fully_map = false;
    //printf("M.size = %d * %d \n", m.K, m.M);
    if(m.size() > ac)
        return false;
    else{
        for(int i = 0; i < XB_vector_premap->size(); ++i){
            if((*XB_vector_premap)[i].compare(m)){
                (*XB_vector_premap)[i].map(m);
                sort(XB_vector_premap->begin(), XB_vector_premap->end());
                is_fully_map = true;
                break;
            }
        }
        if(!is_fully_map){
            Matrix mm = m - (*XB_vector_premap)[XB_vector_premap->size()-1];
            //printf("M.size = %d * %d  %d * %d \n", m.K, m.M, mm.K, mm.M);
            sort(XB_vector_premap->begin(), XB_vector_premap->end());
            bool b = false;
            b = premap(m, depth+1);
            if(!mm.isNull())
                b &= premap(mm, depth+1);
            is_fully_map = b;
        }
    }
    if(depth == 0){
        if(is_fully_map){
            for(int i = 0; i < XB_vector_premap->size(); ++i){
                (*XB_vector)[i].set((*XB_vector_premap)[i]);
                //(*XB_vector)[i].print_size();
            }
        }
        else{
            for(int i = 0; i < XB_vector_premap->size(); ++i){
                //backtrack
                (*XB_vector_premap)[i].set((*XB_vector)[i]);
            }
        }
    }

    return is_fully_map;
}    