/** $lic$
 * Copyright (C) 2012-2015 by Massachusetts Institute of Technology
 * Copyright (C) 2010-2013 by The Board of Trustees of Stanford University
 *
 * This file is part of zsim.
 *
 * zsim is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, version 2.
 *
 * If you use this software in your research, we request that you reference
 * the zsim paper ("ZSim: Fast and Accurate Microarchitectural Simulation of
 * Thousand-Core Systems", Sanchez and Kozyrakis, ISCA-40, June 2013) as the
 * source of the simulator in any publications that use this software, and that
 * you send us a citation of your work.
 *
 * zsim is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "simple_core.h"
#include "filter_cache.h"
#include "zsim.h"

SimpleCore::SimpleCore(FilterCache* _l1i, FilterCache* _l1d, g_string& _name) : Core(_name), l1i(_l1i), l1d(_l1d), instrs(0), curCycle(0), haltedCycles(0) {
}

void SimpleCore::initStats(AggregateStat* parentStat) {
    AggregateStat* coreStat = new AggregateStat();
    coreStat->init(name.c_str(), "Core stats");
    auto x = [this]() -> uint64_t { assert(curCycle >= haltedCycles); return curCycle - haltedCycles; };
    auto cyclesStat = makeLambdaStat(x);
    cyclesStat->init("cycles", "Simulated cycles");
    ProxyStat* instrsStat = new ProxyStat();
    instrsStat->init("instrs", "Simulated instructions", &instrs);
    coreStat->append(cyclesStat);
    coreStat->append(instrsStat);
    parentStat->append(coreStat);
}

uint64_t SimpleCore::getPhaseCycles() const {
    return curCycle % zinfo->phaseLength;
}

void SimpleCore::load(Address addr) {
    curCycle = l1d->load(addr, curCycle);
}

void SimpleCore::store(Address addr) {
    curCycle = l1d->store(addr, curCycle);
}

void SimpleCore::bbl(Address bblAddr, BblInfo* bblInfo) {
    //info("BBL %s %p", name.c_str(), bblInfo);
    //info("%d %d", bblInfo->instrs, bblInfo->bytes);
    instrs += bblInfo->instrs;
    curCycle += bblInfo->instrs;

    Address endBblAddr = bblAddr + bblInfo->bytes;
    for (Address fetchAddr = bblAddr; fetchAddr < endBblAddr; fetchAddr+=(1 << lineBits)) {
        curCycle = l1i->load(fetchAddr, curCycle);
    }
}

void SimpleCore::contextSwitch(int32_t gid) {
    if (gid == -1) {
        l1i->contextSwitch();
        l1d->contextSwitch();
    }
}

void SimpleCore::join() {
    //info("[%s] Joining, curCycle %ld phaseEnd %ld haltedCycles %ld", name.c_str(), curCycle, phaseEndCycle, haltedCycles);
    if (curCycle < zinfo->globPhaseCycles) { //carry up to the beginning of the phase
        haltedCycles += (zinfo->globPhaseCycles - curCycle);
        curCycle = zinfo->globPhaseCycles;
    }
    phaseEndCycle = zinfo->globPhaseCycles + zinfo->phaseLength;
    //note that with long events, curCycle can be arbitrarily larger than phaseEndCycle; however, it must be aligned in current phase
    //info("[%s] Joined, curCycle %ld phaseEnd %ld haltedCycles %ld", name.c_str(), curCycle, phaseEndCycle, haltedCycles);
}


//Static class functions: Function pointers and trampolines

InstrFuncPtrs SimpleCore::GetFuncPtrs() {
    return {LoadFunc, StoreFunc, BblFunc, BranchFunc, PredLoadFunc, PredStoreFunc, ArithmeticFunc, FPTR_ANALYSIS, {0}};
}

void SimpleCore::LoadFunc(THREADID tid, ADDRINT addr) {
    static_cast<SimpleCore*>(cores[tid])->load(addr);
}

void SimpleCore::StoreFunc(THREADID tid, ADDRINT addr) {
    static_cast<SimpleCore*>(cores[tid])->store(addr);
}

void SimpleCore::PredLoadFunc(THREADID tid, ADDRINT addr, BOOL pred) {
    if (pred) static_cast<SimpleCore*>(cores[tid])->load(addr);
}

void SimpleCore::PredStoreFunc(THREADID tid, ADDRINT addr, BOOL pred) {
    if (pred) static_cast<SimpleCore*>(cores[tid])->store(addr);
}

void SimpleCore::BblFunc(THREADID tid, ADDRINT bblAddr, BblInfo* bblInfo) {
    SimpleCore* core = static_cast<SimpleCore*>(cores[tid]);
    core->bbl(bblAddr, bblInfo);

    while (core->curCycle > core->phaseEndCycle) {
        assert(core->phaseEndCycle == zinfo->globPhaseCycles + zinfo->phaseLength);
        core->phaseEndCycle += zinfo->phaseLength;

        uint32_t cid = getCid(tid);
        //NOTE: TakeBarrier may take ownership of the core, and so it will be used by some other thread. If TakeBarrier context-switches us,
        //the *only* safe option is to return inmmediately after we detect this, or we can race and corrupt core state. If newCid == cid,
        //we're not at risk of racing, even if we were switched out and then switched in.
        uint32_t newCid = TakeBarrier(tid, cid);
        if (newCid != cid) break; /*context-switch*/
    }
}

//Add Arithmetic Simulation
void SimpleCore::compute(ArithmeticIns type){
    uint32_t cycles = 0;
    switch(type){
         case ADD:
         case ADC:
         case ADCX:cycles = 1;break;
         case ADDPD:cycles = 7;break;
         case ADDPS:
         case ADDSD:
         case ADDSS:cycles = 5;break;
         case ADDSUBPD:cycles = 7;break;
         case ADDSUBPS:cycles = 5;break;
         case ADOX:
         case DEC:cycles = 1;break;
         case DIV:cycles = 35;break;
         case DIVPD:cycles = 125;break;
         case DIVPS:cycles = 70;break;
         case DIVSD:cycles = 62;break;
         case DIVSS:cycles = 34;break;
         case DPPD:cycles = 12;break;
         case DPPS:cycles = 15;break;
         case FADD:
         case FADDP:cycles = 5;break;
         case FDIV:cycles = 39;break;
         case FIADD:
         case FIMUL:cycles = 11;break;
         case FMUL:
         case FMULP:
         case FSUB:
         case FSUBP:
         case FSUBR:
         case FSUBRP:cycles = 5;break;
         case HADDPD:
         case HADDPS:
         case HSUBPD:
         case HSUBPS:cycles = 9;break;
         case IDIV:cycles = 33;break;
         case IMUL:cycles = 7;break;
         case INC:cycles = 1;break;
         case MUL:cycles = 7;break;
         case MULPD:
         case MULPS:
         case MULSD:cycles = 5;break;
         case MULSS:
         case MULX:cycles = 4;break;
         case PADDB:
         case PADDD:
         case PADDQ:
         case PADDSB:
         case PADDSW:
         case PADDUSB:
         case PADDUSW:
         case PADDW:cycles = 1;break;
         case PHADDD:cycles = 4;break;
         case PHADDSW:cycles = 8;break;
         case PHADDW:cycles = 6;break;
         case PHSUBD:cycles = 4;break;
         case PHSUBSW:cycles = 8;break;
         case PHSUBW:cycles = 6;break;
         case PMADDUBSW:
         case PMADDWD:cycles = 4;break;
         case PMULDQ:cycles = 5;break;
         case PMULHRSW:
         case PMULHUW:
         case PMULHW:cycles = 4;break;
         case PMULLD:cycles = 11;break;
         case PMULLW:cycles = 4;break;
         case PMULUDQ:cycles = 5;break;
         case PSUBB:
         case PSUBD:
         case PSUBQ:
         case PSUBSB:
         case PSUBSW:
         case PSUBUSB:
         case PSUBUSW:
         case PSUBW:
         case SBB:
         case SUB:cycles = 1;break;
         case SUBPD:cycles = 7;break;
         case SUBPS:
         case SUBSD:
         case SUBSS:cycles = 5;break;
         default: cycles = 0;
    }
    curCycle += cycles;
}

void SimpleCore::ArithmeticFunc(THREADID tid, ArithmeticIns type){
    static_cast<SimpleCore*>(cores[tid])->compute(type);
}


