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

#ifndef CORE_H_
#define CORE_H_

#include <stdint.h>
#include "decoder.h"
#include "g_std/g_string.h"
#include "stats.h"
#include "pin.H"

struct BblInfo {
    uint32_t instrs;
    uint32_t bytes;
    DynBbl oooBbl[0]; //0 bytes, but will be 1-sized when we have an element (and that element has variable size as well)
};

enum ArithmeticIns{
    ADD = XED_ICLASS_ADD,
    ADC = XED_ICLASS_ADC,
    ADCX = XED_ICLASS_ADCX,
    ADDPD = XED_ICLASS_ADDPD,
    ADDPS = XED_ICLASS_ADDPS,
    ADDSD = XED_ICLASS_ADDSD,
    ADDSS = XED_ICLASS_ADDSS,
    ADDSUBPD = XED_ICLASS_ADDSUBPD,
    ADDSUBPS = XED_ICLASS_ADDSUBPS,
    ADOX = XED_ICLASS_ADOX,
    DEC = XED_ICLASS_DEC,
    DIV = XED_ICLASS_DIV,
    DIVPD = XED_ICLASS_DIVPD,
    DIVPS = XED_ICLASS_DIVPS,
    DIVSD = XED_ICLASS_DIVSD,
    DIVSS = XED_ICLASS_DIVSS,
    DPPD = XED_ICLASS_DPPD,
    DPPS = XED_ICLASS_DPPS,
    FADD = XED_ICLASS_FADD,
    FADDP = XED_ICLASS_FADDP,
    FDIV = XED_ICLASS_FDIV,
    FIADD = XED_ICLASS_FIADD,
    FIMUL = XED_ICLASS_FIMUL,
    FMUL = XED_ICLASS_FMUL,
    FMULP = XED_ICLASS_FMULP,
    FSUB = XED_ICLASS_FSUB,
    FSUBP = XED_ICLASS_FSUBP,
    FSUBR = XED_ICLASS_FSUBR,
    FSUBRP = XED_ICLASS_FSUBRP,
    HADDPD = XED_ICLASS_HADDPD,
    HADDPS = XED_ICLASS_HADDPS,
    HSUBPD = XED_ICLASS_HSUBPD,
    HSUBPS = XED_ICLASS_HSUBPS,
    IDIV = XED_ICLASS_IDIV,
    IMUL = XED_ICLASS_IMUL,
    INC = XED_ICLASS_INC,
    MUL = XED_ICLASS_MUL,
    MULPD = XED_ICLASS_MULPD,
    MULPS = XED_ICLASS_MULPS,
    MULSD = XED_ICLASS_MULSD,
    MULSS = XED_ICLASS_MULSS,
    MULX = XED_ICLASS_MULX,
    PADDB = XED_ICLASS_PADDB,
    PADDD = XED_ICLASS_PADDD,
    PADDQ = XED_ICLASS_PADDQ,
    PADDSB = XED_ICLASS_PADDSB,
    PADDSW = XED_ICLASS_PADDSW,
    PADDUSB = XED_ICLASS_PADDUSB,
    PADDUSW = XED_ICLASS_PADDUSW,
    PADDW = XED_ICLASS_PADDW,
    PHADDD = XED_ICLASS_PHADDD,
    PHADDSW = XED_ICLASS_PHADDSW,
    PHADDW = XED_ICLASS_PHADDW,
    PHSUBD = XED_ICLASS_PHSUBD,
    PHSUBSW = XED_ICLASS_PHSUBSW,
    PHSUBW = XED_ICLASS_PHSUBW,
    PMADDUBSW = XED_ICLASS_PMADDUBSW,
    PMADDWD = XED_ICLASS_PMADDWD,
    PMULDQ = XED_ICLASS_PMULDQ,
    PMULHRSW = XED_ICLASS_PMULHRSW,
    PMULHUW = XED_ICLASS_PMULHUW,
    PMULHW = XED_ICLASS_PMULHW,
    PMULLD = XED_ICLASS_PMULLD,
    PMULLW = XED_ICLASS_PMULLW,
    PMULUDQ = XED_ICLASS_PMULUDQ,
    PSUBB = XED_ICLASS_PSUBB,
    PSUBD = XED_ICLASS_PSUBD,
    PSUBQ = XED_ICLASS_PSUBQ,
    PSUBSB = XED_ICLASS_PSUBSB,
    PSUBSW = XED_ICLASS_PSUBSW,
    PSUBUSB = XED_ICLASS_PSUBUSB,
    PSUBUSW = XED_ICLASS_PSUBUSW,
    PSUBW = XED_ICLASS_PSUBW,
    SBB = XED_ICLASS_SBB,
    SUB = XED_ICLASS_SUB,
    SUBPD = XED_ICLASS_SUBPD,
    SUBPS = XED_ICLASS_SUBPS,
    SUBSD = XED_ICLASS_SUBSD,
    SUBSS = XED_ICLASS_SUBSS
};
/* Analysis function pointer struct
 * As an artifact of having a shared code cache, we need these to be the same for different core types.
 */
struct InstrFuncPtrs {  // NOLINT(whitespace)
    void (*loadPtr)(THREADID, ADDRINT);
    void (*storePtr)(THREADID, ADDRINT);
    void (*bblPtr)(THREADID, ADDRINT, BblInfo*);
    void (*branchPtr)(THREADID, ADDRINT, BOOL, ADDRINT, ADDRINT);
    // Same as load/store functions, but last arg indicated whether op is executing
    void (*predLoadPtr)(THREADID, ADDRINT, BOOL);
    void (*predStorePtr)(THREADID, ADDRINT, BOOL);

    void (*arithmeticPtr)(THREADID,ArithmeticIns);
    uint64_t type;
    uint64_t pad[1];
    //NOTE: By having the struct be a power of 2 bytes, indirect calls are simpler (w/ gcc 4.4 -O3, 6->5 instructions, and those instructions are simpler)
};


//TODO: Switch type to an enum by using sizeof macros...
#define FPTR_ANALYSIS (0L)
#define FPTR_JOIN (1L)
#define FPTR_NOP (2L)
#define FPTR_RETRY (3L)

//Generic core class

class Core : public GlobAlloc {
    private:
        uint64_t lastUpdateCycles;
        uint64_t lastUpdateInstrs;

    protected:
        g_string name;

    public:
        explicit Core(g_string& _name) : lastUpdateCycles(0), lastUpdateInstrs(0), name(_name) {}

        virtual uint64_t getInstrs() const = 0; // typically used to find out termination conditions or dumps
        virtual uint64_t getPhaseCycles() const = 0; // used by RDTSC faking --- we need to know how far along we are in the phase, but not the total number of phases
        virtual uint64_t getCycles() const = 0;

        virtual void initStats(AggregateStat* parentStat) = 0;
        virtual void contextSwitch(int32_t gid) = 0; //gid == -1 means descheduled, otherwise this is the new gid

        //Called by scheduler on every leave and join action, before barrier methods are called
        virtual void leave() {}
        virtual void join() {}

        virtual InstrFuncPtrs GetFuncPtrs() = 0;
};

#endif  // CORE_H_

