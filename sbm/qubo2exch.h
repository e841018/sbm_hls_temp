#ifndef __QUBO2EXCH_H__
#define __QUBO2EXCH_H__

#include "mm_datatype.h"
#include "global.h"
extern "C" {
void exch_out(bool spin[N],orderEntryOperation_t operations[N]);
}

#endif
