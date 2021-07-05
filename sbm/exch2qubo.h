#ifndef __EXCH2QUBO_H__
#define __EXCH2QUBO_H__

#include "mm_datatype.h"
#include "global.h"


#define Q2I3 // Convert the QUBO formulation into an Ising model, with object function from [3]


extern "C"
{
   void exch_in(orderBookResponse_t* update,float J[N][N], float h[N]);
   void QUBO2Ising(float xrate[n][n], float J[N][N], float h[N]);

}
#endif
