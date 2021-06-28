#ifndef __EXCH_H__
#define __EXCH_H__

#include "mm_datatype.h"
#include "sbm.h"

#define Q2I3 // Convert the QUBO formulation into an Ising model, with object function from [3]

const int n = 5;
const int n_seq = N; // number of xrate updates
const int max_op = N; // output buffer size for storing operations
void exch(orderBookResponse_t update, orderEntryOperation_t operations[N], float x_init[N], float p_init[N]);
void QUBO2Ising(float xrate[n][n], float J[N][N], float h[N]);

#endif