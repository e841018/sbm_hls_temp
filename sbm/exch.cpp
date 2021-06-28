#include "exch.h"
#include "mm_datatype.h"
#include <limits>

/*
Top-level function to export

Parameters:
*/
void exch(orderBookResponse_t update, orderEntryOperation_t operations[N], float x_init[N], float p_init[N]) {
#pragma HLS INTERFACE m_axi port=operations offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=x_init offset=slave bundle=gmem1
#pragma HLS INTERFACE m_axi port=p_init offset=slave bundle=gmem0

#pragma HLS INTERFACE s_axilite port=update
#pragma HLS INTERFACE s_axilite port=operations
#pragma HLS INTERFACE s_axilite port=x_init
#pragma HLS INTERFACE s_axilite port=p_init

    // static matrix keeping track of exchange rate
    static float xrate[n][n];

    // update xrate
    xrate[update.symbolRow][update.symbolRow] = update.askPrice;

    // skip first n_seq - 1 runs
    static int update_count = 0;
    update_count++;
    if (update_count < n_seq)
        return;

    // convert xrate to Ising model
    float J[N][N] = {0};
    float h[N] = {0};
    QUBO2Ising(xrate, J, h);

    // solve Ising model
    bool spin[N] = {0};
    top(J, h, x_init, p_init, spin);

    // output operations
    int count = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (spin[i * n + j]) {
                orderEntryOperation_t operation = operations[count];
                operation.timestamp = 0; // TODO
                operation.symbolRow = i;
                operation.symbolCol = j;
                count++;
            }
        }
    }

    // mark end
    orderEntryOperation_t operation = operations[count];
    operation.timestamp = (uint64_t)(-1);
}

/*
Convert the QUBO formulation into an Ising model

Parameters:
 * xrate: input array, log2(exchange rate)
 * J: output array, coefficient of (s_i * s_j) in Ising model
 * h: output array, coefficient of s_i in Ising model
*/

void QUBO2Ising(float xrate[n][n], float J[N][N], float h[N]) {

    // QUBO formulation
#ifdef Q2I3
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int ij = i * n + j;
            int ji = j * n + i;
            J[ij][ij] -= xrate[i][j] * m_c + 0.5;
            J[ji][ji] -= 0.5;
            J[ij][ji] += 1;
            for (int k = 0; k < n; k++) {
                int ik = i * n + k;
                int ki = k * n + i;
                J[ij][ik] += 1 + 0.5;
                J[ij][ki] -= 1;
                J[ji][ik] -= 1;
                J[ji][ki] += 1 + 0.5;
            }
        }
    }
#else
    // QUBO formulation
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int ij = i * n + j;
            int ji = j * n + i;
            J[ij][ij] -= xrate[i][j] + M2;
            for (int k = 0; k < n; k++) {
                int ik = i * n + k;
                int ki = k * n + i;
                J[ij][ik] += M1 + M2;
                J[ij][ki] -= M1;
                J[ji][ik] -= M1;
                J[ji][ki] += M1;
            }
        }
    }
#endif

    // convert to Ising model (J -> J and h)
    for (int i = 0; i < N; i++) {
        float temp = 0;
        for (int j = 0; j < N; j++) {
            temp += J[i][j] + J[j][i];
        }
        h[i] = temp / 4;
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            J[i][j] /= -2;
        }
    }
}
