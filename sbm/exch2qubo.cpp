#include "exch2qubo.h"

#include "mm_datatype.h"
/*
Take in the maket data update and transform it into QUBO problem

Parameters:
 * orderBookResponse: input update market data
 * x_init: input init x value
 * p_init: input init p value
 * operations: output order operations
p_init

*/


void exch_in(orderBookResponse_t* update,float J[N][N], float h[N]){
#pragma HLS INTERFACE m_axi port=update
#pragma HLS INTERFACE m_axi port=h
#pragma HLS INTERFACE m_axi port=J
#pragma HLS INTERFACE s_axilite port=h
#pragma HLS INTERFACE s_axilite port=J
#pragma HLS INTERFACE s_axilite port=update



    // static matrix keeping track of exchange rate
    static float xrate[n][n] = {0};
    // update xrate
    //xrate[update->symbolRow][update->symbolCol] = update->askPrice;
    //printf("%u %u %f\n", update.symbolRow,update.symbolCol, update.askPrice);
    xrate[update->symbolRow][update->symbolCol] = update->askPrice;

    // convert xrate to Ising model
    QUBO2Ising(xrate, J, h);

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

