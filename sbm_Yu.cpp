#include "sbm.h"
#include <stdlib.h>

/*
Top-level function to export

Parameters:
 * xrate: input array, log2(exchange rate)
    example: xrate[1][3] = log2(2.5) means 1 currency1 = 2.5 currency3
 * activation: output array, activation of links in xrate
*/
void top(float xrate[n][n], bool activation[n][n]) {
    int J[N][N] = {0};
    float h[N] = {0};
    bool spin[N] = {0};
    QUBO2Ising(xrate, J, h);
    SBM(J, h, spin);
    // reshape
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            activation[i][j] = spin[i * n + j];
        }
    }
}

/*
Convert the QUBO formulation into an Ising model

Parameters:
 * xrate: input array, log2(exchange rate)
 * J: output array, coefficient of Ising model
 * h: output array, coefficient of Ising model
*/
void QUBO2Ising(float xrate[n][n], int J[N][N], float h[N]) {
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
    // convert to Ising model (J -> J and h)
    for (int i = 0; i < N; i++) {
        int temp = 0;
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

/*
Simulated bifurcation machine for finding an approximate solution of an Ising model of the form:
  \min_s \left( - \frac{1}{2} \sum_{i=0}^{N-1} \sum_{j=0}^{N-1} {J_{ij} s_i s_j} + \sum_{i=0}^{N-1} {h_i s_i} \right)

Parameters:
 * J: input array, coefficient of (s_i * s_j)
 * h: input array, coefficient of s_i
 * spin: output array, true: s_i=1, false: s_i=-1
*/
void SBM(int J[N][N], float h[N], bool spin[N]) {
    float x_a[N] = {0}; // double buffer
    float x_b[N] = {0}; // double buffer
    for (int step = 0; step < N_step; step++) {
        for (int block=0; block < Pb; block++) {
            int offset = Nb * block; // [Nb * block : Nb * block + Nb]
            if (step % 2 == 0) {
                MMTE(step, J + offset, h + offset, x_a, x_b + offset);
            } else {
                MMTE(step, J + offset, h + offset, x_b, x_a + offset);
            }
        }
    }
    
    for (int i = 0; i < N; i++)
        spin[i] = x_a[i] > 0; // spin = sign(x)
}

/*
Simulate a subset of spins, which is a block of size Nb

Parameters:
 * step: input scalar
 * J_sub: input array, passed to MM
 * h_sub: input array, passed to TE
 * x_prime_in: input array, x_prime from all blocks, passed to MM
 * x_prime_out: output array, x_prime of this block
*/
void MMTE(int step, int J_sub[Nb][N], float h_sub[Nb], float x_prime_in[N], float x_prime_out[Nb]) {
    static float x_sub[Nb];
    static float p_sub[Nb];
    float Delta_P[Nb] = {0};
    
    if (step == 0) {
        rand_init(x_sub);
        rand_init(p_sub);
    }

    MM(step, J_sub, x_prime_in, Delta_P);
    for (int i = 0; i < Nb; i++) {
        p_sub[i] += Delta_P[i];
    }
    TE(step, h_sub, x_sub, p_sub);
    for (int i = 0; i < Nb; i++) {
        x_prime_out[i] = x_sub[i] * Delta_t_times_gamma0;
    }
}

/*
Matrix-vector multiplication, Delta_P = J_sub * x

Parameters:
 * step: input scalar
 * J_sub: input array, cached if step == 0
 * x_prime_in: input array
 * Delta_P: output array
*/
void MM(int step, int J_sub[Nb][N], float x_prime_in[N], float Delta_P[Nb]) {
    static int J_sub_local[Nb][N];
    if (step == 0) {
        for (int i = 0; i < Nb; i++) {
            for (int j = 0; j < N; j++) {
                J_sub_local[i][j] = J_sub[i][j];
            }
        }
    }

    for (int i = 0; i < Nb; i++) {
        int temp = 0;
        for (int j = 0; j < N; j++) {
            temp += J_sub_local[i][j] * x_prime_in[i];
        }
        Delta_P[i] = temp;
    }
}

/*
Initialize with random numbers in [-0.1, 0.1)

Parameters:
 * arr: output array
*/
void rand_init(float arr[Nb]) {
    // temporary implementation on CPU
    for (int i = 0; i < Nb; i++) {
        float r = ((float) rand()) / ((float) RAND_MAX);
        arr[i] = -0.1 + 0.2 * r;
    }
}