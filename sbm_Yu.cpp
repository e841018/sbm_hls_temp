# include "sbm.h"

/*
Simulated bifurcation machine for finding an approximate solution of an Ising problem of the form:
\min_s -1/2 \sum_{i=0}^{N-1} \sum_{j=0}^{N-1} {J_{ij} s_i s_j} + \sum_{i=0}^{N-1} {h_i s_i}

Parameters:
    J: input array, coefficient of (s_i * s_j)
    h: input array, coefficient of s_i
    spin: output array, true: s_i=1, false: s_i=-1
*/
void SBM(int J[N][N], float h[N], bool spin[N]){
    float x[N] = {0}; // x in [1]
    float p[N] = {0}; // y in [1]
    float alpha = alpha_init;
    float eta = alpha; // stated to be dependent on alpha(t) in [2], not presents in [1]

}

// TODOs:
// * calculate std_of_J
// * generate J[N][N] at compile time
// * change datatype of J[N][N]