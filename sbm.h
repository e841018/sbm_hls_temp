#ifndef __SBM_H__
#define __SBM_H__

// problem size
const int n = 3; // number of currencies
const int N = n * n; // number of variables
const int N_step = 40; // number of steps
const int M = 2; // number of substeps in a step

// degree of parallelism
const int Pb = 1; // number of blocks
const int Nb = N / Pb; // the number of spins to be handled in a block

/*
These parameters are set according to p.6 of [2].
The corresponding components in [1] are annotated in comments.
[1] Combinatorial optimization by simulating adiabatic bifurcations in nonlinear Hamiltonian systems
[2] FPGA-based Simulated Bifurcation Machine
*/
const float alpha0 = 1.; // Delta in [1]
const float beta0 = 1.; // K in [1]
const float Delta_t = 0.9;
const float delta_t = Delta_t / M;
const float alpha_init = 0.; // p(t) in [1], changes linearly from 0 to 1
const float Delta_alpha = (1. - alpha_init) / N_step;

const float M1 = 0.2;
const float M2 = 0.2;

class MMTE {
public:
    float J_sub[Nb][N];
    float h_sub[Nb];
    float x_sub[Nb];
    float p_sub[Nb];
    float Delta_t_times_gamma0;
    float alpha = alpha_init;
    float eta = alpha * 0.7 * alpha0 / n; // stated to be dependent on alpha(t) in [2], not presents in [1]
    void init(float J_sub[Nb][N], float h_sub[Nb], float x_sub_init[Nb], float p_sub_init[Nb], float Delta_t_times_gamma0);
    void step(float x_prime_in[N], float x_prime_out[Nb]);
    inline void MM(float x_prime_in[N], float Delta_P[Nb]);
    inline void TE();
    inline float FX(float x_sub_i, float h_sub_i);
    inline float FP(float p_sub_i);
};

void top(float xrate[n][n], float x_init[N], float p_init[N], bool activation[n][n]);
void QUBO2Ising(float xrate[n][n], float J[N][N], float h[N]);
float sd(float J[N][N]);
void SBM(MMTE blocks[Pb], bool spin[N]);

#endif
