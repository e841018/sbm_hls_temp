#ifndef __EXCH_H__
#define __EXCH_H__

// switch between QUBO2Ising and QUBO2Ising_3
#define Q2I3

// problem size
const int n = 5; // number of currencies
const int N = n * n; // number of variables
const int N_step = 400; // number of steps
const int M = 2; // number of substeps in a step

// degree of parallelism
const int Pb = 1; // number of blocks
const int Nb = N / Pb; // the number of spins to be handled in a block

/*
These parameters are set according to p.6 of [2].
The corresponding components in [1] are annotated in comments.
QUBO2Ising_3 is based on the objective function in [3].
[1] Combinatorial optimization by simulating adiabatic bifurcations in nonlinear Hamiltonian systems
[2] FPGA-based Simulated Bifurcation Machine
[3] A Currency Arbitrage Machine based on the Simulated Bifurcation Algorithm for Ultrafast Detection of Optimal Opportunity
*/
const float alpha0 = 1.; // Delta in [1]
const float beta0 = 1.; // K in [1]
const float Delta_t = 0.9;
const float delta_t = Delta_t / M;

// constants in objective function
const float M1 = 0.2;
const float M2 = 0.2;
const float m_c = 1100; // [3]

// initialization
const float x_init_max = 0.1;
const float p_init_max = x_init_max;

void top(float xrate[n][n], float J[N][N], float h[N]);
void QUBO2Ising(float xrate[n][n], float J[N][N], float h[N]);
void QUBO2Ising_3(float xrate[n][n], float J[N][N], float h[N]);

#endif
