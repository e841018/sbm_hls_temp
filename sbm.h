#ifndef __SBM_H__
#define __SBM_H__

#include <iostream>

// problem size
const int n = 4; // number of currencies
const int N = n * n; // number of variables
const int N_step = 40; // number of steps
const int M = 2; // number of substeps in a step

// degree of parallelism
const int Pb = 2; // number of blocks
const int Nb = N / Pb; // the number of spins to be handled in a block
const int Pr = 2; // row size of J submatrix
const int Pc = 2; // column size of J submatrix

// These parameters are set according to p.6 of [2].
// The corresponding components in [1] are annotated in comments.
// [1] Combinatorial optimization by simulating adiabatic bifurcations in nonlinear Hamiltonian systems
// [2] FPGA-based Simulated Bifurcation Machine
const float alpha0 = 1.; // Delta in [1]
const float beta0 = 1.; // K in [1]
const float std_of_J = 1.; // mentioned in [1], haven't calculate yet
const float gamma0 = 0.7 * alpha0 / std_of_J / n; // Xi0 in [1]
const float Delta_t = 0.9;
const float delta_t = Delta_t / M;
const float alpha_init = 0.; // p(t) in [1], changes linearly from 0 to 1
const float Delta_alpha = (1. - alpha_init) / N_step;

void qubo2ising(float xrate[n][n], int J[N][N], float h[N]);
void SBM(int J[N][N], float h[N], bool spin[N]);
void MMTE(float x_prime_in[N], float x_prime_out[Nb]);
void MM(int Ji[Nb][N], float x[N], float Delta_P[Nb]);
void TE(float x_in[Nb], float p_in[Nb], float x_out[Nb], float p_out[Nb]);

#endif