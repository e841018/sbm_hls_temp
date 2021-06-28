#ifndef __SBM_H__
#define __SBM_H__

// problem size
const int N = 1000; // number of variables
const int N_step = 10; // number of steps
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

class MMTE {
public:
    float J_sub[Nb][N];
    float h_sub[Nb];
    float x_sub[Nb];
    float p_sub[Nb];
    float Delta_t_times_gamma0;
    float alpha;
    float eta;
    void init(float J_sub[Nb][N], float h_sub[Nb], float x_sub_init[Nb], float p_sub_init[Nb], float Delta_t_times_gamma0);
    void step(int t, float x_prime_in[N], float x_prime_out[Nb]);
    void update_alpha_eta(int t);
    inline void MM(float x_prime_in[N], float Delta_P[Nb]);
    inline void TE();
    inline float FX(float x_sub_i, float h_sub_i);
    inline float FP(float p_sub_i);
};

void top(float J[N][N], float h[N], float x_init[N], float p_init[N], bool spin[N]);
float sd(float J[N][N]);
void SBM(float J[N][N], float h[N], MMTE blocks[Pb], bool spin[N]);
#ifdef K2000
void evalCost(float J[N][N], float h[N], bool spin[N], int& cut, int& cost);
#endif // K2000

#endif
