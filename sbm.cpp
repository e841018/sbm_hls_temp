#include "sbm.h"
#include <math.h> // sqrt

#define DEBUG_PRINT
#define LOG_X_P

#ifdef DEBUG_PRINT
#include <iostream>
#include <iomanip>

#ifdef LOG_X_P
#include <fstream>
std::ofstream log_file;
#endif

void print_J(float J[N][N]) {
    std::cout << '\n';

    std::cout << "J =          (";
    for (int j = 0; j < N; j++)
        std::cout << std::setw(2) << j / n << ',' << std::setw(2) << j % n << ") (";
    std::cout << "\b \n";

    std::cout << "           +-";
    for (int j = 0; j < N; j++)
        std::cout << "--------";
    std::cout << '\n';

    for (int i = 0; i < N; i++) {
        std::cout << "    (" << std::setw(2) << i / n << ',' << std::setw(2) << i % n << ")| ";
        for (int j = 0; j < N; j++) {
            if (J[i][j] == 0)
                std::cout << std::setw(7) << 0 << " ";
            else
                std::cout << std::setw(7) << J[i][j] << " ";
        }
        std::cout << '\n';
    }
}

void print_h(float h[N]) {
    std::cout << '\n';

    std::cout << "             (";
    for (int j = 0; j < N; j++)
        std::cout << std::setw(2) << j / n << ',' << std::setw(2) << j % n << ") (";
    std::cout << "\b \n";

    std::cout << "h =          ";
    for (int i = 0; i < N; i++)
        std::cout << std::setw(7) << h[i] << " ";
    std::cout << '\n';
}

std::ostream &output() {
#ifdef LOG_X_P
    return log_file;
#else
    return std::cout;
#endif
}

void print_x_p(float x[Nb], float p[Nb]) {
    output() << "x | ";
    for (int i = 0; i < Nb; i++)
        output() << std::setw(7) << x[i] << " | ";
    output() << '\n';

    output() << "p | ";
    for (int i = 0; i < Nb; i++)
        output() << std::setw(7) << p[i] << " | ";
    output() << '\n';
}
#endif

/*
Top-level function to export

Parameters:
 * xrate: input array, log2(exchange rate)
    example: xrate[1][3] = log2(2.5) means 1 currency1 = 2.5 currency3
 * x_init: input array, used to initialize x[N]
 * p_init: input array, used to initialize p[N]
 * activation: output array, activation of links in xrate
*/
void top(float xrate[n][n], float x_init[N], float p_init[N], bool activation[n][n]) {
    float J[N][N] = {0};
    float h[N] = {0};
    bool spin[N] = {0};
    MMTE blocks[Pb];
    QUBO2Ising(xrate, J, h);
    float std_of_J  = sd(J);
    float gamma0 = 0.7 * alpha0 / std_of_J / n; // Xi0 in [1]
    float Delta_t_times_gamma0 = Delta_t * gamma0;
    for (int b=0; b < Pb; b++) {
        int offset = Nb * b; // [Nb * b : Nb * b + Nb]
        blocks[b].init(J + offset, h + offset, x_init + offset, p_init + offset, Delta_t_times_gamma0);
    }
    SBM(blocks, spin);

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
 * J: output array, coefficient of (s_i * s_j) in Ising model
 * h: output array, coefficient of s_i in Ising model
*/
void QUBO2Ising(float xrate[n][n], float J[N][N], float h[N]) {
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

#ifdef DEBUG_PRINT
    print_J(J);
    print_h(h);
#endif
}

/*
calculate population standard deviation

Parameters:
 * J: input array

 Returns:
 * std: output scalar
*/
float sd(float J[N][N]) {
    float mean = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            mean += J[i][j];
        }
    }
    mean /= (N * N);

    float var = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            float diff = J[i][j] - mean;
            var += diff * diff;
        }
    }
    var /= (N * N);

    float std = sqrt(var);
    
#ifdef DEBUG_PRINT
    std::cout << '\n';
    std::cout << "sd(J) = " << std << '\n';
#endif

    return std;
}

/*
Simulated bifurcation machine for finding an approximate solution of an Ising model of the form:
  \min_s \left( - \frac{1}{2} \sum_{i=0}^{N-1} \sum_{j=0}^{N-1} {J_{ij} s_i s_j} + \sum_{i=0}^{N-1} {h_i s_i} \right)

Parameters:
 * blocks: object array, elements are initialized MMTE modules
 * spin: output array, true: s_i=1, false: s_i=-1
*/
void SBM(MMTE blocks[Pb], bool spin[N]) {
    float x_a[N] = {0}; // double buffer
    float x_b[N] = {0}; // double buffer

#ifdef DEBUG_PRINT
#ifdef LOG_X_P
    log_file.open("log.txt");
#else
    std::cout << '\n';
    std::cout << "t=0 (init):\n";
#endif
    for (int b=0; b < Pb; b++)
        print_x_p(blocks[b].x_sub, blocks[b].p_sub);
#endif

    for (int t = 0; t < N_step; t++) {

#ifdef DEBUG_PRINT
#ifndef LOG_X_P
        std::cout << "t=" << t+1 << ":\n";
#endif
#endif

        for (int b=0; b < Pb; b++) {
            int offset = Nb * b; // [Nb * b : Nb * b + Nb]
            if (t % 2 == 0) {
                blocks[b].step(x_a, x_b + offset);
            } else {
                blocks[b].step(x_b, x_a + offset);
            }
        }
    }

#ifdef DEBUG_PRINT
#ifdef LOG_X_P
    log_file.close();
#endif
#endif
    
    for (int i = 0; i < N; i++)
        spin[i] = x_a[i] > 0; // spin = sign(x)
}


/*
Initialize a MMTE block (Nb spins)

Parameters:
 * J_sub: input array
 * h_sub: input array
 * x_sub_init: input array
 * p_sub_init: input array
 * Delta_t_times_gamma0: input scalar
*/
void MMTE::init(float J_sub[Nb][N], float h_sub[Nb], float x_sub_init[Nb], float p_sub_init[Nb], float Delta_t_times_gamma0) {
    for (int i = 0; i < Nb; i++) {
        for (int j = 0; j < N; j++) {
            this->J_sub[i][j] = J_sub[i][j];
        }
        this->h_sub[i] = h_sub[i];
        this->x_sub[i] = x_sub_init[i];
        this->p_sub[i] = p_sub_init[i];
        this->Delta_t_times_gamma0 = Delta_t_times_gamma0;
    }
}

/*
Simulate a time step

Parameters:
 * x_prime_in: input array, x_prime from all blocks, passed to MM
 * x_prime_out: output array, x_prime of this block
*/
void MMTE::step(float x_prime_in[N], float x_prime_out[Nb]) {
    float Delta_P[Nb] = {0};

    MM(x_prime_in, Delta_P);

    for (int i = 0; i < Nb; i++) {
        p_sub[i] += Delta_P[i];
    }

    TE();

    for (int i = 0; i < Nb; i++) {
        x_prime_out[i] = x_sub[i] * Delta_t_times_gamma0;
    }

    alpha += Delta_alpha;
    eta = alpha * 0.7 * alpha0 / n;

#ifdef DEBUG_PRINT
    print_x_p(x_sub, p_sub);
#endif
}

/*
Matrix-vector multiplication, Delta_P = J_sub * x

Parameters:
 * x_prime_in: input array
 * Delta_P: output array
*/
inline void MMTE::MM(float x_prime_in[N], float Delta_P[Nb]) {
    for (int i = 0; i < Nb; i++) {
        float temp = 0;
        for (int j = 0; j < N; j++) {
            temp += J_sub[i][j] * x_prime_in[i];
        }
        Delta_P[i] = temp;
    }
}

/*
Time evolution, updates x_sub and p_sub
*/
inline void MMTE::TE() {
    for (int m = 0; m < M; m++) {
        for (int i = 0; i < Nb; i++) {
            p_sub[i] += FX(x_sub[i], h_sub[i]);
            x_sub[i] += FP(p_sub[i]);
        }
    }
}

/*
FX of aSB

Parameters:
 * x_sub_i: input scalar
 * h_sub_i: input scalar

Returns:
 * delta_p: output scalar
*/
inline float MMTE::FX(float x_sub_i, float h_sub_i) {
    float delta_p =  delta_t * (
        - (alpha0 - alpha + beta0 * x_sub_i * x_sub_i) * x_sub_i
        - eta * h_sub_i
    );
    return delta_p;
}

/*
FP

Parameters:
 * p_sub_i: input scalar

Returns:
 * delta_x: output scalar
*/
inline float MMTE::FP(float p_sub_i) {
    float delta_x = delta_t * p_sub_i;
    return delta_x;
}