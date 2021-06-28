#include "exch.h"
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
 * activation: output array, activation of links in xrate
*/
void top(float xrate[n][n], float J[N][N], float h[N]) {
#ifdef Q2I3
    QUBO2Ising_3(xrate, J, h);
#else
    QUBO2Ising(xrate, J, h);
#endif
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
Convert the QUBO formulation into an Ising model, with object function from [3]

Parameters:
 * xrate: input array, log2(exchange rate)
 * J: output array, coefficient of (s_i * s_j) in Ising model
 * h: output array, coefficient of s_i in Ising model
*/
void QUBO2Ising_3(float xrate[n][n], float J[N][N], float h[N]) {
    // QUBO formulation
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
