#include "sbm.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <cassert>
#include <stdlib.h>

const int n_rep = 1;
const int n = 10;

void rand_init(float arr[N], float low, float high) {
    for (int i = 0; i < N; i++) {
        float r = ((float) rand()) / ((float) RAND_MAX);
        arr[i] = low + r * (high - low);
    }
}

void print_activation(bool activation[n][n]) {
    std::cout << '\n';
    std::cout << "activation =\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << ' ' << activation[i][j];
        }
        std::cout << '\n';
    }
}

void evalCost(float J[N][N], float h[N], bool spin[N]) {
    int cost = 0;
    int cut = 0;
    for (int j=1; j<N; ++j) {
        for (int i=0; i<j; ++i) {
            cost -= J[i][j] * (2*spin[i]-1) * (2*spin[j]-1);
            //cost -= J[i][j] * (spin[i] * spin[j]);
            if (spin[i] != spin[j]) {
                //std::cout << i << " " << j << " " << J[i][j] << std::endl;
                cut += J[i][j];
            }
        }
    }
    std::cout << "spin = ";
    for (int i=0; i<N; ++i) {
        std::cout << spin[i];
    }
    std::cout << "\n";
    std::cout << "cut = " << cut << std::endl;
    std::cout << "cost = " << cost << std::endl;
}

int main(int argc, char *argv[]) {
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " <input graph>\n";
        exit(EXIT_FAILURE);
    }
    std::string filename(argv[1]);
    std::ifstream ifs(filename.c_str());
    if (!ifs) {
        std::cerr << "Error: \"" << filename << "\" does not exist!!" << std::endl;
        return false;
    }

    std::cout << std::fixed << std::setprecision(4);

    float J[N][N] = {0};
    float h[N] = {0};
    int numvar, numline;
    int bufferInt1, bufferInt2, bufferInt3;
    ifs >> numvar >> numline;
    //assert(bufferInt1 == N);
    int size = 100 * 99 / 2;
    for (int i=0; i<numline; ++i) {
        ifs >> bufferInt1 >> bufferInt2 >> bufferInt3;
        //std::cerr << bufferInt1 << " " << bufferInt2 << " " << bufferInt3 << std::endl;
        if (bufferInt2 < 101) {
            J[bufferInt1-1][bufferInt2-1] -= bufferInt3;
        }
    }

    bool activation[n][n] = {0};
    float x_init[N] = {0};
    float p_init[N] = {0};    

    for (int rep=0; rep < n_rep; rep++) {
        // invoke kernel
        rand_init(x_init, -x_init_max, x_init_max);
        rand_init(p_init, -p_init_max, p_init_max);
        
        bool spin[N] = {0};
        top(J, h, x_init, p_init, spin);
        evalCost(J, h, spin);
        // reshape
        //
        //for (int i = 0; i < n; i++) {
        //    for (int j = 0; j < n; j++) {
        //        activation[i][j] = spin[i * n + j];
        //    }
        //}
        //
        //print_activation(activation);
        
        // print_objectives(xrate, activation);
        // brute_force(xrate);
    }

    return 0;
}

/*
TODOs:
 * set eta properly
 * tune parameters
Improvements:
 * is double buffer implemented correctly?
 * generate J[N][N] at compile time
*/
