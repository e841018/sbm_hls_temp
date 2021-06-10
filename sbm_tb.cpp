#include "sbm.h"
#include <iostream>
#include <iomanip>
#include <stdlib.h>

void rand_init(float arr[N], float low, float high) {
    for (int i = 0; i < N; i++) {
        float r = ((float) rand()) / ((float) RAND_MAX);
        arr[i] = low + r * (high - low);
    }
}

float calc_obj0(float xrate[n][n], bool activation[n][n]) {
    float obj0 = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            obj0 += activation[i][j] * xrate[i][j];
    return obj0;
}

float calc_obj1(bool activation[n][n]) {
    float obj1 = 0;
    for (int i = 0; i < n; i++) {
        int temp = 0;
        for (int j = 0; j < n; j++) {
            temp += activation[i][j];
            temp -= activation[j][i];
        }
        obj1 -= M1 * temp * temp;
    }
    return obj1;
}

float calc_obj2(bool activation[n][n]) {
    float obj2 = 0;
    for (int i = 0; i < n; i++) {
        int temp = 0;
        for (int j = 0; j < n; j++) {
            temp += activation[i][j];
        }
        obj2 -= M2 * temp * (temp - 1);
    }
    return obj2;
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

void print_objectives(float xrate[n][n], bool activation[n][n]) {
    float obj0 = calc_obj0(xrate, activation);
    float obj1 = calc_obj1(activation);
    float obj2 = calc_obj2(activation);
    std::cout << '\n';
    std::cout << "maximize {  obj0  +   obj1  +   obj2 }" << '\n';
    std::cout << "          " << std::setw(7) << obj0
                        << " - " << std::setw(7) << - obj1
                        << " - " << std::setw(7) << - obj2
                        << " = " << obj0 + obj1 + obj2 << '\n';
}

bool find_1_on_diagonal(bool activation[n][n]) {
    for (int i = 0; i < n; i++)
        if (activation[i][i] == 1)
            return true;
    return false;
}

void brute_force(float xrate[n][n]) {
    int n_config = 1 << N;
    bool activation[n][n];
    for (int config = 0; config < n_config; config++) {
        int bin = config;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                activation[i][j] = bin & 1;
                bin >>= 1;
            }
        }

        if (find_1_on_diagonal(activation)) continue;
        float obj0 = calc_obj0(xrate, activation);
        if (obj0 <= 0) continue;
        if (calc_obj1(activation) < 0) continue;
        if (calc_obj2(activation) < 0) continue;

        std::cout << '\n';
        std::cout << "brute_force: solution found: \n";
        print_activation(activation);
        std::cout << "obj0 = " << obj0 << '\n';
    }
}

int main(int argc, char *argv[]) {
    std::cout << std::fixed << std::setprecision(4);

    float xrate[n][n] = {
        {0, -0.1203764797, 0.04217729502},
        {0.1203365239, 0, 0.1620771323},
        {-0.042122476, -0.1619456788, 0},
    };

    bool activation[n][n] = {0};
    float x_init[N] = {0};
    float p_init[N] = {0};

    for (int rep=0; rep < 1; rep++) {
        // invoke kernel
        rand_init(x_init, -0.1, 0.1);
        rand_init(p_init, -0.1, 0.1);
        top(xrate, x_init, p_init, activation);
        print_activation(activation);
        print_objectives(xrate, activation);
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
