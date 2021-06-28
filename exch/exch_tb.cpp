#include "exch.h"
#include <iostream>
#include <iomanip>
#include <stdlib.h>

float calc_obj0(float xrate[n][n], bool activation[n][n]) {
    float obj0 = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            obj0 += activation[i][j] * xrate[i][j];
    return obj0;
}

float calc_obj0_3(float xrate[n][n], bool activation[n][n]) {
    return calc_obj0(xrate, activation) * m_c;
}

float calc_obj1(bool activation[n][n]) {
    float obj1 = 0;
    for (int i = 0; i < n; i++) {
        int temp = 0;
        for (int j = 0; j < n; j++) {
            temp += activation[i][j];
            temp -= activation[j][i];
        }
        obj1 -= temp * temp;
    }
    return M1 * obj1;
}

float calc_obj1_3(bool activation[n][n]) {
    return calc_obj1(activation) / M1;
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

float calc_obj2_3(bool activation[n][n]) {
    float obj2 = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int jp = 0; jp < j; jp++)
                obj2 -= activation[i][j] * activation[i][jp] + activation[j][i] * activation[jp][i];
    return obj2;
}

float calc_obj3_3(bool activation[n][n]) {
    float obj3 = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            obj3 -= activation[i][j] * activation[j][i];
    return obj3;
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
#ifdef Q2I3
    float obj0 = calc_obj0_3(xrate, activation);
    float obj1 = calc_obj1_3(activation);
    float obj2 = calc_obj2_3(activation);
    float obj3 = calc_obj3_3(activation);
    std::cout << '\n';
    std::cout << "maximize {  obj0  +   obj1  +   obj2  +   obj3 }" << '\n';
    std::cout << "          " << std::setw(7) << obj0
                        << " - " << std::setw(7) << (obj1 == 0 ? 0. : - obj1)
                        << " - " << std::setw(7) << (obj2 == 0 ? 0. : - obj2)
                        << " - " << std::setw(7) << (obj3 == 0 ? 0. : - obj3)
                        << " = " << obj0 + obj1 + obj2 << '\n';
#else
    float obj0 = calc_obj0(xrate, activation);
    float obj1 = calc_obj1(activation);
    float obj2 = calc_obj2(activation);
    std::cout << '\n';
    std::cout << "maximize {  obj0  +   obj1  +   obj2 }" << '\n';
    std::cout << "          " << std::setw(7) << obj0
                        << " - " << std::setw(7) << (obj1 == 0 ? 0. : - obj1)
                        << " - " << std::setw(7) << (obj2 == 0 ? 0. : - obj2)
                        << " = " << obj0 + obj1 + obj2 << '\n';
#endif
}

int main(int argc, char *argv[]) {
    std::cout << std::fixed << std::setprecision(4);

    //float xrate[n][n] = {
    //    {0, -0.1203764797, 0.04217729502},
    //    {0.1203365239, 0, 0.1620771323},
    //    {-0.042122476, -0.1619456788, 0},
    //};

    float xrate[n][n] = {
        {0,	-0.1203764797,	0.04217729502,	-2.017242085,	-0.8277471774},
        {0.1203365239,	0,	0.1620771323,	-1.897297122,	-0.7078485463},
        {-0.042122476,	-0.1619456788,	0.	-2.059374059,	-0.869869781},
        {2.017276612,	1.897566294,	2.059483515,	0,	1.189565844},
        {0.8278643034,	0.7080542498,	0.8700524427,	-1.189490314,	0}
    };

    float J[N][N] = {0};
    float h[N] = {0};
    // invoke kernel
    top(xrate, J, h);

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
