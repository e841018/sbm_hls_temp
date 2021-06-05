# include "sbm.h"
# include <stdio.h>

int main(int argc, char *argv[]) {
    float xrate[n][n] = {
        {0, -0.1203764797, 0.04217729502, -2.017242085},
        {0.1203365239, 0, 0.1620771323, -1.897297122},
        {-0.042122476, -0.1619456788, 0, -2.059374059},
        {2.017276612, 1.897566294, 2.059483515, 0},
    };
    bool activation [n][n] = {0};

    top(xrate, activation);

    std::cout << "activation =\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << '\t' << activation[i][j];
        }
        std::cout << std::endl;
    }

    return 0;
}

