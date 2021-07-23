#include <iostream>
#include "exch2qubo.h"

using namespace std;

template <int size>
void print_matrix(float m[size][size])
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            cout << m[i][j] << " ";
        }
        cout << endl;
    }
}

int main()
{
    static float p[CUR][CUR] = {{0, 1.3194, 0.90745, 104.05, 6.72585},
                                {0.75799, 0, 0.68853, 78.94, 5.10327},
                                {1.10185, 1.45193, 0, 114.65, 7.41088},
                                {0.00961, 0.01266, 0.00872, 0, 0.06463},
                                {0.14864, 0.19586, 0.13488, 15.47, 0}};
    float penalty = 10;
    float cc = 2; // or 1 in paper
    // total physical bits + ancilla bits for local field
    float Q_Matrix[physical_bits + 1][physical_bits + 1] = {0};
    orderBookResponse_t update;

    for (int i = 0; i < CUR; i++)
    {
        for (int j = 0; j < CUR; j++)
        {
            update.row = i;
            update.col = j;
            update.price = p[i][j];
            ERM(&update, Q_Matrix, penalty, cc, false);
        }
    }

    print_matrix(Q_Matrix);
    return 0;
}
