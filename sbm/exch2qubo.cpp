#include <math.h>
#include "exch2qubo.h"

//#define VERBOSE

#ifdef VERBOSE
#include <iostream>
using namespace std;
#endif

#ifdef VERBOSE
// print the matrix
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
#endif
/***************
 * Magic formulation function
 * Don't know how it works yet
 **************/
template <int size>
inline void f1(int i, int j, int k, float m[size][size], float penalty)
{
    int a1, aa;
    a1 = aa = 0;
    for (int l = 0; l < CUR; l++)
    {
        if (j != l && k != l)
        {
            if (l > k)
                aa = k + 1;
            else
                aa = k;
            a1 = l * (CUR - 1) + aa;
            m[i - 1][a1 - 1] = m[i - 1][a1 - 1] + penalty;
        }
    }
    return;
}
template <int size>
inline void f2(int i, int j, float m[size][size], float penalty)
{
    int a1, aa;
    a1 = aa = 0;
    for (int l = 0; l < CUR; l++)
    {
        if (j != l)
        {
            if (j > l)
                aa = l + 1;
            else if (j < l)
                aa = l;
            a1 = j * (CUR - 1) + aa;
            m[i - 1][a1 - 1] = m[i - 1][a1 - 1] + penalty;
        }
    }
    return;
}
template <int size>
inline void f21(int i, int j, float m[size][size], float penalty)
{
    int aa, a1;
    aa = a1 = 0;
    for (int l = 0; l < CUR; l++)
    {
        if (l != j)
        {
            if (l > j)
                aa = j + 1;
            else
                aa = j;
            a1 = (CUR - 1) * l + aa;
            m[i - 1][a1 - 1] = m[i - 1][a1 - 1] + penalty;
        }
    }
}
template <int size>
inline void f3(int i, int j, float m[size][size], float penalty)
{
    int aa, a1;
    aa = a1 = 0;
    for (int l = 0; l < CUR; l++)
    {
        if (j != l)
        {
            if (j > l)
                aa = l + 1;
            else
                aa = l;
            a1 = (CUR - 1) * j + aa;
            m[i - 1][a1 - 1] = m[i - 1][a1 - 1] - penalty;
        }
    }
    return;
}

template <int size>
inline void f31(int i, int j, float m[size][size], float penalty)
{
    int aa, a1;
    aa = a1 = 0;
    for (int l = 0; l < CUR; l++)
    {
        if (l != j)
        {
            if (l > j)
                aa = j + 1;
            else
                aa = j;
            a1 = (CUR - 1) * l + aa;
            m[i - 1][a1 - 1] = m[i - 1][a1 - 1] - penalty;
        }
    }
    return;
}

template <int size>
void exch2qubo(float p[CUR][CUR], float Q_Matrix[size][size],
               float penalty = 10, float cc = 1)
{

#ifdef VERBOSE
    cout << "matrix after log:" << endl;
    print_matrix(p);
#endif

    int m1, m2 = 0;
    // Start formulation
    // Don't know how
    // TO BE CLEARIFY
    for (int i = 0; i < CUR; i++)
    {
        int aa = 0;
        for (int j = 0; j < CUR; j++)
        {
            if (i != j)
            {
                if (i > j)
                    aa = j + 1;
                else
                    aa = j;
                m2 = (CUR - 1) * i + aa;
                f1(m2, i, j, Q_Matrix, penalty);
            }
        }
    }
    for (int i = 0; i < CUR; i++)
    {
        int aa = 0;
        for (int j = 0; j < CUR; j++)
        {
            if (i != j)
            {
                if (i < j)
                    aa = i + 1;
                else
                    aa = i;
                m2 = (CUR - 1) * j + aa;
                f1(m2, j, i, Q_Matrix, penalty);
            }
        }
    }
    for (int i = 0; i < CUR; i++)
    {
        int aa = 0;
        for (int j = 0; j < CUR; j++)
        {
            if (i != j)
            {
                if (i > j)
                    aa = j + 1;
                else
                    aa = j;
                m2 = (CUR - 1) * i + aa;
                f2(m2, i, Q_Matrix, penalty);
                f31(m2, i, Q_Matrix, penalty);
            }
        }
        for (int j = 0; j < CUR; j++)
        {
            if (i != j)
            {
                if (i < j)
                    aa = i + 1;
                else
                    aa = i;
                m2 = (CUR - 1) * j + aa;
                f21(m2, i, Q_Matrix, penalty);
                f3(m2, i, Q_Matrix, penalty);
            }
        }
    }
    for (int i = 0; i < CUR; i++)
    {
        int aa = 0;
        for (int j = 0; j < CUR; j++)
        {
            if (i != j)
            {
                if (i < j)
                    aa = i + 1;
                else
                    aa = i;
                m2 = (CUR - 1) * j + aa;
                Q_Matrix[m2 - 1][m2 - 1] =
                    Q_Matrix[m2 - 1][m2 - 1] - cc * p[j][i];
            }
        }
    }
    for (int i = 0; i < CUR; i++)
    {
        int aa = 0;
        int aa1 = 0;
        for (int j = 0; j < CUR; j++)
        {
            if (i != j)
            {
                if (i < j)
                {
                    aa = i + 1;
                    aa1 = j;
                }
                else
                {
                    aa = i;
                    aa1 = j + 1;
                }
                m1 = (CUR - 1) * i + aa1;
                m2 = (CUR - 1) * j + aa;
                Q_Matrix[m2 - 1][m1 - 1] = Q_Matrix[m2 - 1][m1 - 1] + penalty;
            }
        }
    }
    // end magic formulation

    // diagonal to off-diagonal matrix element
    // put diagonal elements to ancilla bits
    for (int i = 0; i < physical_bits; i++)
    {
        Q_Matrix[i][physical_bits] = Q_Matrix[i][i];
        Q_Matrix[physical_bits][i] = Q_Matrix[i][i];
        Q_Matrix[i][i] = 0;
    }
#ifdef VERBOSE
    cout << "Q_Matrix:" << endl;
    print_matrix(Q_Matrix);
#endif
}

void ERM(orderBookResponse_t *update, float Q_Matrix[physical_bits + 1][physical_bits + 1], float penalty, float cc, bool logged_10)
{
    // exchange rate matrix
    // 0 in diagnal
    // static float p[CUR][CUR] = {{0, 1.3194, 0.90745, 104.05, 6.72585},
    //                             {0.75799, 0, 0.68853, 78.94, 5.10327},
    //                             {1.10185, 1.45193, 0, 114.65, 7.41088},
    //                             {0.00961, 0.01266, 0.00872, 0, 0.06463},
    //                             {0.14864, 0.19586, 0.13488, 15.47, 0}};
    static float p[CUR][CUR] = {0};
    // take log10
    if (!logged_10 && update->price != 0)
    {
        p[update->row][update->col] = log10(update->price);
    }
    else
    {
        p[update->row][update->col] = update->price;
    }
    // clear previous Q_Matrix
    for (int i = 0; i < physical_bits + 1; i++)
    {
        for (int j = 0; j < physical_bits + 1; j++)
        {
            Q_Matrix[i][j] = 0;
        }
    }
    exch2qubo(p, Q_Matrix, penalty, cc);
}
