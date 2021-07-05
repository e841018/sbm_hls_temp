#include <stdlib.h>

#include <iostream>

#include "exch2qubo.h"
#include "global.h"
#include "mm_datatype.h"

int main() {
    orderBookResponse_t currentorder;
    float J[N][N] = {0};
    float h[N] = {0};
    for (int i = 0; i < N; i++) {
        currentorder.symbolRow = i / n;
        currentorder.symbolCol = i % n;
        currentorder.askPrice = ((float)rand()) / ((float)RAND_MAX);
        exch_in(currentorder, J, h);
    }

    for (int i = 0;i<N*N;i++){
    	std::cout << J[i/N][i%N] << " ";
    }
    std::cout << std::endl;
    for (int i = 0;i<N;i++){
    	std::cout << h[i] << " ";
    }
    std::cout << std::endl;
}
