#include "qubo2exch.h"
/*
Take in the spin that is solved by SBM, and output order operations

Parameters:
 * spin: spin solved by SBM
 * operations: orderEntryOperations, ends with (unsigned)(-1)
*/

void exch_out(bool spin[N],orderEntryOperation_t operations[N]){
     // output operations
    int count = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (spin[i * n + j]) {
                operations[count].timestamp = 0;
                operations[count].symbolRow = i;
                operations[count].symbolCol = j;
                count++;
            }
        }
    }
     // mark end
    if (count < max_op) operations[count].timestamp = (unsigned int)(-1);

}
