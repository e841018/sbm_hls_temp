#pragma once
#include "global.h"
void ERM(orderBookResponse_t *update, float Q_Matrix[physical_bits + 1][physical_bits + 1], float penalty = 10, float cc = 1, bool logged_10 = false);
