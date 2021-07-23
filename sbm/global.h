
#pragma once
#define CUR 5                         // number of currencies
#define physical_bits (CUR - 1) * CUR // total physical bits

typedef struct
{
    float price;
    int col, row;
} orderBookResponse_t;
