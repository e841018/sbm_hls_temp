# include "sbm.h"
# include "math.h"

// maybe you can use these?
float alpha = alpha_init;
float eta = alpha; // stated to be dependent on alpha(t) in [2], not presents in [1]

/*
Time evolution

Parameters:
 * step: input scalar
 * h_sub: input array
 * x_sub: output array
 * p_sub: output array
*/
void TE(int step, float h_sub[Nb], float x_sub[Nb], float p_sub[Nb]) {
    for (int m = 0; m < M; m++) {
        float FX = delta_t * (-(alpha0 - alpha) * x_sub[m] - beta0 * pow(x_sub[m], 3) - eta * h_sub[m]);
        float FP = delta_t * p_sub[m];
        p_sub[m] += FX;
        x_sub[m] += FP;
    }
}
