# include "sbm.h"
# include "math.h"

/*
Time evolution

Parameters:
 * step: input scalar
 * h_sub: input array
 * x_sub: output array
 * p_sub: output array
*/
void TE(int step, float h_sub[Nb], float x_sub[Nb], float p_sub[Nb]) {
    static float alpha;
    static float eta; // stated to be dependent on alpha(t) in [2], not presents in [1]
    static int h_sub_local[Nb];
    if (step == 0) {
        alpha = alpha_init;
        eta = alpha;
        for (int i = 0; i < Nb; i++) {
            h_sub_local[i] = h_sub[i];
        }
    }
    else {
        alpha += Delta_alpha;
        eta = alpha;
    }

    for (int m = 0; m < M; m++) {
        float FX = delta_t * (-(alpha0 - alpha) * x_sub[m] - beta0 * pow(x_sub[m], 3) - eta * h_sub_local[m]);
        p_sub[m] += FX;
        float FP = delta_t * p_sub[m];
        x_sub[m] += FP;
    }
}
