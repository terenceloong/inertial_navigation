#ifndef _TRANSFER_ALIGNMENT_H
#define _TRANSFER_ALIGNMENT_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

#define ALIGNMENT_MODLE     5

#if ALIGNMENT_MODLE == 1
    #define N_STATE     9
    #define N_MEASURE   6
#elif ALIGNMENT_MODLE == 3
    #define N_STATE     15
    #define N_MEASURE   9
#elif ALIGNMENT_MODLE == 5
    #define N_STATE     18
    #define N_MEASURE   9
#endif

typedef struct
{
    double bias[6];
    double assembly[3];
    double X[N_STATE];
    double P[N_STATE][N_STATE];
    double Q[N_STATE];
    double R[N_MEASURE];
}ST_TA_FILTER;

void initialize_ta_filter(ST_TA_FILTER *filter, double *P, double *Q, double *R);
uint8_t transfer_alignment(double qvp[10], double acc[3], double *measure, ST_TA_FILTER *filter, double dt);

#ifdef __cplusplus
}
#endif

#endif
