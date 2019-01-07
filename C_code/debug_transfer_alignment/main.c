#include <stdio.h>
#include <stdint.h>
#include "inertial_solve.h"
#include "transfer_alignment.h"

int main()
{
    FILE *imu_file, *measure_file;
    double imu[7], measure[10];

    double qvp[10], dt, P[N_STATE], Q[N_STATE], R[N_MEASURE];
    ST_TA_FILTER ta_filter;
    uint8_t state;
    int k = 0;

    dt = 0.02;
    qvp[0] = 1.0;
    qvp[1] = 0.0;
    qvp[2] = 0.0;
    qvp[3] = 0.0;
    qvp[4] = 0.0;
    qvp[5] = 0.0;
    qvp[6] = 0.0;
    qvp[7] = 0.523598775598299;
    qvp[8] = 2.09439510239320;
    qvp[9] = 300;

    P[0] = 0.000304617419786709;
    P[1] = 0.000304617419786709;
    P[2] = 0.000304617419786709;
    P[3] = 1.0;
    P[4] = 1.0;
    P[5] = 1.0;
    P[6] = 6.14543064411833e-13;
    P[7] = 6.14543064411833e-13;
    P[8] = 25;

    Q[0] = 4.69946644021501e-17;
    Q[1] = 4.69946644021501e-17;
    Q[2] = 4.69946644021501e-17;
    Q[3] = 1.99939600000000e-10;
    Q[4] = 1.99939600000000e-10;
    Q[5] = 1.99939600000000e-10;
    Q[6] = 4.91485977925104e-28;
    Q[7] = 4.91485977925104e-28;
    Q[8] = 1.99939600000000e-14;

    R[0] = 0.0001;
    R[1] = 0.0001;
    R[2] = 0.0001;
    R[3] = 2.45817225764733e-16;
    R[4] = 2.45817225764733e-16;
    R[5] = 0.01;
    initialize_ta_filter(&ta_filter, P, Q, R);

    imu_file = fopen("../imu.bin", "rb");
    measure_file = fopen("../measure.bin", "rb");

    while(fread(imu, 8, 7, imu_file) == 7)
    {
        imu_compensate(&imu[1], ta_filter.bias);
        state = inertial_slove_rate(qvp, &imu[1], dt);
        if(state)
        {
            fread(measure, 8, 10, measure_file);
            k++;
        }
        if(state==1 && k>1)
        {
            transfer_alignment(qvp, &imu[4], &measure[1], &ta_filter, dt);
        }
    }

    fclose(imu_file);
    fclose(measure_file);

    return 0;
}
