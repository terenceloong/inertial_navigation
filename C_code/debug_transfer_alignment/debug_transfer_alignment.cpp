// debug_transfer_alignment.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <stdio.h>
#include <stdint.h>
#include "inertial_solve.h"
#include "transfer_alignment.h"

int main()
{
    FILE *imu_file, *measure_file;
    double imu[7], measure[10];

    double qvp[10], dt, P[N_STATE], Q[N_STATE], R[N_MEASURE], acc1[3];
    ST_TA_FILTER ta_filter;
    uint8_t state;
    int k = 0;

    dt = 0.02;
    qvp[0] = 1.0;
    qvp[1] = -1.24900090270330e-16;
    qvp[2] = 8.32667268468867e-17;
    qvp[3] = 1.66533453693773e-16;
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

    fopen_s(&imu_file, "../../imu.bin", "rb");
    fopen_s(&measure_file, "../../measure.bin", "rb");

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
            if(!transfer_alignment(qvp, acc1, &measure[1], &ta_filter, dt))
                return 1;
        }
        acc1[0] = imu[4];
        acc1[1] = imu[5];
        acc1[2] = imu[6];
    }

    fclose(imu_file);
    fclose(measure_file);

    return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
