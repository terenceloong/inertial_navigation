#include "mex.h"
#include "inertial_solve.h"
#include "transfer_alignment.h"
#include <string.h>
#include <math.h>

#define R2D     57.295779513082323

extern uint8_t flag; //##########

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    size_t Am;
    int Nimu, Nnav;
    double *imu_data, *measure_data, *qvp0, *P, *Q, *R; //input
    double *nav_result, *bias_esti, *assembly_esti, *filter_P; //output
    double imu[6], measure[9], qvp[10], dt, acc1[3];
    ST_TA_FILTER ta_filter;
    double angle[3];
    uint8_t state;
    int i, k, j;

    /* Size */
    Am = mxGetM(prhs[0]);
    Nimu = Am;
    Nnav = (Nimu-1)/2 + 1;

    /* Apply for output space */
    plhs[0] = mxCreateDoubleMatrix(Nnav, 10, mxREAL); //nav_result
    plhs[1] = mxCreateDoubleMatrix(Nnav, 7, mxREAL);  //bias_esti
    plhs[2] = mxCreateDoubleMatrix(Nnav, 4, mxREAL);  //assembly_esti
    plhs[3] = mxCreateDoubleMatrix(Nnav, sizeof(ta_filter.X)/8, mxREAL); //filter_P

    /* Input */
    imu_data = mxGetDoubles(prhs[0]);
    measure_data = mxGetDoubles(prhs[1]);
    qvp0 = mxGetDoubles(prhs[2]);
    P = mxGetDoubles(prhs[3]);
    Q = mxGetDoubles(prhs[4]);
    R = mxGetDoubles(prhs[5]);

    /* Output */
    nav_result = mxGetDoubles(plhs[0]);
    bias_esti = mxGetDoubles(plhs[1]);
    assembly_esti = mxGetDoubles(plhs[2]);
    filter_P = mxGetDoubles(plhs[3]);

    /* Initialize */
    memcpy(qvp, qvp0, sizeof(qvp));
    initialize_ta_filter(&ta_filter, P, Q, R);

    dt = (imu_data[1]-imu_data[0]) * 2.0;

    flag = 0; //##########
    k = 0;
    for(i=0; i<Nimu; i++)
    {
        imu[0] = imu_data[  Nimu+i];
        imu[1] = imu_data[2*Nimu+i];
        imu[2] = imu_data[3*Nimu+i];
        imu[3] = imu_data[4*Nimu+i];
        imu[4] = imu_data[5*Nimu+i];
        imu[5] = imu_data[6*Nimu+i];

        imu_compensate(imu, ta_filter.bias);
        state = inertial_slove_rate(qvp, imu, dt);

        if(state==1 && k>0)
        {
            measure[0] = measure_data[  Nnav+k];
            measure[1] = measure_data[2*Nnav+k];
            measure[2] = measure_data[3*Nnav+k];
            measure[3] = measure_data[4*Nnav+k];
            measure[4] = measure_data[5*Nnav+k];
            measure[5] = measure_data[6*Nnav+k];
            measure[6] = measure_data[7*Nnav+k];
            measure[7] = measure_data[8*Nnav+k];
            measure[8] = measure_data[9*Nnav+k];
            transfer_alignment(qvp, acc1, measure, &ta_filter, dt);
        }
        acc1[0] = imu[3];
        acc1[1] = imu[4];
        acc1[2] = imu[5];

        if(state)
        {
            quat2angle(angle, qvp);
            nav_result[       k] = k*dt; //t, s
            nav_result[  Nnav+k] = qvp[7] * R2D; //lat, deg
            nav_result[2*Nnav+k] = qvp[8] * R2D; //lon, deg
            nav_result[3*Nnav+k] = qvp[9]; //h, m
            nav_result[4*Nnav+k] = qvp[4]; //vn, m/s
            nav_result[5*Nnav+k] = qvp[5]; //ve, m/s
            nav_result[6*Nnav+k] = qvp[6]; //vd, m/s
            nav_result[7*Nnav+k] = angle[0] * R2D; //yaw, deg
            nav_result[8*Nnav+k] = angle[1] * R2D; //pitch, deg
            nav_result[9*Nnav+k] = angle[2] * R2D; //roll, deg

            bias_esti[       k] = k*dt;
            bias_esti[  Nnav+k] = ta_filter.bias[0] * R2D;
            bias_esti[2*Nnav+k] = ta_filter.bias[1] * R2D;
            bias_esti[3*Nnav+k] = ta_filter.bias[2] * R2D;
            bias_esti[4*Nnav+k] = ta_filter.bias[3];
            bias_esti[5*Nnav+k] = ta_filter.bias[4];
            bias_esti[6*Nnav+k] = ta_filter.bias[5];

            assembly_esti[       k] = k*dt;
            assembly_esti[  Nnav+k] = ta_filter.assembly[0] * R2D;
            assembly_esti[2*Nnav+k] = ta_filter.assembly[1] * R2D;
            assembly_esti[3*Nnav+k] = ta_filter.assembly[2] * R2D;

            for(j=0; j<sizeof(ta_filter.X)/8; j++)
                filter_P[j*Nnav+k] = sqrt(ta_filter.P[j][j]);

            k++;
        }
    }
}
