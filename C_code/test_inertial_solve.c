#include "mex.h"
#include "inertial_solve.h"
#include <string.h>

extern uint8_t flag; //##########

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    size_t Am;
    int Nimu, Nnav;
    double *imu_data, *qvp0, *h, *nav_result;
    int i, k;
    double imu[6], qvp[10], dt;
    double angle[3];
    uint8_t state;

    /* Size */
    Am = mxGetM(prhs[0]);
    Nimu = Am;
    if(Nimu%2 == 1)
        Nnav = (Nimu-1)/2 + 1;
    else
        Nnav = Nimu/2 + 1;

    /* Apply for output space */
    plhs[0] = mxCreateDoubleMatrix(Nnav, 10, mxREAL);

    /* Input */
    imu_data = mxGetDoubles(prhs[0]);
    qvp0 = mxGetDoubles(prhs[1]);
    h = mxGetDoubles(prhs[2]);

    /* Output */
    nav_result = mxGetDoubles(plhs[0]);

    /* Initialize */
    memcpy(qvp, qvp0, 80);

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

        // state = inertial_slove_rate(qvp, imu, dt);
        state = inertial_slove_delta(qvp, imu, dt);

        if(state)
        {
            // qvp[9] = h[k]; //height locking

            quat2angle(angle, qvp);
            nav_result[k] = k*dt; //t, s
            nav_result[  Nnav+k] = qvp[7] * 57.295779513082323; //lat, deg
            nav_result[2*Nnav+k] = qvp[8] * 57.295779513082323; //lon, deg
            nav_result[3*Nnav+k] = qvp[9]; //h, m
            nav_result[4*Nnav+k] = qvp[4]; //vn, m/s
            nav_result[5*Nnav+k] = qvp[5]; //ve, m/s
            nav_result[6*Nnav+k] = qvp[6]; //vd, m/s
            nav_result[7*Nnav+k] = angle[0] * 57.295779513082323; //yaw, deg
            nav_result[8*Nnav+k] = angle[1] * 57.295779513082323; //pitch, deg
            nav_result[9*Nnav+k] = angle[2] * 57.295779513082323; //roll, deg

            k++;
        }
    }
}
