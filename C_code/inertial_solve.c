#include "inertial_solve.h"
#include "axes.h"
#include "matrix_operation.h"
#include <math.h>

void imu_compensate(double imu[6], double bias[6])
{
    imu[0] -= bias[0];
    imu[1] -= bias[1];
    imu[2] -= bias[2];
    imu[3] -= bias[3];
    imu[4] -= bias[4];
    imu[5] -= bias[5];
}

uint8_t flag = 0;
/* Define 'flag' as global vaiable. In Matlab, after executing mex function, memory isn't
   cleared. When executing the mex function again, static variable can't be refreshed, it
   remain the last value, so that causing calculation error. Using directive 'clear mex'
   can slove this problem, but it will slow the execution. So it's better not to use static
   variable. In other environment, define 'flag' as static vaiable.*/

uint8_t inertial_slove_rate(double qvp[10], double imu[6], double dt)
{
    // static uint8_t flag = 0;
    static double wibb0[3], wibb1[3], acc0[3], acc1[3];
    double sin_lat, cos_lat, Rm, Rn, Rd;
    double wien[3], wenn[3], winn[3], winb[3], wnbb0[3], wnbb1[3], wnbb2[3];
    double Cnb[3][3];
    double dt_2, dt_4, dt_6;
    double q1, q2, q3, q4, q;
    double K1[4], K2[4], K3[4], K4[4];
    double dtheta1[3], dtheta2[3], dv1[3], dv2[3], dvt[3], dv[3], v0[3], dp[3];

    if(flag == 0)
    {
        wibb0[0] = imu[0];
        wibb0[1] = imu[1];
        wibb0[2] = imu[2];
        acc0[0]  = imu[3];
        acc0[1]  = imu[4];
        acc0[2]  = imu[5];
        flag = 1;
        return 1;
    }

    if(flag == 1)
    {
        wibb1[0] = imu[0];
        wibb1[1] = imu[1];
        wibb1[2] = imu[2];
        acc1[0]  = imu[3];
        acc1[1]  = imu[4];
        acc1[2]  = imu[5];
        flag = 2;
        return 0;
    }

    dt_2 = dt / 2.0;
    dt_4 = dt / 4.0;
    dt_6 = dt / 6.0;

    sin_lat = sin(qvp[7]);
    cos_lat = cos(qvp[7]);

    Rd = 1.0-(2.0-EARTH_F)*EARTH_F*sin_lat*sin_lat;
    Rm = (1.0-EARTH_F)*(1.0-EARTH_F)*EARTH_A / sqrt(Rd*Rd*Rd) + qvp[9];
    Rn = EARTH_A / sqrt(Rd) + qvp[9];

    quat2dcm(Cnb, qvp);

    wien[0] =  EARTH_W*cos_lat;
    wien[1] =  0.0;
    wien[2] = -EARTH_W*sin_lat;

    wenn[0] =  qvp[5] / Rn;
    wenn[1] = -qvp[4] / Rm;
    wenn[2] = -qvp[5] / Rn * tan(qvp[7]);

    winn[0] = wien[0] + wenn[0];
    winn[1] = wien[1] + wenn[1];
    winn[2] = wien[2] + wenn[2];

    winb[0] = Cnb[0][0]*winn[0] + Cnb[0][1]*winn[1] + Cnb[0][2]*winn[2];
    winb[1] = Cnb[1][0]*winn[0] + Cnb[1][1]*winn[1] + Cnb[1][2]*winn[2];
    winb[2] = Cnb[2][0]*winn[0] + Cnb[2][1]*winn[1] + Cnb[2][2]*winn[2];

    wnbb0[0] = wibb0[0] - winb[0];
    wnbb0[1] = wibb0[1] - winb[1];
    wnbb0[2] = wibb0[2] - winb[2];
    wnbb1[0] = wibb1[0] - winb[0];
    wnbb1[1] = wibb1[1] - winb[1];
    wnbb1[2] = wibb1[2] - winb[2];
    wnbb2[0] = imu[0]   - winb[0];
    wnbb2[1] = imu[1]   - winb[1];
    wnbb2[2] = imu[2]   - winb[2];

    /*----------------------------Attitude update----------------------------*/
    /* Differential equations of quaternions */
    q1 = qvp[0];
    q2 = qvp[1];
    q3 = qvp[2];
    q4 = qvp[3];
    K1[0] = (-wnbb0[0]*q2 - wnbb0[1]*q3 - wnbb0[2]*q4) * 0.5;
    K1[1] = ( wnbb0[0]*q1 + wnbb0[2]*q3 - wnbb0[1]*q4) * 0.5;
    K1[2] = ( wnbb0[1]*q1 - wnbb0[2]*q2 + wnbb0[0]*q4) * 0.5;
    K1[3] = ( wnbb0[2]*q1 + wnbb0[1]*q2 - wnbb0[0]*q3) * 0.5;
    q1 = qvp[0] + K1[0] * dt_2;
    q2 = qvp[1] + K1[1] * dt_2;
    q3 = qvp[2] + K1[2] * dt_2;
    q4 = qvp[3] + K1[3] * dt_2;
    K2[0] = (-wnbb1[0]*q2 - wnbb1[1]*q3 - wnbb1[2]*q4) * 0.5;
    K2[1] = ( wnbb1[0]*q1 + wnbb1[2]*q3 - wnbb1[1]*q4) * 0.5;
    K2[2] = ( wnbb1[1]*q1 - wnbb1[2]*q2 + wnbb1[0]*q4) * 0.5;
    K2[3] = ( wnbb1[2]*q1 + wnbb1[1]*q2 - wnbb1[0]*q3) * 0.5;
    q1 = qvp[0] + K2[0] * dt_2;
    q2 = qvp[1] + K2[1] * dt_2;
    q3 = qvp[2] + K2[2] * dt_2;
    q4 = qvp[3] + K2[3] * dt_2;
    K3[0] = (-wnbb1[0]*q2 - wnbb1[1]*q3 - wnbb1[2]*q4) * 0.5;
    K3[1] = ( wnbb1[0]*q1 + wnbb1[2]*q3 - wnbb1[1]*q4) * 0.5;
    K3[2] = ( wnbb1[1]*q1 - wnbb1[2]*q2 + wnbb1[0]*q4) * 0.5;
    K3[3] = ( wnbb1[2]*q1 + wnbb1[1]*q2 - wnbb1[0]*q3) * 0.5;
    q1 = qvp[0] + K3[0] * dt;
    q2 = qvp[1] + K3[1] * dt;
    q3 = qvp[2] + K3[2] * dt;
    q4 = qvp[3] + K3[3] * dt;
    K4[0] = (-wnbb2[0]*q2 - wnbb2[1]*q3 - wnbb2[2]*q4) * 0.5;
    K4[1] = ( wnbb2[0]*q1 + wnbb2[2]*q3 - wnbb2[1]*q4) * 0.5;
    K4[2] = ( wnbb2[1]*q1 - wnbb2[2]*q2 + wnbb2[0]*q4) * 0.5;
    K4[3] = ( wnbb2[2]*q1 + wnbb2[1]*q2 - wnbb2[0]*q3) * 0.5;
    q1 = qvp[0] + (K1[0] + 2.0*K2[0] + 2.0*K3[0] + K4[0]) * dt_6;
    q2 = qvp[1] + (K1[1] + 2.0*K2[1] + 2.0*K3[1] + K4[1]) * dt_6;
    q3 = qvp[2] + (K1[2] + 2.0*K2[2] + 2.0*K3[2] + K4[2]) * dt_6;
    q4 = qvp[3] + (K1[3] + 2.0*K2[3] + 2.0*K3[3] + K4[3]) * dt_6;
    /* Quaternion normalization */
    q = sqrt(q1*q1 + q2*q2 + q3*q3 + q4*q4);
    qvp[0] = q1 / q;
    qvp[1] = q2 / q;
    qvp[2] = q3 / q;
    qvp[3] = q4 / q;

    /*----------------------------Velocity update----------------------------*/
    v0[0] = qvp[4];
    v0[1] = qvp[5];
    v0[2] = qvp[6];

    dtheta1[0] = (wnbb0[0]+wnbb1[0]) * dt_4;
    dtheta1[1] = (wnbb0[1]+wnbb1[1]) * dt_4;
    dtheta1[2] = (wnbb0[2]+wnbb1[2]) * dt_4;
    dtheta2[0] = (wnbb1[0]+wnbb2[0]) * dt_4;
    dtheta2[1] = (wnbb1[1]+wnbb2[1]) * dt_4;
    dtheta2[2] = (wnbb1[2]+wnbb2[2]) * dt_4;
    dv1[0] = (acc0[0]+acc1[0]) * dt_4;
    dv1[1] = (acc0[1]+acc1[1]) * dt_4;
    dv1[2] = (acc0[2]+acc1[2]) * dt_4;
    dv2[0] = (acc1[0]+ imu[3]) * dt_4;
    dv2[1] = (acc1[1]+ imu[4]) * dt_4;
    dv2[2] = (acc1[2]+ imu[5]) * dt_4;

    dvt[0] = 0.5*(dtheta1[1]*dv1[2]-dtheta1[2]*dv1[1]) + 7.0/6.0*(dtheta1[1]*dv2[2]-dtheta1[2]*dv2[1]) + 
             0.5*(dtheta2[1]*dv2[2]-dtheta2[2]*dv2[1]) - 1.0/6.0*(dtheta2[1]*dv1[2]-dtheta2[2]*dv1[1]) + dv1[0] + dv2[0];
    dvt[1] = 0.5*(dtheta1[2]*dv1[0]-dtheta1[0]*dv1[2]) + 7.0/6.0*(dtheta1[2]*dv2[0]-dtheta1[0]*dv2[2]) + 
             0.5*(dtheta2[2]*dv2[0]-dtheta2[0]*dv2[2]) - 1.0/6.0*(dtheta2[2]*dv1[0]-dtheta2[0]*dv1[2]) + dv1[1] + dv2[1];
    dvt[2] = 0.5*(dtheta1[0]*dv1[1]-dtheta1[1]*dv1[0]) + 7.0/6.0*(dtheta1[0]*dv2[1]-dtheta1[1]*dv2[0]) + 
             0.5*(dtheta2[0]*dv2[1]-dtheta2[1]*dv2[0]) - 1.0/6.0*(dtheta2[0]*dv1[1]-dtheta2[1]*dv1[0]) + dv1[2] + dv2[2];
    
    winn[0] += wien[0]; // winn[0] = 2.0*wien[0] + wenn[0];
    winn[1] += wien[1]; // winn[1] = 2.0*wien[1] + wenn[1];
    winn[2] += wien[2]; // winn[2] = 2.0*wien[2] + wenn[2];

    dv[0] = Cnb[0][0]*dvt[0] + Cnb[1][0]*dvt[1] + Cnb[2][0]*dvt[2] - dt*(winn[1]*v0[2]-winn[2]*v0[1]);
    dv[1] = Cnb[0][1]*dvt[0] + Cnb[1][1]*dvt[1] + Cnb[2][1]*dvt[2] - dt*(winn[2]*v0[0]-winn[0]*v0[2]);
    dv[2] = Cnb[0][2]*dvt[0] + Cnb[1][2]*dvt[1] + Cnb[2][2]*dvt[2] - dt*(winn[0]*v0[1]-winn[1]*v0[0]) + dt*9.8;

    // if((dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2]) > 1e-20) //can be commented
    {
        qvp[4] += dv[0];
        qvp[5] += dv[1];
        qvp[6] += dv[2];
    }

    /*----------------------------Position update----------------------------*/
    dp[0] = (v0[0]+qvp[4]) * dt_2;
    dp[1] = (v0[1]+qvp[5]) * dt_2;
    dp[2] = (v0[2]+qvp[6]) * dt_2;

    // if((dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2]) > 1e-20) //can be commented
    {
        qvp[7] += dp[0] / Rm;
        qvp[8] += dp[1] / (Rn*cos_lat);
        qvp[9] -= dp[2];
    }

    wibb0[0] = imu[0];
    wibb0[1] = imu[1];
    wibb0[2] = imu[2];
    acc0[0]  = imu[3];
    acc0[1]  = imu[4];
    acc0[2]  = imu[5];

    flag= 1;
    return 1;
}

uint8_t inertial_slove_delta(double qvp[10], double imu[6], double dt)
{
    // static uint8_t flag = 0;
    static double dtheta1[3], dv1[3];
    double dtheta2[3], dv2[3];
    double dtheta1_c[3], dtheta2_c[3];
    double dt_2;
    double sin_lat, cos_lat, Rm, Rn, Rd;
    double wien[3], wenn[3], winn[3], winb[3];
    double Cnb[3][3];
    double v0[3], dvt[3], dv[3], p0[3], dp[3];
    double Phi[3], phi, qc[4], q[4], qm;
    double Cn0b[3][3], Ce0n0[3][3], Cee0[3][3], Cen[3][3], Cne[3][3], Cne0[3][3], Cnn0[3][3];

    if(flag == 0)
    {
        dtheta1[0] = imu[0];
        dtheta1[1] = imu[1];
        dtheta1[2] = imu[2];
        dv1[0]     = imu[3];
        dv1[1]     = imu[4];
        dv1[2]     = imu[5];
        flag = 2;
        return 1;
    }

    if(flag == 1)
    {
        dtheta1[0] = imu[0];
        dtheta1[1] = imu[1];
        dtheta1[2] = imu[2];
        dv1[0]     = imu[3];
        dv1[1]     = imu[4];
        dv1[2]     = imu[5];
        flag = 2;
        return 0;
    }

    dtheta2[0] = imu[0];
    dtheta2[1] = imu[1];
    dtheta2[2] = imu[2];
    dv2[0]     = imu[3];
    dv2[1]     = imu[4];
    dv2[2]     = imu[5];

    dt_2 = dt / 2.0;

    sin_lat = sin(qvp[7]);
    cos_lat = cos(qvp[7]);

    Rd = 1.0-(2.0-EARTH_F)*EARTH_F*sin_lat*sin_lat;
    Rm = (1.0-EARTH_F)*(1.0-EARTH_F)*EARTH_A / sqrt(Rd*Rd*Rd) + qvp[9];
    Rn = EARTH_A / sqrt(Rd) + qvp[9];

    quat2dcm(Cnb, qvp);

    wien[0] =  EARTH_W*cos_lat;
    wien[1] =  0.0;
    wien[2] = -EARTH_W*sin_lat;

    wenn[0] =  qvp[5] / Rn;
    wenn[1] = -qvp[4] / Rm;
    wenn[2] = -qvp[5] / Rn * tan(qvp[7]);

    winn[0] = wien[0] + wenn[0];
    winn[1] = wien[1] + wenn[1];
    winn[2] = wien[2] + wenn[2];

    winb[0] = Cnb[0][0]*winn[0] + Cnb[0][1]*winn[1] + Cnb[0][2]*winn[2];
    winb[1] = Cnb[1][0]*winn[0] + Cnb[1][1]*winn[1] + Cnb[1][2]*winn[2];
    winb[2] = Cnb[2][0]*winn[0] + Cnb[2][1]*winn[1] + Cnb[2][2]*winn[2];

    /*----------------------------Velocity update----------------------------*/
    v0[0] = qvp[4];
    v0[1] = qvp[5];
    v0[2] = qvp[6];
    
    dtheta1_c[0] = dtheta1[0] - winb[0]*dt_2;
    dtheta1_c[1] = dtheta1[1] - winb[1]*dt_2;
    dtheta1_c[2] = dtheta1[2] - winb[2]*dt_2;
    dtheta2_c[0] = dtheta2[0] - winb[0]*dt_2;
    dtheta2_c[1] = dtheta2[1] - winb[1]*dt_2;
    dtheta2_c[2] = dtheta2[2] - winb[2]*dt_2;

    dvt[0] = 0.5*(dtheta1_c[1]*dv1[2]-dtheta1_c[2]*dv1[1]) + 7.0/6.0*(dtheta1_c[1]*dv2[2]-dtheta1_c[2]*dv2[1]) + 
             0.5*(dtheta2_c[1]*dv2[2]-dtheta2_c[2]*dv2[1]) - 1.0/6.0*(dtheta2_c[1]*dv1[2]-dtheta2_c[2]*dv1[1]) + dv1[0] + dv2[0];
    dvt[1] = 0.5*(dtheta1_c[2]*dv1[0]-dtheta1_c[0]*dv1[2]) + 7.0/6.0*(dtheta1_c[2]*dv2[0]-dtheta1_c[0]*dv2[2]) + 
             0.5*(dtheta2_c[2]*dv2[0]-dtheta2_c[0]*dv2[2]) - 1.0/6.0*(dtheta2_c[2]*dv1[0]-dtheta2_c[0]*dv1[2]) + dv1[1] + dv2[1];
    dvt[2] = 0.5*(dtheta1_c[0]*dv1[1]-dtheta1_c[1]*dv1[0]) + 7.0/6.0*(dtheta1_c[0]*dv2[1]-dtheta1_c[1]*dv2[0]) + 
             0.5*(dtheta2_c[0]*dv2[1]-dtheta2_c[1]*dv2[0]) - 1.0/6.0*(dtheta2_c[0]*dv1[1]-dtheta2_c[1]*dv1[0]) + dv1[2] + dv2[2];

    winn[0] += wien[0]; // winn[0] = 2.0*wien[0] + wenn[0];
    winn[1] += wien[1]; // winn[1] = 2.0*wien[1] + wenn[1];
    winn[2] += wien[2]; // winn[2] = 2.0*wien[2] + wenn[2];

    dv[0] = Cnb[0][0]*dvt[0] + Cnb[1][0]*dvt[1] + Cnb[2][0]*dvt[2] - dt*(winn[1]*v0[2]-winn[2]*v0[1]);
    dv[1] = Cnb[0][1]*dvt[0] + Cnb[1][1]*dvt[1] + Cnb[2][1]*dvt[2] - dt*(winn[2]*v0[0]-winn[0]*v0[2]);
    dv[2] = Cnb[0][2]*dvt[0] + Cnb[1][2]*dvt[1] + Cnb[2][2]*dvt[2] - dt*(winn[0]*v0[1]-winn[1]*v0[0]) + dt*9.8;

    if((dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2]) > 1e-20) //can be commented
    {
        qvp[4] += dv[0];
        qvp[5] += dv[1];
        qvp[6] += dv[2];
    }

    /*----------------------------Position update----------------------------*/
    p0[0] = qvp[7];
    p0[1] = qvp[8];
    p0[2] = qvp[9];

    dp[0] = (v0[0]+qvp[4]) * dt_2;
    dp[1] = (v0[1]+qvp[5]) * dt_2;
    dp[2] = (v0[2]+qvp[6]) * dt_2;

    if((dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2]) > 1e-20) //can be commented
    {
        qvp[7] += dp[0] / Rm;
        qvp[8] += dp[1] / (Rn*cos_lat);
        qvp[9] -= dp[2];
    }

    /*----------------------------Attitude update----------------------------*/
    q[0] = qvp[0];
    q[1] = qvp[1];
    q[2] = qvp[2];
    q[3] = qvp[3];
    
    Phi[0] = dtheta1[0] + dtheta2[0] + 2.0/3.0*(dtheta1[1]*dtheta2[2] - dtheta1[2]*dtheta2[1]);
    Phi[1] = dtheta1[1] + dtheta2[1] + 2.0/3.0*(dtheta1[2]*dtheta2[0] - dtheta1[0]*dtheta2[2]);
    Phi[2] = dtheta1[2] + dtheta2[2] + 2.0/3.0*(dtheta1[0]*dtheta2[1] - dtheta1[1]*dtheta2[0]);

    phi = sqrt(Phi[0]*Phi[0] + Phi[1]*Phi[1] + Phi[2]*Phi[2]);
    if(phi > 1e-10)
    {
        qc[0] = sin(phi/2.0)/phi;
        qc[1] = Phi[0]*qc[0];
        qc[2] = Phi[1]*qc[0];
        qc[3] = Phi[2]*qc[0];
        qc[0] = cos(phi/2.0);
        quatmul(q, qc);
    }
    quat2dcm(Cn0b, q);
    dcmecef2ned(Ce0n0, p0[0], p0[1]);
    Cee0[0][0] =  cos(EARTH_W*dt);
    Cee0[0][1] = -sin(EARTH_W*dt);
    Cee0[0][2] = 0.0;
    Cee0[1][0] = -Cee0[0][1];
    Cee0[1][1] =  Cee0[0][0];
    Cee0[1][2] = 0.0;
    Cee0[2][0] = 0.0;
    Cee0[2][1] = 0.0;
    Cee0[2][2] = 1.0;
    dcmecef2ned(Cen, qvp[7], qvp[8]);
    mat3tran(Cne, Cen);
    mat3mul(Cne0, Cee0, Cne);
    mat3mul(Cnn0, Ce0n0, Cne0);
    mat3mul(Cnb, Cn0b, Cnn0);
    dcm2quat(q, Cnb);
    qm = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
    qvp[0] = q[0] / qm;
    qvp[1] = q[1] / qm;
    qvp[2] = q[2] / qm;
    qvp[3] = q[3] / qm;

    flag = 1;
    return 1;
}
