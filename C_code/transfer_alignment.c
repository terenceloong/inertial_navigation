#include "transfer_alignment.h"
#include "inertial_solve.h"
#include "axes.h"
#include "matrix_operation.h"
#include <string.h>
#include <math.h>

void initialize_ta_filter(ST_TA_FILTER *filter, double *P, double *Q, double *R)
{
    int i;

    memset(filter->bias, 0, sizeof(filter->bias));
    memset(filter->assembly, 0, sizeof(filter->assembly));
    memset(filter->X, 0, sizeof(filter->X));
    memset(filter->P, 0, sizeof(filter->P));
    for(i=0; i<sizeof(filter->X)/8; i++)
        filter->P[i][i] = P[i];
    memcpy(filter->Q, Q, sizeof(filter->Q));
    memcpy(filter->R, R, sizeof(filter->R));
}

// measure = [lat, lon, h, vn, ve, vd, psi, theta, gamma], deg
uint8_t transfer_alignment(double qvp[10], double acc[3], double *measure, ST_TA_FILTER *filter, double dt)
{
    double sin_lat, cos_lat, Rm, Rn, Rd;
    double Cnb[3][3];
    double wien[3], wenn[3], winn[3], fn[3];

    double Phi[N_STATE][N_STATE], Z[N_MEASURE], H[N_MEASURE][N_STATE];
    double P1[N_STATE][N_STATE], P2[N_STATE][N_STATE];
    double PH1[N_STATE][N_MEASURE], PH2[N_MEASURE][N_MEASURE], PH3[N_MEASURE][N_MEASURE], K[N_STATE][N_MEASURE], dX[N_STATE];

    double phi, qc[4];
    double Cab[3][3], att[3], Cna0[3][3], Cna[3][3], Cbn[3][3], C[3][3], Cab1[3][3];

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

    memset(Phi, 0, sizeof(Phi));
    Phi[0][1] =  winn[2] * dt;
    Phi[0][2] = -winn[1] * dt;
    Phi[1][0] = -winn[2] * dt;
    Phi[1][2] =  winn[0] * dt;
    Phi[2][0] =  winn[1] * dt;
    Phi[2][1] = -winn[0] * dt;

    fn[0] = Cnb[0][0]*acc[0] + Cnb[1][0]*acc[1] + Cnb[2][0]*acc[2];
    fn[1] = Cnb[0][1]*acc[0] + Cnb[1][1]*acc[1] + Cnb[2][1]*acc[2];
    fn[2] = Cnb[0][2]*acc[0] + Cnb[1][2]*acc[1] + Cnb[2][2]*acc[2];

    Phi[3][1] = -fn[2] * dt;
    Phi[3][2] =  fn[1] * dt;
    Phi[4][0] =  fn[2] * dt;
    Phi[4][2] = -fn[0] * dt;
    Phi[5][0] = -fn[1] * dt;
    Phi[5][1] =  fn[0] * dt;

    winn[0] = 2.0*wien[0] + wenn[0];
    winn[1] = 2.0*wien[1] + wenn[1];
    winn[2] = 2.0*wien[2] + wenn[2];

    Phi[3][4] =  winn[2] * dt;
    Phi[3][5] = -winn[1] * dt;
    Phi[4][3] = -winn[2] * dt;
    Phi[4][5] =  winn[0] * dt;
    Phi[5][3] =  winn[1] * dt;
    Phi[5][4] = -winn[0] * dt;

    Phi[6][3] =  dt / Rm;
    Phi[7][4] =  dt / (Rn*cos_lat);
    Phi[8][5] = -dt;

#if ALIGNMENT_MODLE == 1
    Z[0] = qvp[4] - measure[3]; //dvn
    Z[1] = qvp[5] - measure[4]; //dve
    Z[2] = qvp[6] - measure[5]; //dvd
    Z[3] = qvp[7] - measure[0] * 0.017453292519943; //dlat, rad
    Z[4] = qvp[8] - measure[1] * 0.017453292519943; //dlon, rad
    Z[5] = qvp[9] - measure[2]; //dh
    memset(H, 0, sizeof(H));
    H[0][3] = 1.0;
    H[1][4] = 1.0;
    H[2][5] = 1.0;
    H[3][6] = 1.0;
    H[4][7] = 1.0;
    H[5][8] = 1.0;
#elif ALIGNMENT_MODLE == 3
    Phi[0][12] = -Cnb[0][0] * dt;
    Phi[0][13] = -Cnb[1][0] * dt;
    Phi[0][14] = -Cnb[2][0] * dt;
    Phi[1][12] = -Cnb[0][1] * dt;
    Phi[1][13] = -Cnb[1][1] * dt;
    Phi[1][14] = -Cnb[2][1] * dt;
    Phi[2][12] = -Cnb[0][2] * dt;
    Phi[2][13] = -Cnb[1][2] * dt;
    Phi[2][14] = -Cnb[2][2] * dt;

    att[0] = measure[6] * 0.017453292519943;
    att[1] = measure[7] * 0.017453292519943;
    att[2] = measure[8] * 0.017453292519943;
    angle2dcm(Cna0, att);
    angle2dcm(Cab, filter->assembly);
    mat3mul(Cna, Cab, Cna0);
    mat3tran(Cbn, Cnb);
    mat3mul(C, Cbn, Cna);

    Z[0] = qvp[4] - measure[3]; //dvn
    Z[1] = qvp[5] - measure[4]; //dve
    Z[2] = qvp[6] - measure[5]; //dvd
    Z[3] = qvp[7] - measure[0] * 0.017453292519943; //dlat, rad
    Z[4] = qvp[8] - measure[1] * 0.017453292519943; //dlon, rad
    Z[5] = qvp[9] - measure[2]; //dh
    Z[6] = C[1][2];
    Z[7] = C[2][0];
    Z[8] = C[0][1];

    memset(H, 0, sizeof(H));
    H[0][3] = 1.0;
    H[1][4] = 1.0;
    H[2][5] = 1.0;
    H[3][6] = 1.0;
    H[4][7] = 1.0;
    H[5][8] = 1.0;
    H[6][0] = 1.0;
    H[7][1] = 1.0;
    H[8][2] = 1.0;
    H[6][9]  = Cbn[1][2]*Cna[1][2] - Cbn[1][1]*Cna[2][2];
    H[6][10] = Cbn[1][0]*Cna[2][2] - Cbn[1][2]*Cna[0][2];
    H[6][11] = Cbn[1][1]*Cna[0][2] - Cbn[1][0]*Cna[1][2];
    H[7][9]  = Cbn[2][2]*Cna[1][0] - Cbn[2][1]*Cna[2][0];
    H[7][10] = Cbn[2][0]*Cna[2][0] - Cbn[2][2]*Cna[0][0];
    H[7][11] = Cbn[2][1]*Cna[0][0] - Cbn[2][0]*Cna[1][0];
    H[8][9]  = Cbn[0][2]*Cna[1][1] - Cbn[0][1]*Cna[2][1];
    H[8][10] = Cbn[0][0]*Cna[2][1] - Cbn[0][2]*Cna[0][1];
    H[8][11] = Cbn[0][1]*Cna[0][1] - Cbn[0][0]*Cna[1][1];
#elif ALIGNMENT_MODLE == 5
    Phi[0][12] = -Cnb[0][0] * dt;
    Phi[0][13] = -Cnb[1][0] * dt;
    Phi[0][14] = -Cnb[2][0] * dt;
    Phi[1][12] = -Cnb[0][1] * dt;
    Phi[1][13] = -Cnb[1][1] * dt;
    Phi[1][14] = -Cnb[2][1] * dt;
    Phi[2][12] = -Cnb[0][2] * dt;
    Phi[2][13] = -Cnb[1][2] * dt;
    Phi[2][14] = -Cnb[2][2] * dt;

    Phi[3][15] = Cnb[0][0] * dt;
    Phi[3][16] = Cnb[1][0] * dt;
    Phi[3][17] = Cnb[2][0] * dt;
    Phi[4][15] = Cnb[0][1] * dt;
    Phi[4][16] = Cnb[1][1] * dt;
    Phi[4][17] = Cnb[2][1] * dt;
    Phi[5][15] = Cnb[0][2] * dt;
    Phi[5][16] = Cnb[1][2] * dt;
    Phi[5][17] = Cnb[2][2] * dt;

    att[0] = measure[6] * 0.017453292519943;
    att[1] = measure[7] * 0.017453292519943;
    att[2] = measure[8] * 0.017453292519943;
    angle2dcm(Cna0, att);
    angle2dcm(Cab, filter->assembly);
    mat3mul(Cna, Cab, Cna0);
    mat3tran(Cbn, Cnb);
    mat3mul(C, Cbn, Cna);

    Z[0] = qvp[4] - measure[3]; //dvn
    Z[1] = qvp[5] - measure[4]; //dve
    Z[2] = qvp[6] - measure[5]; //dvd
    Z[3] = qvp[7] - measure[0] * 0.017453292519943; //dlat, rad
    Z[4] = qvp[8] - measure[1] * 0.017453292519943; //dlon, rad
    Z[5] = qvp[9] - measure[2]; //dh
    Z[6] = C[1][2];
    Z[7] = C[2][0];
    Z[8] = C[0][1];

    memset(H, 0, sizeof(H));
    H[0][3] = 1.0;
    H[1][4] = 1.0;
    H[2][5] = 1.0;
    H[3][6] = 1.0;
    H[4][7] = 1.0;
    H[5][8] = 1.0;
    H[6][0] = 1.0;
    H[7][1] = 1.0;
    H[8][2] = 1.0;
    H[6][9]  = Cbn[1][2]*Cna[1][2] - Cbn[1][1]*Cna[2][2];
    H[6][10] = Cbn[1][0]*Cna[2][2] - Cbn[1][2]*Cna[0][2];
    H[6][11] = Cbn[1][1]*Cna[0][2] - Cbn[1][0]*Cna[1][2];
    H[7][9]  = Cbn[2][2]*Cna[1][0] - Cbn[2][1]*Cna[2][0];
    H[7][10] = Cbn[2][0]*Cna[2][0] - Cbn[2][2]*Cna[0][0];
    H[7][11] = Cbn[2][1]*Cna[0][0] - Cbn[2][0]*Cna[1][0];
    H[8][9]  = Cbn[0][2]*Cna[1][1] - Cbn[0][1]*Cna[2][1];
    H[8][10] = Cbn[0][0]*Cna[2][1] - Cbn[0][2]*Cna[0][1];
    H[8][11] = Cbn[0][1]*Cna[0][1] - Cbn[0][0]*Cna[1][1];
#endif

    matrixadddiagone((double*)Phi, N_STATE);
    matrixmul((double*)P1, (double*)Phi, (double*)filter->P, N_STATE, N_STATE, N_STATE); //P1 = Phi*P
    matrixmultran((double*)P2, (double*)P1, (double*)Phi, N_STATE, N_STATE, N_STATE); //P2 = P1*Phi' = Phi*P*Phi'
    matrixadddiag((double*)P2, filter->Q, N_STATE); //P2 = Phi*P*Phi'+Q

    matrixmultran((double*)PH1, (double*)P2, (double*)H, N_STATE, N_STATE, N_MEASURE); //PH1 = P2*H' = P*H'
    matrixmul((double*)PH2, (double*)H, (double*)PH1, N_MEASURE, N_STATE, N_MEASURE); //PH2 = H*PH1 = H*P*H'
    matrixadddiag((double*)PH2, filter->R, N_MEASURE); //PH2 = H*P*H'+R
    if(!matrixinv((double*)PH3, (double*)PH2, N_MEASURE)) //PH3 = PH2^-1 = (H*P*H'+R)^-1
        return 0;
    matrixmul((double*)K, (double*)PH1, (double*)PH3, N_STATE, N_MEASURE, N_MEASURE); //K = PH1*PH3 = P*H'/(H*P*H'+R)
    matrixmul((double*)P1, (double*)K, (double*)H, N_STATE, N_MEASURE, N_STATE); //P1 = K*H
    matrixmulconst((double*)P1, (double*)P1, -1.0, N_STATE, N_STATE); //P1 = -K*H
    matrixadddiagone((double*)P1, N_STATE); //P1 = I-K*H
    matrixmul((double*)filter->P, (double*)P1, (double*)P2, N_STATE, N_STATE, N_STATE); //P = (I-K*H)*P
    matrixsym((double*)filter->P, N_STATE);
    matrixmul(dX, (double*)K, Z, N_STATE, N_MEASURE, 1); //dX = K*Z

    /*----------------------------Adjust----------------------------*/
    phi = sqrt(dX[0]*dX[0] + dX[1]*dX[1] + dX[2]*dX[2]);
    if(phi > 0)
    {
        qc[0] = sin(phi/2.0)/phi;
        qc[1] = dX[0]*qc[0];
        qc[2] = dX[1]*qc[0];
        qc[3] = dX[2]*qc[0];
        qc[0] = cos(phi/2.0);
        quatmul(qc, qvp);
        qvp[0] = qc[0];
        qvp[1] = qc[1];
        qvp[2] = qc[2];
        qvp[3] = qc[3];
    }
    qvp[4] -= dX[3];
    qvp[5] -= dX[4];
    qvp[6] -= dX[5];
    qvp[7] -= dX[6];
    qvp[8] -= dX[7];
    qvp[9] -= dX[8];

#if ALIGNMENT_MODLE == 3
    phi = sqrt(dX[9]*dX[9] + dX[10]*dX[10] + dX[11]*dX[11]);
    if(phi > 0)
    {
        qc[0] = sin(phi/2.0)/phi;
        qc[1] = dX[9] *qc[0];
        qc[2] = dX[10]*qc[0];
        qc[3] = dX[11]*qc[0];
        qc[0] = cos(phi/2.0);
        quat2dcm(C, qc);
        mat3mul(Cab1, C, Cab);
        dcm2angle(filter->assembly, Cab1);
    }
    filter->bias[0] += dX[12];
    filter->bias[1] += dX[13];
    filter->bias[2] += dX[14];
#elif ALIGNMENT_MODLE == 5
    phi = sqrt(dX[9]*dX[9] + dX[10]*dX[10] + dX[11]*dX[11]);
    if(phi > 0)
    {
        qc[0] = sin(phi/2.0)/phi;
        qc[1] = dX[9] *qc[0];
        qc[2] = dX[10]*qc[0];
        qc[3] = dX[11]*qc[0];
        qc[0] = cos(phi/2.0);
        quat2dcm(C, qc);
        mat3mul(Cab1, C, Cab);
        dcm2angle(filter->assembly, Cab1);
    }
    filter->bias[0] += dX[12];
    filter->bias[1] += dX[13];
    filter->bias[2] += dX[14];
    filter->bias[3] += dX[15];
    filter->bias[4] += dX[16];
    filter->bias[5] += dX[17];
#endif

    return 1;
}
