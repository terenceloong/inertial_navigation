#ifndef __AXES_H
#define __AXES_H

#ifdef __cplusplus
extern "C" {
#endif

void angle2dcm(double t[3][3], double a[3]);
void angle2quat(double q[4], double a[3]);
void quat2dcm(double t[3][3], double q[4]);
void quat2angle(double a[3], double q[4]);
void dcm2quat(double q[4], double t[3][3]);
void dcm2angle(double a[3], double t[3][3]);

void quatmul(double a[4], double b[4]);

void dcmecef2ned(double t[3][3], double lat, double lon);

#ifdef __cplusplus
}
#endif

#endif
