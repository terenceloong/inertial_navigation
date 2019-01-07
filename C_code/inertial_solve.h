#ifndef _INERTIAL_SOLVE_H
#define _INERTIAL_SOLVE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

#define EARTH_A      6378137
#define EARTH_F      0.003352810664747
#define EARTH_W      7.292115e-5

void imu_compensate(double imu[6], double bias[6]);
uint8_t inertial_slove_rate(double qvp[10], double imu[6], double dt);
uint8_t inertial_slove_delta(double qvp[10], double imu[6], double dt);

#ifdef __cplusplus
}
#endif

#endif
