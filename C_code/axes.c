#include <math.h>

#define PI_2     1.5707963267949

void angle2dcm(double t[3][3], double a[3])
{
	double cx, cy, cz, sx, sy, sz;

	cx = cos(a[2]);
	cy = cos(a[1]);
	cz = cos(a[0]);
	sx = sin(a[2]);
	sy = sin(a[1]);
	sz = sin(a[0]);
	t[0][0] = cy*cz;
	t[0][1] = cy*sz;
	t[0][2] = -sy;
	t[1][0] = sy*sx*cz - sz*cx;
	t[1][1] = sy*sx*sz + cz*cx;
	t[1][2] = cy*sx;
	t[2][0] = sy*cx*cz + sz*sx;
	t[2][1] = sy*cx*sz - cz*sx;
	t[2][2] = cy*cx;
}

void angle2quat(double q[4], double a[3])
{
	double cx, cy, cz, sx, sy, sz;

	cx = cos(a[2]/2);
	cy = cos(a[1]/2);
	cz = cos(a[0]/2);
	sx = sin(a[2]/2);
	sy = sin(a[1]/2);
	sz = sin(a[0]/2);
	q[0] = cz*cy*cx + sz*sy*sx;
	q[1] = cz*cy*sx - sz*sy*cx;
	q[2] = cz*sy*cx + sz*cy*sx;
	q[3] = sz*cy*cx - cz*sy*sx;
}

void quat2dcm(double t[3][3], double q[4])
{
	double mod_q, q0, q1, q2, q3;
	
	mod_q = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
	q0 = q[0] / mod_q;
	q1 = q[1] / mod_q;
	q2 = q[2] / mod_q;
	q3 = q[3] / mod_q;
	t[0][0] = q0*q0 + q1*q1 - q2*q2 - q3*q3;
	t[0][1] = 2 * (q1*q2 + q0*q3);
	t[0][2] = 2 * (q1*q3 - q0*q2);
	t[1][0] = 2 * (q1*q2 - q0*q3);
	t[1][1] = q0*q0 - q1*q1 + q2*q2 - q3*q3;
	t[1][2] = 2 * (q2*q3 + q0*q1);
	t[2][0] = 2 * (q1*q3 + q0*q2);
	t[2][1] = 2 * (q2*q3 - q0*q1);
	t[2][2] = q0*q0 - q1*q1 - q2*q2 + q3*q3;
}

void quat2angle(double a[3], double q[4])
{
	double mod_q, q0, q1, q2, q3;
	double t00, t01, t02, t12, t22;

	mod_q = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
	q0 = q[0] / mod_q;
	q1 = q[1] / mod_q;
	q2 = q[2] / mod_q;
	q3 = q[3] / mod_q;
	t02 = 2 * (q0*q2 - q1*q3);
	if (t02 >= 1.0)
	{
		t00 = q0*q0 - q1*q1 + q2*q2 - q3*q3;
		t01 = 2 * (q0*q3 - q1*q2);
		a[0] = atan2(t01, t00);
		a[1] = PI_2;
		a[2] = 0.0;
	}
	else if (t02 <= -1.0)
	{
		t00 = q0*q0 - q1*q1 + q2*q2 - q3*q3;
		t01 = 2 * (q0*q3 - q1*q2);
		a[0] = atan2(t01, t00);
		a[1] = -PI_2;
		a[2] = 0.0;
	}
	else
	{
		t00 = q0*q0 + q1*q1 - q2*q2 - q3*q3;
		t01 = 2 * (q1*q2 + q0*q3);
		t12 = 2 * (q2*q3 + q0*q1);
		t22 = q0*q0 - q1*q1 - q2*q2 + q3*q3;
		a[0] = atan2(t01, t00);
		a[1] = asin(t02);
		a[2] = atan2(t12, t22);
	}
}

void dcm2quat(double q[4], double t[3][3])
{
	if (t[0][0] >= t[1][1]+t[2][2])
	{
		q[1] = 0.5 * sqrt(1.0 + t[0][0] - t[1][1] - t[2][2]);
		q[0] = 0.25 * (t[1][2] - t[2][1]) / q[1];
		q[2] = 0.25 * (t[0][1] + t[1][0]) / q[1];
		q[3] = 0.25 * (t[2][0] + t[0][2]) / q[1];
	}
	else if (t[1][1] >= t[0][0]+t[2][2])
	{
		q[2] = 0.5 * sqrt(1.0 - t[0][0] + t[1][1] - t[2][2]);
		q[0] = 0.25 * (t[2][0] - t[0][2]) / q[2];
		q[1] = 0.25 * (t[0][1] + t[1][0]) / q[2];
		q[3] = 0.25 * (t[1][2] + t[2][1]) / q[2];
	}
	else if (t[2][2] >= t[0][0]+t[1][1])
	{
		q[3] = 0.5 * sqrt(1.0 - t[0][0] - t[1][1] + t[2][2]);
		q[0] = 0.25 * (t[0][1] - t[1][0]) / q[3];
		q[1] = 0.25 * (t[2][0] + t[0][2]) / q[3];
		q[2] = 0.25 * (t[1][2] + t[2][1]) / q[3];
	}
	else
	{
		q[0] = 0.5 * sqrt(1.0 + t[0][0] + t[1][1] + t[2][2]);
		q[1] = 0.25 * (t[1][2] - t[2][1]) / q[0];
		q[2] = 0.25 * (t[2][0] - t[0][2]) / q[0];
		q[3] = 0.25 * (t[0][1] - t[1][0]) / q[0];
	}
}

// 0.999999999 avoids the problem at calculating atan2, the two input arguments can't be zero.
void dcm2angle(double a[3], double t[3][3])
{
	if (t[0][2] >= 0.9999999999)
	{
		a[0] = atan2(-t[1][0], t[1][1]);
		a[1] = -PI_2;
		a[2] = 0.0;
	}
	else if (t[0][2] <= -0.9999999999)
	{
		a[0] = atan2(-t[1][0], t[1][1]);
		a[1] = PI_2;
		a[2] = 0.0;
	}
	else
	{
		a[0] = atan2(t[0][1], t[0][0]);
		a[1] = asin(-t[0][2]);
		a[2] = atan2(t[1][2], t[2][2]);
	}
}

void quatmul(double a[4], double b[4])
{
	double a0, a1, a2, a3;
	
	a0 = a[0];
	a1 = a[1];
	a2 = a[2];
	a3 = a[3];
	a[0] = b[0]*a0 - b[1]*a1 - b[2]*a2 - b[3]*a3;
	a[1] = b[0]*a1 + b[1]*a0 - b[2]*a3 + b[3]*a2;
	a[2] = b[0]*a2 + b[1]*a3 + b[2]*a0 - b[3]*a1;
	a[3] = b[0]*a3 - b[1]*a2 + b[2]*a1 + b[3]*a0;
}

void dcmecef2ned(double t[3][3], double lat, double lon)
{
	double sin_lat, cos_lat, sin_lon, cos_lon;

	sin_lat = sin(lat);
	cos_lat = cos(lat);
	sin_lon = sin(lon);
	cos_lon = cos(lon);

	t[0][0] = -sin_lat * cos_lon;
	t[0][1] = -sin_lat * sin_lon;
	t[0][2] =  cos_lat;
	t[1][0] = -sin_lon;
	t[1][1] =  cos_lon;
	t[1][2] =  0.0;
	t[2][0] = -cos_lat * cos_lon;
	t[2][1] = -cos_lat * sin_lon;
	t[2][2] = -sin_lat;
}
