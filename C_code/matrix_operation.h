#ifndef __MATRIX_OPERATION_H
#define __MATRIX_OPERATION_H

#ifdef __cplusplus
extern "C" {
#endif

void mat3mul(double s[3][3], double a[3][3], double b[3][3]);
void mat3tran(double s[3][3], double a[3][3]);
int mat3inv(double s[3][3], double a[3][3]);
double mat3det(double a[3][3]);

void matrixadd(double *s, double *a, double *b, int m, int n);
void matrixsub(double *s, double *a, double *b, int m, int n);
void matrixmul(double *s, double *a, double *b, int l, int m, int n);
int matrixinv(double *s, double *a, int n); //will destroy a.
void matrixtran(double *s, double* a, int m, int n);
double matrixdet(double *a, int n); //will destroy a.

void matrixadddiag(double *a, double *b, int n);
void matrixadddiagone(double *a, int n);
void matrixmulconst(double *s, double *a, double b, int m, int n);
void matrixmuldiag(double *s, double *a, double *b, int m, int n);
void matrixmultran(double *s, double *a, double *b, int l, int m, int n);
void matrixsym(double *a, int n);

void matrixcpy(double *dst, double *src, int m, int n);
void matrixzeros(double *a, int m, int n); //set a to zero matrix.
void matrixeye(double *a, int n); //set a to identity matrix.
void matrixprint(double *a, int m, int n);

#ifdef __cplusplus
}
#endif

#endif
