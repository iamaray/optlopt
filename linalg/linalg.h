#ifndef _LINALG_H_
#define _LINALG_H_

#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

#include <math.h>
#include "matutils.h"

void gaussJElim(float **a, int n, float **b, int m);
float **luDecomp(float **a, int n, int *indx, float *d);
float **luBackSub(float **a, int n, int *indx, float b[]);
float **inv(float **a);
float **pinv(float ** a);
float **svd(float **a, int m, int n, float w[], float **v);

#endif