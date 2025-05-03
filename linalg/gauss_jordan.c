#include "linalg.h"

void searchPiv(
    float **a, int n, 
    int icol, int irow, 
    int i, int j, 
    int k, int *ipiv, 
    float big) {

    big = 0.0;
    for (j=1; j<=n; j++)
        if (ipiv[j] != 1)
            for (k=1; k<=n; k++) {
                if (ipiv[k] == 0) {
                    if (fabs(a[j][k]) >= big) {
                        big = fabs(a[j][k]);
                        irow = j;
                        icol = k;
                    }
                }
            }
    ++(ipiv[icol]);
}

void pivotAndReduce(
    float **a, int icol, 
    float pivinv, int n, 
    int m, int l, 
    int ll, float dum, 
    float **b) {
    
    if (a[icol][icol] == 0.0) nrerror("GaussJ: Singular matrix encountered.");
    pivinv = 1.0/a[icol][icol];
    a[icol][icol] = 1.0;

    for (l=1;l<=n;l++) a[icol][l] *= pivinv;
    for (l=1;l<=m;l++) b[icol][l] *= pivinv;

    for (ll=1;ll<=n;ll++)
        if (ll != icol) {
            dum = a[ll][icol];
            a[ll][icol] = 0.0;

            for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
            for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
        }
}

void gaussJElim(
    float **a, 
    int n,
    int m,
    float **b) {

    int *indxc, *indxr, *ipiv;
    float big, dum, pivinv, temp;
    int i,icol,irow,j,k,l,ll;

    indxc = ivector(1, n);
    indxr = ivector(1, n);
    ipiv = ivector(1, n);
    for (j=1; j<=n; j++) ipiv[j]=0;
    
    for (i=1; i<=n; i++) {
        searchPiv(
            a, n, icol, irow, i, j, k, ipiv, big);

        if (irow != icol) {
            for (l=1;l<=n;l++) SWAP(a[irow][l], a[icol][l])
            for (l=1;l<=m;l++) SWAP(b[irow][l], b[icol][l])
        }   
        indxr[i] = irow;
        indxc[i] = icol;

        pivotAndReduce(
            a, icol, pivinv, n, m, l, ll, dum, b);
    }
    
    for (l=n;l>=1;l--) {
        if (indxr[l] != indxc[l])
            for (k=1;k<=n;k++)
                SWAP(a[k][indxr[l]], a[k][indxc[l]]);
    }

    free_ivector(ipiv,1,n);
    free_ivector(indxr,1,n);
    free_ivector(indxc,1,n);
}