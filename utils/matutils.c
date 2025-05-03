#include "matutils.h"

#define END 1
#define FREE_ARG char*

void nrerror(char error_text[])
{
    fprintf(stderr,"Run-time error...\n");
    fprintf(stderr,"%s\n",error_text);
    fprintf(stderr,"...now exiting to system...\n");
    exit(1);
}

float *vector(long nl, long nh) {
    float *v;
    v = (float *)malloc((size_t) ((nh - nl+1+END)*sizeof(float)));

    if (!v) nrerror("allocation failure in vector()");
    return v-nl+END;
}

int *ivector(long nl, long nh) {
    float *v;
    v = (float *)malloc((size_t) ((nh - nl+1+END)*sizeof(int)));

    if (!v) nrerror("allocation failure in ivector()");
    return v-nl+END;
}

char *cvector(long nl, long nh) {
    float *v;
    v = (float *)malloc((size_t) ((nh - nl+1+END)*sizeof(char)));

    if (!v) nrerror("allocation failure in cvector()");
    return v-nl+END;
}

unsigned long *lvector(long nl, long nh)
{
    unsigned long *v;
    v=(unsigned long *)malloc((size_t) ((nh-nl+1+END)*sizeof(long)));

    if (!v) nrerror("allocation failure in lvector()");
    return v-nl+END;
}

double *dvector(long nl, long nh)
{
    double *v;
    v=(double *)malloc((size_t) ((nh-nl+1+END)*sizeof(double)));

    if (!v) nrerror("allocation failure in dvector()");
    return v-nl+END;
}

float **matrix(long nrl, long nrh, long ncl, long nch)
{
    long i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
    float **m;

    m = (float **) malloc((size_t)((nrow+END)*sizeof(float*)));
    if (!m) nrerror("allocation failure 1 in matrixA()");

    m += END;
    m -= nrl;

    m[nrl]=(float *) malloc((size_t)((nrow*ncol+END)*sizeof(float)));
    if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
    m[nrl] += END;
    m[nrl] -= ncl;

    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

    return m;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
{
    long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
    double **m;

    m=(double **) malloc((size_t)((nrow+END)*sizeof(double*)));
    if (!m) nrerror("allocation failure 1 in matrix()");
    m += END;
    m -= nrl;

    m[nrl]=(double *) malloc((size_t)((nrow*ncol+END)*sizeof(double)));
    if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
    m[nrl] += END;
    m[nrl] -= ncl;
    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

    return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
{
    long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
    int **m;

    m=(int **) malloc((size_t)((nrow+END)*sizeof(int*)));
    if (!m) nrerror("allocation failure 1 in matrix()");
    m += END;
    m -= nrl;

    m[nrl]=(int *) malloc((size_t)((nrow*ncol+END)*sizeof(int)));
    if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
    m[nrl] += END;
    m[nrl] -= ncl;
    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

    return m;
}

float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch, long newrl, long newcl)
{
    long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
    float **m;

    m=(float **) malloc((size_t) ((nrow+END)*sizeof(float*)));
    if (!m) nrerror("allocation failure in submatrix()");
    m += END;
    m -= newrl;

    for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

    return m;
}

float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch)
{
    long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
    float **m;

    m=(float **) malloc((size_t) ((nrow+END)*sizeof(float*)));
    if (!m) nrerror("allocation failure in convert_matrix()");
    m += END;
    m -= nrl;

    m[nrl]=a-ncl;
    for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;

    return m;
}

void free_vector(float *v, long nl, long nh)
{
    free((FREE_ARG) (v+nl-END));
}

void free_ivector(int *v, long nl, long nh)
{
    free((FREE_ARG) (v+nl-END));
}

void free_cvector(unsigned char *v, long nl, long nh)
{
    free((FREE_ARG) (v+nl-END));
}

void free_lvector(unsigned long *v, long nl, long nh)
{
    free((FREE_ARG) (v+nl-END));
}

void free_dvector(double *v, long nl, long nh)
{
    free((FREE_ARG) (v+nl-END));
}

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
{
    free((FREE_ARG) (m[nrl]+ncl-END));
    free((FREE_ARG) (m+nrl-END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
    free((FREE_ARG) (m[nrl]+ncl-END));
    free((FREE_ARG) (m+nrl-END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)

{
    free((FREE_ARG) (m[nrl]+ncl-END));
    free((FREE_ARG) (m+nrl-END));
}
void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch)
{
    free((FREE_ARG) (b+nrl-END));
}
void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch)
{
    free((FREE_ARG) (b+nrl-END));
}