#include <acml.h>

#define CblasRowMajor 101
#define CblasColMajor 102

#define CblasNoTrans    'N'
#define CblasTrans      'T'
#define CblasConjTrans  'C'

//extern void cblas_sgemm(char layout, char transa, char transb, int m, int n, int k, float alpha, float *a, int lda, float *b, int ldb, float beta, float *c, int ldc);

void cblas_sgemm(char layout, char transa, char transb, int m, int n, int k, float alpha, float *a, int lda, float *b, int ldb, float beta, float *c, int ldc);