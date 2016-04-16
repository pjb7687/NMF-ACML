#include "blas.h"

void cblas_sgemm(char layout, char transa, char transb, int m, int n, int k, float alpha, float *a, int lda, float *b, int ldb, float beta, float *c, int ldc) {
	if (layout == CblasColMajor)
		sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
	else if (layout == CblasRowMajor)
		sgemm(transb, transa, n, m, k, alpha, b, ldb, a, lda, beta, c, ldc);
}