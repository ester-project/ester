#ifndef _LAPACK_H_
#define _LAPACK_H_

extern "C" {
void dgetrf_( int* m, int* n, double* a, int* lda, int* ipiv, int* info );
void dgetrs_( char* trans, int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info );
void dgetri_( int* n, double* a, int* lda, int* ipiv, double* work, int* lwork, int* info );
void dsterf_( int* n, double* d, double* e, int* info );
void dgeequ_( int* m, int* n, double* a, int* lda, double* r, double* c, double* rowcnd, double* colcnd, double* amax, int* info );
void dgecon_( char* norm, int* n, double* a, int* lda, double* anorm, double* rcond, double* work, int* iwork, int* info );
double dlange_( char* norm, int* m, int* n, double* a, int* lda, double* work );
void dgels_( char* trans, int* m, int* n, int* nrhs, double* a, int* lda, double* b, int* ldb, double* work, int* lwork, int* info );
}
#endif
