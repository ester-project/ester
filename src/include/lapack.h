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

void sgetrf_( int* m, int* n, float* a, int* lda, int* ipiv, int* info );
void sgetrs_( char* trans, int* n, int* nrhs, float* a, int* lda, int* ipiv, float* b, int* ldb, int* info );
void sgetri_( int* n, float* a, int* lda, int* ipiv, float* work, int* lwork, int* info );
void ssterf_( int* n, float* d, float* e, int* info );
void sgeequ_( int* m, int* n, float* a, int* lda, float* r, float* c, float* rowcnd, float* colcnd, float* amax, int* info );
void sgecon_( char* norm, int* n, float* a, int* lda, float* anorm, float* rcond, float* work, int* iwork, int* info );
float slange_( char* norm, int* m, int* n, float* a, int* lda, float* work );
void sgels_( char* trans, int* m, int* n, int* nrhs, float* a, int* lda, float* b, int* ldb, float* work, int* lwork, int* info );

}
#endif
