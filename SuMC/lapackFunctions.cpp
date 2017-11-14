//
// Created by lukasz struski on 24.06.16.
//

#include "subspaceClustering.h"

/** 
 * update covariance matrix C := alpha*x**T*x + beta*C
 * 
 * @param n (<i><b>int</i></b>) - size of symmetric matrix
 * @param C (<i><b>double*</i></b>) - symmetric matrix
 * @param point (<i><b>double*</i></b>) - vector, which update matrix
 * @param alpha (<i><b>double</i></b>) - sclar
 * @param beta (<i><b>double</i></b>) - sclar
 */
void updateCOV(int n, double *C, double *point, double alpha, double beta) {
    cblas_dsyrk(CblasRowMajor, CblasUpper, CblasNoTrans, n, 1, alpha,
                point, 1, beta, C, n);
}


/**
 * calculate all eigenvalues and optional eigenvectors
 * 
 * @param dim (<i><b>int</i></b>) - size of symmetric matrix n x n
 * @param n_select (<i><b>int</i></b>) - if n_select>0 then calculate number of smallest eigenvalues; otherwise n_select largest eigenvalues to be returned
 * @param A (<i><b>double*</i></b>) - symmetric matrix
 * @param eigenvalues (<i><b>double*</i></b>) - eigenvalues of A (vector with size 'dim')
 * @param eigenvectors (<i><b>double*</i></b>) - eigrnvectors of A (matrix with size [n_select * dim])
 * @param calculateEigenvectors (<i><b>bool</i></b>) - whether calculate eigenvectors (=true)
 */
int
eigensystem(int dim, int n_select, double *A, double *eigenvalues, double *eigenvectors, bool calculateEigenvectors) {
    int n = dim, il, iu, m, lda = dim, ldz = n_select, info;
    double abstol, vl, vu;
    char jobz = 'N', range = 'I';
    /*
     jobz
            Must be 'N' or 'V'.
            If jobz = 'N', then only eigenvalues are computed.
            If jobz = 'V', then eigenvalues and eigenvectors are computed.
    range
            Must be 'A', 'V', or 'I'.
            If range = 'A', all eigenvalues will be found.
            If range = 'V', all eigenvalues in the half-open interval (vl, vu] will be found.
            If range = 'I', the eigenvalues with indices il through iu will be found.
     */

    /* Local arrays */
    int *isuppz;

    if (n_select == dim)
        range = 'A';
    else
        /* Set il, iu to compute NSELECT smallest eigenvalues */
    if (n_select > 0) {
        il = 1;
        iu = n_select;
    } else {
        ldz = -n_select;
        il = dim + n_select + 1;
        iu = dim;
    }
    if (calculateEigenvectors) {
        jobz = 'V';
        isuppz = new int[2 * n];
    }

    /* Negative abstol means using the default value */
    abstol = -1.0;

    /* Solve eigenproblem */
    info = LAPACKE_dsyevr(LAPACK_ROW_MAJOR, jobz, range, 'U', n, A, lda,
                          vl, vu, il, iu, abstol, &m, eigenvalues, eigenvectors, ldz, isuppz);

    if (calculateEigenvectors)
        delete[] isuppz;
    /* Check for convergence */
    if (info > 0) {
        printf("The algorithm failed to compute eigenvalues.\n");
        exit(1);
    } else
        return info;
}