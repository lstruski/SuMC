//
// Created by lukasz on 12.11.17.
//

#ifndef SUMC_SUBSPACECLUSTERING_H_H
#define SUMC_SUBSPACECLUSTERING_H_H


#include <vector>
#include <lapacke.h>
#include <cblas.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * DCOPY copies a vector 'x' to a vector 'y'.
 */
void cblas_dcopy(const int N, const double *X, const int incX, double *Y, const int incY);


/*
 * DSYRK  performs one of the symmetric rank k operations
 *     C := alpha*A*A**T + beta*C,
 * or
 *     C := alpha*A**T*A + beta*C,
 * where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
 * and  A  is an  n by k  matrix in the first case and a  k by n  matrix
 * in the second case.
 */
void cblas_dsyrk(CBLAS_LAYOUT layout, CBLAS_UPLO Uplo,
                 CBLAS_TRANSPOSE Trans, const int N, const int K,
                 const double alpha, const double *A, const int lda,
                 const double beta, double *C, const int ldc);

/*
 * computes all eigenvalues and, optionally,
 * eigenvectors of a real symmetric matrix A.
 */
lapack_int LAPACKE_dsyevr(int matrix_order, char jobz, char range, char uplo,
                          lapack_int n, double *a, lapack_int lda, double vl, double vu,
                          lapack_int il, lapack_int iu, double abstol, lapack_int *m,
                          double *w, double *z, lapack_int ldz, lapack_int *isuppz);

lapack_int LAPACKE_dsyevx(int matrix_layout, char jobz, char range, char uplo,
                          lapack_int n, double *a, lapack_int lda, double vl,
                          double vu, lapack_int il, lapack_int iu,
                          double abstol, lapack_int *m, double *w, double *z,
                          lapack_int ldz, lapack_int *ifail);

#ifdef __cplusplus
}
#endif


/*
 * Auxiliary routines prototypes
 */
void updateCOV(int n, double *C, double *point, double alpha, double beta);

//int eigensystem(int n, double* A, double* eigenvalues, double* eigenvectors,
//        bool calculateEigenvectors = false);
int eigensystem(int dim, int n_select, double *A, double *eigenvalues,
                double *eigenvectors, bool calculateEigenvectors = false);

/*
 * This function copies a vector 'x' to a vector 'y'.
 */
void copy_vec(const int dim, const double *x, double *y);

double *readData(std::string filename, size_t &n_samples, size_t &dim, const char separator = ' ');

/**
 *
 */
namespace subspaceClustering {

    class Cluster;

    /**
     * @class ContainerClusters
     *
     */
    class ContainerClusters {
        size_t size;
        double error;
        std::vector<Cluster *> tableClusters;

    public:
        ContainerClusters();

        ContainerClusters(const ContainerClusters &orig);

        virtual ~ContainerClusters();

        size_t getSize() const;

        double getError() const;

        void resetError();

        std::vector<Cluster *> getContainer() const;

        void clean();

        ContainerClusters &operator=(const ContainerClusters &orig);

        void stepHartigan(size_t *grups, std::vector<size_t> &activeClusters, const size_t size[], const double *data,
                          const double allMemory, const size_t bits = 0);

        void Hartigan(size_t *grups, const size_t size[], const double *data, const size_t howClusters,
                      const double allMemory, const size_t iteration = 5, const size_t bits = 0);

        void compression(double *data, const size_t size[], const size_t *grups);

        void run(double *data, size_t *grups, const size_t size[], const size_t howClusters,
                 const double degreeOfCompression, const size_t iteration = 5, const size_t bits = 0);


    private:
        void createCluster(const size_t size[], const double *data, const size_t *grups, const size_t howClusters,
                           const double allMemory, const size_t bits);

        void errorsTWOclusters(double *array, const size_t idCluster1, const size_t idCluster2, const int totalWeight,
                               const size_t bits = 0) const;

        size_t
        updateDim(const size_t idCluster1, const size_t idCluster2, const int totalWeight, const size_t bits = 0);
    };

    /**
     * @class Cluster
     *
     */
    class Cluster {
        double dim;
        double *mean;
        double *cov; // row-major order
        size_t weight;
        double memory;
        double *eigenvalues;
        double *eigenvectors;

        friend class ContainerClusters;

        friend std::ostream &operator<<(std::ostream &out, const Cluster &c);

    public:
        static size_t N;

        Cluster(size_t N);

        Cluster(const Cluster &orig);

        virtual ~Cluster();

        double getDim() const;

        double *getMean() const;

        double *getCov() const;

        int getWeight() const;

        double getMemory() const;

        void compression(double *point) const;

    private:
        void changePoint(const double *point, const int weightPoint, double *temp);

        Cluster &operator=(const Cluster &orig);

        double err(const double factor) const;

        double errorFUN(const int totalWeight, const size_t bits = 0) const;

        bool unassign(const int totalWeight, const double toleranceFactor = 0.02) const;
    };

    std::ostream &operator<<(std::ostream &out, const Cluster &c);
} // namespace subspaceClustering


#endif //SUMC_SUBSPACECLUSTERING_H_H
