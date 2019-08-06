//
// Created by lukasz struski on 12.11.17.
//

#ifndef SOURCE_SUBSPACECLUSTERING_H
#define SOURCE_SUBSPACECLUSTERING_H

#include <vector>
#include <lapacke.h>
#include <cblas.h>


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
void copy_vec(int dim, const double *x, double *y);

void data_param(std::string filename, int &n_samples, int &dim, char separator = ',');

void readfile(double *data, std::string filename, char separator = ',');

double *readData(std::string filename, int &n_samples, int &dim, char separator = ',');

/**
 *
 */
namespace subspaceClustering {
    const int bound = 2;

    class Cluster;

    /**
     * @class ContainerClusters
     *
     */
    class ContainerClusters {
        int size;
        double error;
        std::vector<Cluster *> tableClusters;

    public:
        ContainerClusters();

        ContainerClusters(const ContainerClusters &orig);

        ContainerClusters &operator=(const ContainerClusters &orig);

        virtual ~ContainerClusters();

        int getSize() const;

        double getError() const;

        void resetError();

        std::vector<Cluster *> getContainer() const;

        void clean();

        void stepHartigan(int *grups, std::vector<int> &activeClusters, const int size[], const double *data,
                          double allMemory, int bits = 0);

        void Hartigan(int *grups, const int size[], const double *data, int howClusters,
                      double allMemory, int iteration = 5, int bits = 0);

        void Hartigan_parallel(int *group, const int size[], const double *data, unsigned int howClusters,
                               double allMemory, unsigned int n_threads, int iteration = 5, int bits = 0);

        void compression(double *data, const int size[], const int *grups);

        void run(double *data, int *grups, const int size[], int howClusters,
                 double degreeOfCompression, int iteration = 5, int bits = 0);


//    private:
        void createCluster(const int size[], const double *data, const int *group, int howClusters,
                           double allMemory, int bits);

        void errorsTWOclusters(double *array, int idCluster1, int idCluster2, int totalWeight,
                               int bits = 0) const;

        int
        updateDim(int idCluster1, int idCluster2, int totalWeight, int bits = 0);
    };

    /**
     * @class Cluster
     *
     */
    class Cluster {
        double dim;
        double *mean;
        double *cov; // row-major order
        int weight;
        double memory;
        double *eigenvalues;
        double *eigenvectors;

        friend class ContainerClusters;

        friend std::ostream &operator<<(std::ostream &out, const Cluster &c);

    public:
        static int N;

        explicit Cluster(int N);

        Cluster(const Cluster &orig);

        Cluster &operator=(const Cluster &orig);

        virtual ~Cluster();

        double getDim() const;

        double *getMean() const;

        double *getCov() const;

        double *getEigenvalues() const;

        int getWeight() const;

        double getMemory() const;

        void compression(double *point) const;

//    private:  // todo: odznacz
        void changePoint(const double *point, int weightPoint, double *temp);

        double err(double factor) const;

        double errorFUN(int totalWeight, int bits = 0) const;

        bool unassign(int totalWeight, double toleranceFactor = 0.02) const;
    };

    std::ostream &operator<<(std::ostream &out, const Cluster &c);
} // namespace subspaceClustering

#endif //SOURCE_SUBSPACECLUSTERING_H
