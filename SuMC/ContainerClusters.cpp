//
// Created by lukasz struski on 24.06.16.
//

#include <iostream>
#include <random>
#include "subspaceClustering.h"


using namespace subspaceClustering;

/*
 * This function calculates index for upper matrix without diagonal (indices (i, j) used as in classical 2-dimensional matrix).
 */
inline size_t index(size_t i, size_t j, size_t dim) {
    size_t res;
    if (i == j)
        res = -1;
    else if (i < j)
        res = (dim - 1) * i - (i * (i - 1)) / 2 + j - i - 1;
    else
        res = (dim - 1) * j - (j * (j - 1)) / 2 + i - j - 1;
    return res;
}

/**
 * Constructor
 */
ContainerClusters::ContainerClusters() {
    size = 0;
    error = std::numeric_limits<double>::infinity();
}

/**
 * Copy constructor
 */
ContainerClusters::ContainerClusters(const ContainerClusters &orig) {
    this->size = orig.size;
    this->error = orig.error;
    this->tableClusters.reserve(this->size);
    for (size_t i = 0; i < this->size; i++)
        this->tableClusters.push_back(new Cluster(*orig.tableClusters[i]));
}

/**
 * Destructor
 */
ContainerClusters::~ContainerClusters() {
    std::vector<Cluster *>::iterator ite;
    for (ite = tableClusters.begin(); ite != tableClusters.end(); ite++)
        delete *ite;
    tableClusters.clear();
}

/**
 * @return size of Container
 */
size_t ContainerClusters::getSize() const {
    return this->size;
}

/**
 * @return error of Container
 */
double ContainerClusters::getError() const {
    return this->error;
}

/**
 * reset error clustring
 */
void ContainerClusters::resetError() {
    this->error = std::numeric_limits<double>::infinity();
}

/**
 * @return the container of clusters
 */
std::vector<Cluster *> ContainerClusters::getContainer() const {
    return this->tableClusters;
}

/**
 * Cleaning the contents of the container
 */
void ContainerClusters::clean() {
    this->size = 0;
    this->error = std::numeric_limits<double>::infinity();
    for (std::vector<Cluster *>::iterator ite = this->tableClusters.begin(); ite != this->tableClusters.end(); ite++)
        delete *ite;
    this->tableClusters.clear();
}

/**
 * The (overload operator =) function copies the values of the 'orig' cluster container, deleting the previous data
 * 
 * @param orig (<i><b>ContainerClusters&</b></i>) - the container of clusters
 * @return copy 'orig'
 */
ContainerClusters &ContainerClusters::operator=(const ContainerClusters &orig) {
    std::vector<Cluster *>::reverse_iterator ite;
    for (ite = this->tableClusters.rbegin(); ite != this->tableClusters.rend(); ++ite)
        delete *ite;
    this->tableClusters.clear();

    this->size = orig.size;
    this->error = orig.error;
    this->tableClusters.reserve(this->size);
    for (size_t i = 0; i < this->size; i++)
        this->tableClusters.push_back(new Cluster(*orig.tableClusters[i]));
    return *this;
}

/**
 * Create clusters, randomly linking points to clusters
 * 
 * @param size (<i><b>size_t[]</b></i>) - size of number of samples and number of features
 * @param data (<i><b>double*</b></i>) - 2-dimensional array, which save data <i>row-major order</i>
 * @param groups (<i><b>size_t*</b></i>) - list of id clusters (values: 0,1,...,howClusters-1;  <i>row-major order</i>)
 * @param howClusters (<i><b>size_t</b></i>) - determines how much we want to create clusters
 * @param allMemory (<i><b>double</b></i>) - general memory we have at our disposal
 */
void
ContainerClusters::createCluster(const size_t size[], const double *data, const size_t *grups, const size_t howClusters,
                                 const double allMemory, const size_t bits) {
    this->size = howClusters;
    this->tableClusters.reserve(howClusters);
    size_t i, j;

    for (i = 0; i < howClusters; i++)
        this->tableClusters.push_back(new Cluster(size[1]));

    double *temp = new double[size[1]];

    for (i = 0; i < size[0]; i++) {
        j = grups[i];
        this->tableClusters[j]->changePoint(data + i * Cluster::N, 1, temp);
    }

    delete[] temp;
    double *matrix = new double[Cluster::N * Cluster::N];
    int info;
    double *eigenvec;
    int n_select;

    for (i = 0; i < howClusters; i++) {
        tableClusters[i]->memory = allMemory * (double) tableClusters[i]->weight / size[0];

        if (bits == 0)
            this->tableClusters[i]->dim = this->tableClusters[i]->memory /
                                          this->tableClusters[i]->weight;
        else
            this->tableClusters[i]->dim = (this->tableClusters[i]->memory +
                                           this->tableClusters[i]->weight *
                                           std::log2((double) this->tableClusters[i]->weight / size[0])) /
                                          (bits * this->tableClusters[i]->weight);

        for (j = 0; j < Cluster::N * Cluster::N; j++)
            matrix[j] = this->tableClusters[i]->cov[j];

        n_select = Cluster::N - (int) std::floor(this->tableClusters[i]->dim - 2);
        n_select = (n_select >= Cluster::N) ? Cluster::N : n_select;
        if (n_select > 0)
            try {
                info = eigensystem(Cluster::N, n_select, matrix, this->tableClusters[i]->eigenvalues, eigenvec);
                if (info != 0) throw "Calculating eigenvalues failed!";
            } catch (const char *str) {
                std::cerr << "Exception: " << str << "\n";
            }
    }
    delete[] matrix;
}

/**
 * The function calculates the minimum energy of two clusters and their memory
 * 
 * @param array (<i><b>double[2]</i></b>) -  array, which saves cost function and memory of second cluster
 * @param idCluster1 (<i><b>Cluster &</i></b>) - index of the first cluster
 * @param idCluster2 (<i><b>Cluster &</i></b>) - index of the second cluster
 * @param totalWeight (<i><b>int</i></b>) - the total number of data
 * @param bits (<i><b>size_t</i></b>) - number of bits needed to memorize one scalar
 */
void ContainerClusters::errorsTWOclusters(double *array, const size_t idCluster1, const size_t idCluster2,
                                          const int totalWeight, const size_t bits) const {
    double ilewsp;

    if (bits == 0)
        ilewsp = this->tableClusters[idCluster1]->memory + this->tableClusters[idCluster2]->memory;
    else
        ilewsp = (this->tableClusters[idCluster1]->memory + this->tableClusters[idCluster2]->memory +
                  std::log2((double) this->tableClusters[idCluster1]->weight / totalWeight) *
                  this->tableClusters[idCluster1]->weight +
                  std::log2((double) this->tableClusters[idCluster2]->weight / totalWeight) *
                  this->tableClusters[idCluster2]->weight) / bits;

    //error two clusters
    array[0] = std::numeric_limits<double>::infinity();
    //memory of cluster 'idCluster2'
    array[1] = this->tableClusters[idCluster2]->memory;

    if (ilewsp > std::numeric_limits<double>::epsilon()) {
        double dim1, dim2;
        double *temp = new double[4];

        // old dim1, dim2
        dim1 = this->tableClusters[idCluster1]->memory / this->tableClusters[idCluster1]->weight;
        dim2 = this->tableClusters[idCluster2]->memory / this->tableClusters[idCluster2]->weight;

        dim1 = std::round(dim1);
        dim2 = std::round(dim2);
        for (int i = -2; i <= 2; ++i) {
            temp[0] = (dim1 + i) * this->tableClusters[idCluster1]->weight;
            temp[2] = (dim2 + i) * this->tableClusters[idCluster2]->weight;

            if (temp[0] <= ilewsp) {
                temp[1] = ilewsp - temp[0];
                temp[3] = tableClusters[idCluster1]->err(dim1 + i) * tableClusters[idCluster1]->weight +
                          tableClusters[idCluster2]->err(temp[1] / tableClusters[idCluster2]->weight) *
                          tableClusters[idCluster2]->weight;
                if (temp[3] < array[0]) {
                    array[1] = temp[1];
                    array[0] = temp[3];
                }
            }

            if (temp[2] <= ilewsp) {
                temp[1] = ilewsp - temp[2];
                temp[3] = this->tableClusters[idCluster1]->err(temp[1] / this->tableClusters[idCluster1]->weight) *
                          this->tableClusters[idCluster1]->weight +
                          this->tableClusters[idCluster2]->err(dim2 + i) * this->tableClusters[idCluster2]->weight;

                if (temp[3] < array[0]) {
                    array[1] = temp[2];
                    array[0] = temp[3];
                }
            }
        }

        // new dim1, dim2
        dim1 = (ilewsp - array[1]) / this->tableClusters[idCluster1]->weight;
        dim2 = array[1] / this->tableClusters[idCluster2]->weight;

        this->tableClusters[idCluster1]->dim = dim1;
        this->tableClusters[idCluster2]->dim = dim2;

        delete[] temp;

        if (bits != 0)
            array[1] = -std::log2((double) this->tableClusters[idCluster2]->weight / totalWeight) *
                       this->tableClusters[idCluster2]->weight + bits * array[1];
        // new memory of cluster 'idCluster1' = old memories Cluster1+Cluster2 - new memory Cluster2
    }
}

/**
 * The function improves the size and memory of two clusters
 * 
 * @param idCluster1 (<i><b>Cluster &</i></b>) - index of the first cluster
 * @param idCluster2 (<i><b>Cluster &</i></b>) - index of the second cluster
 * @param totalWeight (<i><b>int</i></b>) - the total number of data
 * @param bits (<i><b>size_t</i></b>) - number of bits needed to memorize one scalar
 * 
 * @return returns 0 or 1, 0 when cluster 'idCluster1' has an integer dimension, 1 for the second case
 */
size_t ContainerClusters::updateDim(const size_t idCluster1, const size_t idCluster2, const int totalWeight,
                                    const size_t bits) {
    size_t ret = 0;
    double ilewsp;
    if (bits == 0)
        ilewsp = this->tableClusters[idCluster1]->memory + this->tableClusters[idCluster2]->memory;
    else
        ilewsp = (this->tableClusters[idCluster1]->memory + this->tableClusters[idCluster2]->memory +
                  std::log2((double) this->tableClusters[idCluster1]->weight / totalWeight) *
                  this->tableClusters[idCluster1]->weight +
                  std::log2((double) this->tableClusters[idCluster2]->weight / totalWeight) *
                  this->tableClusters[idCluster2]->weight) / bits;

    double tempM, blad = std::numeric_limits<double>::infinity();
    if (bits == 0) {
        this->tableClusters[idCluster1]->dim =
                this->tableClusters[idCluster1]->memory / this->tableClusters[idCluster1]->weight;
        this->tableClusters[idCluster2]->dim =
                this->tableClusters[idCluster2]->memory / this->tableClusters[idCluster2]->weight;
    } else {
        this->tableClusters[idCluster1]->dim = (this->tableClusters[idCluster1]->memory +
                                                this->tableClusters[idCluster1]->weight * std::log2(
                                                        (double) this->tableClusters[idCluster1]->weight /
                                                        totalWeight)) /
                                               (bits * this->tableClusters[idCluster1]->weight);
        this->tableClusters[idCluster2]->dim = (this->tableClusters[idCluster2]->memory +
                                                this->tableClusters[idCluster2]->weight * std::log2(
                                                        (double) this->tableClusters[idCluster2]->weight /
                                                        totalWeight)) /
                                               (bits * this->tableClusters[idCluster2]->weight);
    }

    if (this->tableClusters[idCluster1]->dim - std::floor(this->tableClusters[idCluster1]->dim) >
        std::numeric_limits<double>::epsilon())
        ret = 1;

    if (ilewsp > std::numeric_limits<double>::epsilon()) {
        double dim1, dim2;

        double *temp = new double[4];

        dim1 = std::round(this->tableClusters[idCluster1]->dim);
        dim2 = std::round(this->tableClusters[idCluster2]->dim);
        for (int i = -2; i <= 2; ++i) {
            temp[0] = (dim1 + i) * this->tableClusters[idCluster1]->weight;
            temp[2] = (dim2 + i) * this->tableClusters[idCluster2]->weight;

            if (temp[0] <= ilewsp) {
                temp[1] = ilewsp - temp[0];
                temp[3] = this->tableClusters[idCluster1]->err(dim1 + i) * this->tableClusters[idCluster1]->weight +
                          this->tableClusters[idCluster2]->err(temp[1] / this->tableClusters[idCluster2]->weight) *
                          this->tableClusters[idCluster2]->weight;
                if (temp[3] < blad) {
                    tempM = temp[1];
                    blad = temp[3];
                    this->tableClusters[idCluster1]->dim = dim1 + i;
                    this->tableClusters[idCluster2]->dim = temp[1] / this->tableClusters[idCluster2]->weight;
                    ret = 0;
                }
            }

            if (temp[2] <= ilewsp) {
                temp[1] = ilewsp - temp[2];
                temp[3] = this->tableClusters[idCluster1]->err(temp[1] / this->tableClusters[idCluster1]->weight) *
                          this->tableClusters[idCluster1]->weight +
                          this->tableClusters[idCluster2]->err(dim2 + i) * this->tableClusters[idCluster2]->weight;

                if (temp[3] < blad) {
                    tempM = temp[2];
                    blad = temp[3];
                    this->tableClusters[idCluster1]->dim = temp[1] / this->tableClusters[idCluster1]->weight;
                    this->tableClusters[idCluster2]->dim = dim2 + i;
                    ret = 1;
                }
            }
        }

        delete[] temp;

        this->tableClusters[idCluster1]->memory =
                this->tableClusters[idCluster1]->memory + this->tableClusters[idCluster2]->memory;
        if (bits == 0)
            this->tableClusters[idCluster2]->memory = tempM;
        else
            this->tableClusters[idCluster2]->memory =
                    -std::log2((double) this->tableClusters[idCluster2]->weight / totalWeight) *
                    this->tableClusters[idCluster2]->weight + bits * tempM;
        this->tableClusters[idCluster1]->memory -= this->tableClusters[idCluster2]->memory;
    }
    return ret;
}

/**
 * Hartigan iteration
 * 
 * @param size (<i><b>size_t[]</b></i>) - size of data
 * @param data (<i><b>double*</b></i>) - 2-dimensional array, which save data <i>row-major order</i>
 * @param groups (<i><b>size_t*</b></i>) - list of id clusters (values: 0,1,...,howClusters-1;  <i>row-major order</i>)
 * @param howClusters (<i><b>size_t</b></i>) - determines how much we want to create clusters
 * @param allMemory (<i><b>double</b></i>) -  general memory we have at our disposal
 * @param bits (<i><b>size_t</i></b>) - number of bits needed to memorize one scalar
 */
void ContainerClusters::stepHartigan(size_t *grups, std::vector<size_t> &activeClusters, const size_t size[],
                                     const double *data, const double allMemory, const size_t bits) {
    bool outside, T = false;
    size_t howClusters = activeClusters.back() + 1;
    this->createCluster(size, data, grups, howClusters, allMemory, bits);

    std::vector<size_t>::iterator it, min_it;

    size_t l, j, m, i, min_weight, id_i, n_select;
    double temp, c, tempMemory;
    double *point = new double[size[1]];
    double *tempVec = new double[size[1]];
    double array[2];
    double *matrix = new double[Cluster::N * Cluster::N];
    int info;
    double *eigenvec;

    ContainerClusters *copy = new ContainerClusters();

    size_t stop = size[0];
    bool check = true;

    size_t items = 0;

    double error2cluster;
    double *errorsTWOclustersArray = new double[(howClusters * (howClusters - 1)) / 2];
    std::fill(errorsTWOclustersArray, errorsTWOclustersArray + (howClusters * (howClusters - 1)) / 2, -1.0);

    while (!T && items < 50) {
        items++;
        T = true;
        for (m = 0; m < size[0] && !(m == stop && check); m++) {
            check = true;
            copy_vec(size[1], data + m * size[1], point);

            l = grups[m];

            if (bits != 0) {
                outside = true;
                for (it = activeClusters.begin(); it != activeClusters.end(); it++)
                    if (l == *it) {
                        outside = false;
                        break;
                    }

                if (outside) {
                    temp = this->tableClusters[l]->memory / this->tableClusters[l]->weight;
                    l = activeClusters.front();
                    grups[m] = l;
                    this->tableClusters[l]->changePoint(point, 1, tempVec);
                    this->tableClusters[l]->memory += temp;

                    for (i = 0; i < Cluster::N * Cluster::N; i++)
                        matrix[i] = this->tableClusters[l]->cov[i];

                    if (bits == 0)
                        this->tableClusters[l]->dim = this->tableClusters[l]->memory /
                                                      this->tableClusters[l]->weight;
                    else
                        this->tableClusters[l]->dim = (this->tableClusters[l]->memory +
                                                       this->tableClusters[l]->weight *
                                                       std::log2((double) this->tableClusters[l]->weight / size[0])) /
                                                      (bits * this->tableClusters[l]->weight);
                    n_select = Cluster::N - (int) std::floor(this->tableClusters[l]->dim - 2);
                    n_select = (n_select >= Cluster::N) ? Cluster::N : n_select;
                    if (n_select > 0)
                        try {
                            info = eigensystem(Cluster::N, n_select, matrix, this->tableClusters[l]->eigenvalues,
                                               eigenvec);
                            if (info != 0) throw "Calculating eigenvalues failed!";
                        } catch (const char *str) {
                            std::cerr << "Exception: " << str << "\n";
                        }

                    T = false;
                    stop = m + 1;
                    check = false;
                }
            }

            *copy = *this;
            copy->tableClusters[l]->changePoint(point, -1, tempVec);

            for (i = 0; i < Cluster::N * Cluster::N; i++)
                matrix[i] = copy->tableClusters[l]->cov[i];

            n_select = Cluster::N - (int) std::floor(copy->tableClusters[l]->dim - 2);
            n_select = (n_select >= Cluster::N) ? Cluster::N : n_select;
            if (n_select > 0)
                try {
                    info = eigensystem(Cluster::N, n_select, matrix, copy->tableClusters[l]->eigenvalues, eigenvec);
                    if (info != 0) throw "Calculating eigenvalues failed!";
                } catch (const char *str) {
                    std::cerr << "Exception: " << str << "\n";
                }

            j = l;
            c = 0.0;

            for (it = activeClusters.begin(); it != activeClusters.end(); it++) {
                if (*it != l) {
                    id_i = index(l, *it, howClusters);
                    if (errorsTWOclustersArray[id_i] == -1) {
                        this->errorsTWOclusters(array, l, *it, size[0], bits);
                        temp = array[0];
                        errorsTWOclustersArray[id_i] = temp;
                    } else
                        temp = errorsTWOclustersArray[id_i];

                    copy->tableClusters[*it]->changePoint(point, 1, tempVec);

                    for (i = 0; i < Cluster::N * Cluster::N; i++)
                        matrix[i] = copy->tableClusters[*it]->cov[i];

                    n_select = Cluster::N - (int) std::floor(copy->tableClusters[*it]->dim - 2);
                    n_select = (n_select >= Cluster::N) ? Cluster::N : n_select;
                    if (n_select > 0)
                        try {
                            info = eigensystem(Cluster::N, n_select, matrix, copy->tableClusters[*it]->eigenvalues,
                                               eigenvec);
                            if (info != 0) throw "Calculating eigenvalues failed!";
                        } catch (const char *str) {
                            std::cerr << "Exception: " << str << "\n";
                        }

                    copy->errorsTWOclusters(array, l, *it, size[0], bits);

                    if (array[0] < (c + temp)) {
                        j = *it;
                        c = array[0] - temp;
                        tempMemory = array[1];
                        error2cluster = array[0];
                    }
                }
                if (c == -std::numeric_limits<double>::infinity()) break;
            }

            if (j != l) {
                T = false;
                stop = m + 1;
                check = false;
                grups[m] = j;

                for (int kk = 0; kk < howClusters; ++kk) {
                    if (l != kk) {
                        id_i = index(l, kk, howClusters);
                        errorsTWOclustersArray[id_i] = -1.0;
                    }
                    if (j != kk) {
                        id_i = index(j, kk, howClusters);
                        errorsTWOclustersArray[id_i] = -1.0;
                    }
                }
                id_i = index(l, j, howClusters);
                errorsTWOclustersArray[id_i] = error2cluster;

                this->tableClusters[l]->weight = copy->tableClusters[l]->weight;
                this->tableClusters[j]->weight = copy->tableClusters[j]->weight;
                for (size_t ii = 0; ii < size[1]; ii++) {
                    this->tableClusters[l]->mean[ii] = copy->tableClusters[l]->mean[ii];
                    this->tableClusters[j]->mean[ii] = copy->tableClusters[j]->mean[ii];
                    this->tableClusters[l]->eigenvalues[ii] = copy->tableClusters[l]->eigenvalues[ii];
                    this->tableClusters[j]->eigenvalues[ii] = copy->tableClusters[j]->eigenvalues[ii];
                    for (size_t jj = 0; jj < size[1]; jj++) {
                        this->tableClusters[l]->cov[ii * size[1] + jj] = copy->tableClusters[l]->cov[ii * size[1] + jj];
                        this->tableClusters[j]->cov[ii * size[1] + jj] = copy->tableClusters[j]->cov[ii * size[1] + jj];
                    }
                }

                this->tableClusters[l]->memory =
                        this->tableClusters[l]->memory + this->tableClusters[j]->memory - tempMemory;
                this->tableClusters[j]->memory = tempMemory;

                this->tableClusters[l]->dim = copy->tableClusters[l]->dim;
                this->tableClusters[j]->dim = copy->tableClusters[j]->dim;

            }
        }

        // Delete cluster which have a small number of points
        if (bits != 0) {
            min_weight = size[0];
            for (it = activeClusters.begin(); it != activeClusters.end(); it++)
                if (min_weight >= this->tableClusters[*it]->weight) {
                    min_it = it;
                    min_weight = this->tableClusters[*it]->weight;
                }
            if (this->tableClusters[*min_it]->unassign(size[0])) {
                activeClusters.erase(min_it);
                T = false;
            }
            if (activeClusters.size() == 0) {
                std::cout << "Delete all clusters. Small number of points in relation to the dimension of the data.\t";
                T = true;
            }
        }


    }
    delete copy;
    delete[] point;
    delete[] tempVec;
    delete[] matrix;
    delete[] errorsTWOclustersArray;

    //update dimensions and memory of clusters
    if (activeClusters.size() == 1) {
        this->tableClusters[activeClusters.front()]->dim =
                this->tableClusters[activeClusters.front()]->memory / size[0];
        if (bits != 0)
            this->tableClusters[activeClusters.front()]->dim /= bits;
    } else {
        std::vector<size_t> vec;
        vec.reserve(activeClusters.size());
        for (it = activeClusters.begin(); it != activeClusters.end(); ++it)
            vec.push_back(*it);
        while (vec.size() > 1) {
            it = vec.begin();
            i = this->updateDim(*it, *(it + 1), size[0], bits);
            vec.erase(it + i);
        }
        vec.clear();
        vec.resize(0);
    }

    this->error = 0.0;
    if (bits != 0) {
        for (int i = this->size - 1; i >= 0; i--) {
            outside = true;
            for (std::vector<size_t>::reverse_iterator rit = activeClusters.rbegin();
                 rit != activeClusters.rend(); ++rit)
                if (i == *rit) {
                    outside = false;
                    break;
                }
            if (outside) {
                delete this->tableClusters[i];
                this->tableClusters.erase(this->tableClusters.begin() + i);
                this->size -= 1;
            } else {
                this->error += this->tableClusters[i]->errorFUN(size[0], bits);
            }
        }
    } else {
        for (int i = this->size - 1; i >= 0; i--)
            this->error += this->tableClusters[i]->errorFUN(size[0], bits);
    }
}

/**
 *
 * @param groups (<i><b>size_t*</b></i>) - list of id clusters (values: 0,1,...,howClusters-1;  <i>row-major order</i>)
 * @param size (<i><b>size_t[]</b></i>) - size of data
 * @param data (<i><b>double*</b></i>) - 2-dimensional array, which save data <i>row-major order</i>
 * @param howClusters (<i><b>size_t</b></i>) - determines how much we want to create clusters
 * @param allMemory (<i><b>double</b></i>) -  general memory we have at our disposal
 * @param bits (<i><b>size_t</i></b>) - number of bits needed to memorize one scalar
 * @param iteration (<i><b>size_t</i></b>) - number of iterations
 *
 */
void ContainerClusters::Hartigan(size_t *grups, const size_t size[], const double *data, const size_t howClusters,
                                 const double allMemory, const size_t iteration, const size_t bits) {

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, howClusters - 1);

    size_t *TEMPgrups = new size_t[size[0]];
    ContainerClusters *temp = new ContainerClusters();

    std::vector<size_t> *TEMPactiveClusters = new std::vector<size_t>();
    TEMPactiveClusters->reserve(howClusters);

    std::vector<size_t> *activeClusters;
    if (bits != 0) {
        activeClusters = new std::vector<size_t>();
        activeClusters->reserve(howClusters);
    }

    size_t i, j;

    for (i = 0; i < iteration; i++) {
        for (j = 0; j < howClusters; j++)
            TEMPactiveClusters->push_back(j);
        for (j = 0; j < size[0]; j++)
            TEMPgrups[j] = dis(gen);

        temp->stepHartigan(TEMPgrups, *TEMPactiveClusters, size, data, allMemory, bits);

        if (temp->error < this->error) {
            *this = *temp;
            for (j = 0; j < size[0]; j++)
                grups[j] = TEMPgrups[j];
            if (bits != 0) activeClusters->assign(TEMPactiveClusters->begin(), TEMPactiveClusters->end());
        }

        TEMPactiveClusters->clear();
        temp->clean();
    }

    if (bits != 0 && this->size != howClusters) {
        j = 0;
        for (i = 0; i < size[0]; i++) {
            for (std::vector<size_t>::iterator it = activeClusters->begin(); it != activeClusters->end(); it++) {
                if (grups[i] == *it) {
                    grups[i] = j;
                    break;
                }
                j++;
            }
            j = 0;
        }
    }

    delete temp;
    delete[] TEMPgrups;
    TEMPactiveClusters->clear();
    delete TEMPactiveClusters;
    if (bits != 0) {
        activeClusters->clear();
        delete activeClusters;
    }
}

/**
 * 
 * function compresses data - changes "data" array
 *
 * @param size (<i><b>size_t[]</b></i>) - size of data
 * @param data (<i><b>double*</b></i>) - 2-dimensional array, which save data <i>row-major order</i>
 * @param groups (<i><b>size_t*</b></i>) - list of id clusters (values: 0,1,...,howClusters-1;  <i>row-major order</i>)
 *
 */
void ContainerClusters::compression(double *data, const size_t size[], const size_t *grups) {
    size_t i, j;

    // create array cluster eigenvectors
    double *matrix = new double[Cluster::N * Cluster::N];
    double *temp = new double[Cluster::N];
    int info, n_select;

    try {
        for (j = 0; j < this->size; j++) {
            n_select = Cluster::N - (int) std::floor(this->tableClusters[j]->dim);
            this->tableClusters[j]->eigenvectors = new double[Cluster::N * n_select];
            for (i = 0; i < Cluster::N * Cluster::N; i++)
                matrix[i] = this->tableClusters[j]->cov[i];
            // one eigenvector in one column (stored columnwise)
            //            info = eigensystem(Cluster::N, matrix, temp, this->tableClusters[j]->eigenvectors, true);
            info = eigensystem(Cluster::N, -n_select, matrix, temp, this->tableClusters[j]->eigenvectors, true);
            if (info != 0) throw "Calculating eigenvectors failed!";
        }
    } catch (const char *str) {
        std::cerr << "Exception: " << str << "\n";
    }

    delete[] matrix;
    delete[] temp;

    for (i = 0; i < size[0]; i++) {
        j = grups[i];
        this->tableClusters[j]->compression(data + i * size[1]);
    }

    // delete eigenvectors because I do not cleare memory in destructors
    for (i = 0; i < this->size; i++)
        delete[] this->tableClusters[i]->eigenvectors;
}

/**
 *
 * @param data (<i><b>double*</b></i>) - 2-dimensional array, which save data <i>row-major order</i>
 * @param groups (<i><b>size_t*</b></i>) - list of id clusters (values: 0,1,...,howClusters-1;  <i>row-major order</i>)
 * @param size (<i><b>size_t[]</b></i>) - size of data
 * @param howClusters (<i><b>size_t</b></i>) - determines how much we want to create clusters
 * @param iteration (<i><b>size_t</i></b>) - number of iterations
 * @param bits (<i><b>size_t</i></b>) - number of bits needed to memorize one scalar
 * @param allMemory (<i><b>double</b></i>) -  general memory we have at our disposal
 * @param degreeOfCompression - Compression ratio
 *
 */
void ContainerClusters::run(double *data, size_t *grups, const size_t size[], const size_t howClusters,
                            const double degreeOfCompression, const size_t iteration, const size_t bits) {
    double allMemory = (bits == 0) ? degreeOfCompression * size[0] * size[1] : degreeOfCompression * bits * size[0] *
                                                                               size[1];
    this->Hartigan(grups, size, data, howClusters, allMemory, iteration, bits);
    this->compression(data, size, grups);
}