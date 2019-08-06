//
// Created by lukasz struski on 24.06.16.
//


#include <iostream>
#include <fstream>
#include <algorithm>
#include "subspaceClustering.h"

using namespace subspaceClustering;

/*
 * This function copies a vector 'x' to a vector 'y'.
 */
void copy_vec(const int dim, const double *x, double *y) {
    cblas_dcopy(dim, x, 1, y, 1);
}

void data_param(std::string filename, int &n_samples, int &dim, const char separator) {
    std::ifstream in(filename, std::ios::in);
    std::string line;
    n_samples = 0;

    if (!in) {
        std::cout << "Cannot open input file.\n";
        exit(1);
    }

    std::getline(in, line);
    dim = (int) std::count(line.begin(), line.end(), separator);
    ++dim;
    ++n_samples;
    while (std::getline(in, line))
        ++n_samples;

    in.clear();
    in.close();

    if (!in.good()) {
        std::cout << "A file error occurred.\n";
        exit(1);
    }
}


void readfile(double *data, std::string filename, char separator) {
    std::ifstream in(filename, std::ios::in);
    std::string line, item;
    int i;

    if (!in) {
        std::cout << "Cannot open input file.\n";
        exit(1);
    }

    i = 0;
    std::stringstream linestream;
    while (std::getline(in, line)) {
        linestream << line;
        while (std::getline(linestream, item, separator))
            if (line.length() > 0)
                data[i++] = std::stod(item);
        linestream.clear();
    }

    in.clear();
    in.close();

    if (!in.good()) {
        std::cout << "A file error occurred.\n";
        exit(1);
    }
}


/*
 * This function reads data from file.
 */
double *readData(std::string filename, int &n_samples, int &dim, const char separator) {
    std::ifstream in(filename, std::ios::in);
    std::string line, item;
    double *data;
    int i;

    n_samples = 0;

    if (!in) {
        std::cout << "Cannot open input file.\n";
        exit(1);
    }

    std::getline(in, line);
    dim = (int) std::count(line.begin(), line.end(), separator);
    ++dim;
    ++n_samples;
    while (std::getline(in, line))
        ++n_samples;

    in.clear();
    in.seekg(0, std::ios::beg);

    i = 0;
    data = new double[n_samples * dim];
    std::stringstream linestream;
    while (std::getline(in, line)) {
        linestream << line;
        while (std::getline(linestream, item, separator))
            if (line.length() > 0)
                data[i++] = std::stod(item);
        linestream.clear();
    }

    in.clear();
    in.close();

    if (!in.good()) {
        std::cout << "A file error occurred.\n";
        exit(1);
    }
    return data;
}


int Cluster::N;

/**
 * @param N (<i><b>int</b></i>) - size of cluster (what dimension are the data represented)
 */
Cluster::Cluster(int N) {
    Cluster::N = N;
    this->dim = 0.0;
    this->mean = new double[this->N];
    this->eigenvalues = new double[this->N];
    this->cov = new double[this->N * this->N];
    for (int i = 0; i < this->N; i++) {
        this->mean[i] = 0;
        for (int j = 0; j < this->N; j++)
            this->cov[i * this->N + j] = 0;
    }
    this->weight = 0;
    this->memory = 0.0;
}

/**
 * @param orig (<i><b>Cluster</b></i>)
 */
Cluster::Cluster(const Cluster &orig) {
    this->N = orig.N;
    this->dim = orig.dim;
    this->mean = new double[this->N];
    this->eigenvalues = new double[this->N];
    this->cov = new double[this->N * this->N];
    for (int i = 0; i < orig.N; i++) {
        this->mean[i] = orig.mean[i];
        this->eigenvalues[i] = orig.eigenvalues[i];
        for (int j = 0; j < orig.N; j++)
            this->cov[i * orig.N + j] = orig.cov[i * orig.N + j];
    }
    this->weight = orig.weight;
    this->memory = orig.memory;
}

/**
 * Destructor
 */
Cluster::~Cluster() {
    delete[] mean;
    delete[] cov;
    delete[] eigenvalues;
}

/**
 * @return dim of Cluster
 */
double Cluster::getDim() const {
    return this->dim;
}

/**
 * @return mean of Cluster
 */
double *Cluster::getMean() const {
    return this->mean;
}

/**
 * @return cov of Cluster
 */
double *Cluster::getCov() const {
    return this->cov;
}

/**
 * @return eigenvalues of Cluster
 */
double *Cluster::getEigenvalues() const {
    return this->eigenvalues;
}

/**
 * @return weight of Cluster
 */
int Cluster::getWeight() const {
    return this->weight;
}

/**
 * @return memory of Cluster
 */
double Cluster::getMemory() const {
    return this->memory;
}

/**
 * The function changes all (except 'N', 'dim', memory ',' error ') data in a cluster caused by adding or removing a point to/from a cluster.
 *
 * @param point (<i><b>double*</i></b>) - point, which will be added to the cluster
 * @param weightPoint (<i><b>int</i></b>) - weight of this point: '1' when we add, '-1' when we subtract 'point' from the cluster
 * @param temp (<i><b>double*</i></b>) - temporary vector of point size
 */
void Cluster::changePoint(const double *point, const int weightPoint, double *temp) {
    double pu, pv;
    pu = (double) this->weight / (this->weight + weightPoint);
    pv = (double) weightPoint / (this->weight + weightPoint);

    for (int i = 0; i < this->N; i++) {
        temp[i] = point[i] - this->mean[i];
        this->mean[i] = pu * this->mean[i] + pv * point[i];
    }

    pv *= pu;
    updateCOV(this->N, this->cov, temp, pv, pu);

    this->weight += weightPoint;
}

/**
 * The function compresses the point from a given cluster
 *
 * @param point (<i><b>double*</i></b>) - the coordinates of the point to be compressed
 */
void Cluster::compression(double *point) const {
    auto *temp = new double[this->N];
    int i, j, n_select = this->N - (int) std::floor(this->dim);

    double s = 0.0;
    for (i = 0; i < this->N; i++) {
        temp[i] = point[i] - this->mean[i];
        s += temp[i] * this->eigenvectors[i * n_select + n_select - 1];
    }

    for (i = 0; i < this->N; i++)
        point[i] = (this->dim - (int) std::floor(this->dim)) * s *
                   this->eigenvectors[i * n_select + n_select - 1];

    if (this->dim >= 1) {
        for (i = 0; i <= (int) std::floor(this->dim) - 1; i++) {
            s = 0.0;
            for (j = 0; j < this->N; j++)
                s += temp[j] * this->eigenvectors[j * n_select + this->N - 1 - i];

            for (j = 0; j < this->N; j++)
                point[j] += s * this->eigenvectors[j * n_select + this->N - 1 - i];
        }
    }

    for (i = 0; i < this->N; i++)
        point[i] += this->mean[i];

    delete[] temp;
}

/**
 * The function (overloading operator =) copies the contents of the 'orig'
 *
 * @param orig (<i><b>Cluster&</b></i>) - cluster
 * @return copy of 'orig'
 */
Cluster &Cluster::operator=(const Cluster &orig) {
    try {
        if (this->N != orig.N) throw "Those clusters have different size!!!";
        this->dim = orig.dim;
        for (int i = 0; i < orig.N; i++) {
            this->mean[i] = orig.mean[i];
            this->eigenvalues[i] = orig.eigenvalues[i];
            for (int j = 0; j < N; j++)
                this->cov[i * orig.N + j] = orig.cov[i * orig.N + j];
        }
        this->weight = orig.weight;
        this->memory = orig.memory;
    } catch (const char *str) {
        std::cerr << "Exception: " << str << "\n";
    }
    return *this;
}

/**
 * @param factor (<i><b>double</i></b>) - ratio of memory to number of points (weight)
 * @return accuracy/error with number of parameters
 */
double Cluster::err(const double factor) const {
    double s = 0.0;
    if (factor >= (double) this->N) return s;
    else if (factor <= std::numeric_limits<double>::epsilon()) {
        for (int i = 0; i < this->N; i++)
            s += this->cov[i * this->N + i];
    } else if (factor > this->N / 2) {
        for (int i = 0; i < this->N - std::ceil(factor); i++)
            s += this->eigenvalues[i];
        s += (std::ceil(factor) - factor) * this->eigenvalues[(int) (this->N - 1 - std::floor(factor))];
    } else {
        for (int i = 0; i < this->N; i++)
            s += this->cov[i * this->N + i];
        for (int i = this->N - std::floor(factor); i < this->N; i++)
            s -= this->eigenvalues[i];
        s -= (factor - std::floor(factor)) * this->eigenvalues[(int) (this->N - 1 - std::floor(factor))];
//        for (int i = 1; i < static_cast<int>(floor(factor) + 1); i++)
//            s -= this->eigenvalues[i];
//        s -= (factor - std::floor(factor)) * this->eigenvalues[0];
    }
    return s;
}

/**
 * Cost function
 *
 * @param totalWeight (<i><b>int</i></b>) - the total number of data
 * @param bits (<i><b>int</i></b>) - number of bits needed to memorize one scalar
 * @return cost function for one cluster
 */
double Cluster::errorFUN(const int totalWeight, const int bits) const {
    double factor;
    if (bits == 0)
        factor = this->memory / this->weight;
    else
        factor = (this->memory + std::log2((double) this->weight / totalWeight) * this->weight) / (bits * this->weight);
    return this->err(factor) * this->weight;
}

/**
 * The function checks whether the cluster is deleted or not
 *
 * @param totalWeight (<i><b>int</i></b>) - the total number of data
 * @param toleranceFactor (<i><b>double</i></b>) - the coefficient defining the minimum number of points in relation to the total number of data sets
 * @return whether to remove the cluster (true) or not (false)
 */
bool Cluster::unassign(const int totalWeight, const double toleranceFactor) const {
    bool ret = false;
    if (this->weight <= std::max(toleranceFactor * totalWeight, 2.0 * this->N)) ret = true;
    return ret;
}

/**
 * Overloaded exit operator for 'Cluster' class
 *
 * @return Prints cluster data: dimension, memory, weight, center, covariance matrix
 */
std::ostream &subspaceClustering::operator<<(std::ostream &out, const Cluster &c) {
    int i, j;
    out << "\nDim:\t" << c.dim;
    out << "\nMemory:\t" << c.memory;
    out << "\nWeight:\t" << c.weight;
    out << "\nMean:\t(";
    for (i = 0; i < c.N - 1; i++)
        out << c.mean[i] << ",";
    out << c.mean[c.N - 1] << ")";

    out << "\nCovariance matrix:\n";
    for (i = 0; i < c.N; i++) {
        for (j = 0; j < c.N; j++)
            if (i <= j) {
                out << c.cov[i * c.N + j] << " ";
            } else {
                out << c.cov[j * c.N + i] << " ";
            }
        out << "\n";
    }

//    out << "\nEigenvalues: ";
//    for (i = 0; i < c.N; i++)
//        out << c.eigenvalues[i] << ", ";
//    out << "\n";
    return out;
}
