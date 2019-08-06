//
// Created by lukasz struski on 24.06.16.
//


#include <getopt.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include "subspaceClustering.h"


using namespace std;
using namespace subspaceClustering;

std::string inputfile, outputDir = "../results_SuMC/";
char delimiter = ' ';
int k = 10, iters = 5, bits = 0;
double comp_ratio;

void PrintHelp() {
    std::cout <<
              "--input -i:          Input file\n"
              "--delimiter -d:      The char used to separate values\n"
              "--k -k:              Number of clusters\n"
              "--comp_ratio -c:     Compression ratio\n"
              "--iters -t:          Number of iterations\n"
              "--bits -b:           Number of bits\n"
              "--output -o:         Path to output directory\n";
    exit(1);
}

void ProcessArgs(int argc, char **argv) {
    const char *const short_opts = "i:d:k:c:t:b:o:h";
    const option long_opts[] = {
            {"input",      1, nullptr, 'i'},
            {"delimiter",  1, nullptr, 'd'},
            {"k",          1, nullptr, 'k'},
            {"comp_ratio", 1, nullptr, 'c'},
            {"iters",      0, nullptr, 't'},
            {"bits",       0, nullptr, 'b'},
            {"output",     0, nullptr, 'o'},
            {"help",       0, nullptr, 'h'},
            {nullptr,      0, nullptr, 0}
    };

    while (true) {
        const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

        if (-1 == opt)
            break;

        switch (opt) {
            case 'i':
                inputfile = std::string(optarg);
                std::cout << "Input file: " << inputfile << std::endl;
                break;

            case 'd':
                delimiter = *optarg;
                std::cout << "Delimiter: " << delimiter << std::endl;
                break;

            case 'k':
                k = std::stoi(optarg);
                std::cout << "Number of clusters:" << k << std::endl;
                break;

            case 't':
                iters = std::stoi(optarg);
                std::cout << "Number of iterations: " << iters << std::endl;
                break;

            case 'c':
                comp_ratio = std::stof(optarg);
                std::cout << "Compression ratio: " << comp_ratio << std::endl;
                break;

            case 'b':
                bits = std::stoi(optarg);
                std::cout << "Number of bits:" << bits << std::endl;
                break;

            case 'o':
                outputDir = std::string(optarg);
                std::cout << "Path to output directory: " << outputDir << std::endl;
                break;

            case 'h': // -h or --help
            default:
                PrintHelp();
                break;
        }
    }
}


/*
 *
 */
int main(int argc, char **argv) {
    ProcessArgs(argc, argv);


//    // create directory
//    inputfile = "mkdir -p " + outputDir;
//    const int dir_err = system(inputfile.c_str());
//    if (-1 == dir_err) {
//        printf("Error creating directory!n");
//        exit(1);
//    }

    std::chrono::high_resolution_clock::time_point start, finish;
    std::chrono::duration<double> elapsed;

    ContainerClusters c;
    int size[2];

    double *data = readData(inputfile, *size, *(size + 1), delimiter); // n_samples = size[0], dim = size[1]
    auto *grups = new int[size[0]];

    double p_Memory;

    p_Memory = (bits == 0) ? comp_ratio * size[0] * size[1] : comp_ratio * bits * size[0] * size[1];

    start = std::chrono::high_resolution_clock::now();
    c.Hartigan(grups, size, data, k, p_Memory, iters, bits); // calculate time of working this function
    finish = std::chrono::high_resolution_clock::now();
    elapsed = finish - start;
    printf("Elapsed time: %.5G sec\n", elapsed.count());


    if (c.getSize() != 0) {

        ofstream plik;

        string output_file = outputDir + "/SuMC_labels.txt";
        plik.open(output_file.c_str(), ios::out);
        if (plik.good()) {
            for (int i = 0; i < size[0]; i++)
                plik << grups[i] + 1 << "\n";
            plik.close();
        } else cout << "Access to file \"" << output_file << "\" has been forbidden!" << endl;


        output_file = outputDir + "/SuMC_description.txt";
        plik.open(output_file.c_str(), ios::out);
        if (plik.good()) {
            plik << "Number of clusters: " << c.getSize() << "\n";
            plik << "Memory:\n";
            for (int i = 0; i < c.getSize(); i++)
                plik << c.getContainer()[i]->getMemory() << "\n";

            plik << "\nDimension:\n";
            for (int i = 0; i < c.getSize(); i++)
                plik << c.getContainer()[i]->getDim() << "\n";

            plik << "\nWeight:\n";
            for (int i = 0; i < c.getSize(); i++)
                plik << c.getContainer()[i]->getWeight() << "\n";

            plik << "\n============================================\nResults:\n";
            for (int i = 0; i < c.getSize(); i++)
                plik << *c.getContainer()[i];

            plik << "\n\n============================================\nE = ";
            plik << c.getError() << "\n";
            plik.close();
        } else cout << "Access to file \"" << output_file << "\" has been forbidden!" << endl;

        double bladSreKwad = 0.0;

        for (int i = 0; i < size[0]; i++)
            for (int j = 0; j < size[1]; j++)
                bladSreKwad += pow(abs(data[i * size[1] + j] - c.getContainer()[grups[i]]->getMean()[j]), 2);

        cout << " Cost function: " << c.getError() / size[0] << "\t" << c.getError() / bladSreKwad << "\n";
    }

    delete[] grups;
    delete[] data;
    return 0;
}
