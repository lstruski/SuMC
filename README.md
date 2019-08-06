# SuMC

Implementation of method from paper ''[*Lossy compression approach to subspace clustering*](https://www.sciencedirect.com/science/article/pii/S0020025516311628
)''.
## Requirements:
* cmake (3.8 or higher),
* gcc, g++ (5.4.0 or higher),
* lapack and blas libraries.

## Build

Go to directory [SuMC/SuMC/](https://github.com/struski2/SuMC/tree/master/SuMC) and run the following instruction in terminal:

```
mkdir build
cd build

cmake -DCMAKE_BUILD_TYPE=Release -G "CodeBlocks - Unix Makefiles" ..

cd ..

cmake --build ./build/ --target SuMC -- -j 2
```
You run 
```
cmake --build ./build/ --target clean -- -j 2
```
if you want to clean it.

## Demo

In folder [SuMC/demo/](https://github.com/struski2/SuMC/tree/master/demo) I putted sample data together [results](https://github.com/struski2/SuMC/tree/master/demo/res) of SuMC method. The following image presents dataset, which are consisted 5 (with three 2-dimensional and two 1-dimensional subspaces) subspaces.

<p align="center">
<img src="https://github.com/struski2/SuMC/blob/master/demo/orig_data.gif" width="500" height="500" />
</p>

#### How to run SuMC on this data?
```
./SuMC -i SuMC/demo/data.txt -k 5 -t 10 -c 0.5835 -o SuMC/demo/res/
```

Legend of command line interpreter:
--input -i:          Input file
--delimiter -d:      The char used to separate values
--k -k:              Number of clusters
--comp_ratio -c:     Compression ratio
--iters -t:          Number of iterations
--bits -b:           Number of bits
--output -o:         Path to output directory

