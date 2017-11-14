# SuMC

## Requirements:
* cmake (3.8 or higher),
* gcc, g++ (5.4.0 or higher),
* lapack and blas libraries.

## Build

Go to directory [SuMC/SuMC/](https://github.com/struski2/SuMC/tree/master/SuMC) and run the following instruction in terminal:

```
mkdir build
cd build

cmake -DCMAKE_BUILD_TYPE=Debug -G "CodeBlocks - Unix Makefiles" ..

cd ..

cmake --build ./build/ --target SuMC -- -j 2
```
You run 
```
/home/lukasz/clion/bin/cmake/bin/cmake --build ./build/ --target clean -- -j 2
```
if you want to clean it.

## Demo

In folder [SuMC/demo/](https://github.com/struski2/SuMC/tree/master/demo) I putted sample data together [results](https://github.com/struski2/SuMC/tree/master/demo/res) of SuMC method. The following image presents dataset, which are consisted 5 (with three 2-dimensional and two 1-dimensional subspaces) subspaces.

<p align="center">
<img src="https://github.com/struski2/SuMC/blob/master/demo/orig_data.gif" width="500" height="500" />
</p>

#### How to run SuMC on this data?
```
./SuMC -i SuMC/demo/data.txt -k 5 -t 10 -c 0.66 -o SuMC/demo/res/
```

