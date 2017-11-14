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

![data](SuMC/demo/orig_data.gif)
