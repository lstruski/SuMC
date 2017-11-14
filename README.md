# SuMC

mkdir build
cd build

/home/lukasz/clion/bin/cmake/bin/cmake -DCMAKE_BUILD_TYPE=Debug -G "CodeBlocks - Unix Makefiles" ..
cd ..
/home/lukasz/clion/bin/cmake/bin/cmake --build ./build/ --target SuMC -- -j 2 


/home/lukasz/clion/bin/cmake/bin/cmake --build ./build/ --target clean -- -j 2
