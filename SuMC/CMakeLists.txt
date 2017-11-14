cmake_minimum_required(VERSION 3.8)
project(SuMC)

#set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -O2 -std=c++11")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -m64 -O2 -std=c++11")

find_package(BLAS REQUIRED)
if (BLAS_FOUND)
    #    message(WARNING "BLAS_LIB FOUND at " ${BLAS_LIBRARIES} " include files ")
    link_directories(${BLAS_LIBRARIES})
else (BLAS_FOUND)
    message(FATAL_ERROR "BLAS NOT FOUND")
endif (BLAS_FOUND)

find_package(LAPACK REQUIRED)
if (LAPACK_FOUND)
    #    message(WARNING "LAPACK_LIB FOUND at " ${LAPACK_LIBRARIES} ", include files ")
    link_directories(${LAPACK_LIBRARIES})
else (LAPACK_FOUND)
    message(FATAL_ERROR "LAPACK NOT FOUND")
endif (LAPACK_FOUND)

find_library(LAPACKE_LIB lapacke /usr/lib/)

set(SOURCE_FILES ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} main.cpp Cluster.cpp ContainerClusters.cpp lapackFunctions.cpp subspaceClustering.h)
add_executable(SuMC ${SOURCE_FILES})
target_link_libraries(SuMC ${BLAS_LIBRARIES}; ${LAPACK_LIBRARIES}; ${LAPACKE_LIB})
