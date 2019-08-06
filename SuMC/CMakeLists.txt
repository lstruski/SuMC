cmake_minimum_required(VERSION 3.14)
project(SuMC)

set(CMAKE_CXX_STANDARD 14)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -O2 -std=c++14")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -O3 -std=c++14 -pthread")

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

#set(CMAKE_BUILD_TYPE Release)
set(SOURCE_FILES ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} main.cpp cluster.cpp containerClusters.cpp lapackFunctions.cpp subspaceClustering.h)
#set(SOURCE_FILES ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} main_test.cpp cluster.cpp containerClusters.cpp lapackFunctions.cpp subspaceClustering.h)
add_executable(SuMC ${SOURCE_FILES})
#target_link_libraries(SuMC ${BLAS_LIBRARIES}; ${LAPACK_LIBRARIES}; ${LAPACKE_LIB})
target_link_libraries(SuMC ${BLAS_LIBRARIES}; ${LAPACK_LIBRARIES}; ${LAPACKE_LIB}; pthread)
