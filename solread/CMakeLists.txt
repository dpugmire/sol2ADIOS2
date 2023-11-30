cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

project(solread C CXX)

find_package(MPI REQUIRED)
find_package(ADIOS2 REQUIRED)
find_package(HDF5 REQUIRED)

# We are not using the C++ API of MPI, this will stop the compiler look for it
add_definitions(-DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX)

add_executable(solread-adios solread-adios.cpp)
target_link_libraries(solread-adios adios2::cxx11_mpi MPI::MPI_C)


add_executable(solread-hdf5 solread-hdf5.cpp)
target_link_libraries(solread-hdf5 adios2::cxx11_mpi hdf5::hdf5 MPI::MPI_C)
