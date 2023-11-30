# sol2ADIOS2

Converts GE's datafiles from CGNS to ADIOS2


Build on whoopingcough.

cmake \
 -DCMAKE_BUILD_TYPE=Debug \
 -DADIOS2_DIR=/apps/ADIOS/install_par/lib/cmake/adios2 \
 -DHDF5_DIR=/apps/hdf5/install \
../solread

Build on Andes:
