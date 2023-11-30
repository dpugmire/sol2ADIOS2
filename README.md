# sol2ADIOS2

Converts GE's datafiles from CGNS to ADIOS2

solread-hdf5 will read the CGNS data files and output to an ADIOS2 file.


# Building
Build on whoopingcough.

cmake \
 -DCMAKE_BUILD_TYPE=Debug \
 -DADIOS2_DIR=/apps/ADIOS/install_par/lib/cmake/adios2 \
 -DHDF5_DIR=/apps/hdf5/install \
../solread

Build on Andes:



# Running

mpirun -np 1 ./build/solread-hdf5 <nZones> sol.cgns

nZones can be found by running 'h5ls -r <cgns-file>' and looking at how many zones are in the file.
For the small case that Norbert gave me, (sol.cgns), nZones = 5
