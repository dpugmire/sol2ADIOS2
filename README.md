# sol2ADIOS2

Converts GE's datafiles from CGNS to ADIOS2

solread-hdf5 will read the CGNS data files and output to an ADIOS2 file.


# Building
Build on whoopingcough.

cmake \
 -DCMAKE_BUILD_TYPE=Debug \
 -DADIOS2_DIR=/apps/ADIOS/install_par/lib/cmake/adios2 \
 -DHDF5_DIR=/apps/hdf5/install \
../src

Build on Andes:

module load hdf5

cmake \
 -DCMAKE_C_COMPILER:FILEPATH=/sw/andes/gcc/9.3.0/bin/gcc \
 -DCMAKE_CXX_COMPILER:FILEPATH=/sw/andes/gcc/9.3.0/bin/g++ \
 -DADIOS2_DIR=/ccs/home/pugmire/software/andes/adios2/2.9.2/install/lib64/cmake/adios2 \
../src


# Running

mpirun -np 1 ./build/solread-hdf5 <nZones> sol.cgns

nZones can be found by running 'h5ls -r <cgns-file>' and looking at how many zones are in the file.
For the small case that Norbert gave me, (sol.cgns), nZones = 5


Andes:

#!/bin/bash
#SBATCH -A CFD164
#SBATCH -J cgns_convert
#SBATCH -N 16
#SBATCH -t 0:20:00

date
srun -n 16 ./cgns2adios --dumpFides output_vki --output output_vki --partTextFile ./vki.txt



## vti.txt contains the parts files, 1 per line.
