#!/bin/bash

MPIRUN="mpiexec -np ${NSLOTS}"
REPLICAS=1

#PBS  -l nodes=1:ppm=16

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mem-per-cpu=4751
#SBATCH --time=48:00:00

module purge
module load GCC/7.3.0-2.30 OpenMPI/3.1.1 LibTIFF/4.0.9 Python/3.6.6

cd Base
${MPIRUN} ./runderyaSEwrite --replicate $REPLICAS 1
cd ..

cd Fits
${MPIRUN} ./runderyaSEfit --replicate $REPLICAS 1
cd ..
