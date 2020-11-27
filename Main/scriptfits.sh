#!/bin/bash
# Delete the lines that do not correspond to your cluster manager
# and alter according to your needs

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mem-per-cpu=4571
#SBATCH --time=48:00:00

#PBS  -l nodes=1:ppm=16

module purge
module load GCC/7.3.0-2.30 OpenMPI/3.1.1 LibTIFF/4.0.9 Python/3.6.6

MPIRUN="mpiexec -np ${NSLOTS}"

REPLICAS=1

cd Base
${MPIRUN} ./runderyaSEwrite --replicate $REPLICAS 1
cd ..

cd Fits
${MPIRUN} ./runderyaSEfit --replicate $REPLICAS 1
cd ..
