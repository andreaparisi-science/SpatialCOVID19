#!/bin/bash

MPIRUN="mpiexec -np ${NSLOTS}"

#PBS  -l nodes=1:ppm=16

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mem-per-cpu=4751
#SBATCH --time=48:00:00

module purge
module load GCC/7.3.0-2.30 OpenMPI/3.1.1 LibTIFF/4.0.9 Python/3.6.6

cd Base
${MPIRUN} ./runderyaSEwrite 1
cd ..


for jj in {1..100}
do
	if [ ! -e ./EXEC ]; then
		echo File EXEC not found. Last attempted seed [${jj}]
		echo Use 'touch EXEC' to allow execution.
		exit 1
	fi
	printf -v kk "%03d" $jj
	mkdir -p Run-${kk}
	cd Run-${kk}
	${MPIRUN} ../runderyaSEread ${kk}
	cd ..
done

