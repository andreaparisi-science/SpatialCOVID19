#!/bin/python3

import os
import os.path

os.system( "cd Base; mpirun -np 12 ./runderyaSEwrite 1; cd .." )
for jj in range(0,100):
	kk = jj+1
	if not os.path.isfile("EXEC"):
		print( "File EXEC not found. Last known seed [" + str(kk-1) + "]\n" )
		print( "Use 'touch EXEC' to allow execution.\n" )
		break
	os.system( "mkdir -p Run-{:03d}; cd ..".format(kk) )
#	os.system( "cd Run-{:03d}; srun ../runderyaSEread {:d}; cd ../..".format(kk, kk) )
	os.system( "cd Run-{:03d}; mpirun -np 12 ../runderyaSEread {:d}; cd ../..".format(kk, kk) )
