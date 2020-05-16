#!/bin/python3

import os
import os.path

os.system( "cd Base; mpirun -np 12 ./runderyaSEwrite 1; cd .." )
os.system( "cd Fits;  mpirun -np 12 ./runderyaSEfit --replicate 1 1; cd .." )
