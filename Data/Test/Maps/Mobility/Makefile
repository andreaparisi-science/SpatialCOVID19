CORES:=12
GRIDRES:=20
COUNTRY:=Test

all:
	deryaSE --define GRIDRES=${GRIDRES} CoV.lse
	mpirun -np ${CORES} ./runderyaSE 1

clean:
	\rm -f runderyaSE
	\rm -f storage*
	\rm -f pdf.dat
	\rm -f histo*.dat
	\rm -f deryaSE.log
	\rm -f DeryaSE-summary*

