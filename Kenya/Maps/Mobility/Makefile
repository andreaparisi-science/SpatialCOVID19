CORES := 12

all:
	deryaSE CoV.lse
	mpirun -np 12 ./runderyaSE 1

clean:
	\rm -f runderyaSE
	\rm -f storage*
	\rm -f pdf.dat
	\rm -f histo*.dat
	\rm -f deryaSE.log
	\rm -f DeryaSE-summary*

