CORES:=12
GRIDRES:=20
COUNTRY:=Test
PROCS:=012
 
map:
	\rm -f lse-userdefined.*
	\ln -s ../../../Main/Setup/lse-userdefined.cpp .
	\ln -s ../../../Main/Setup/lse-userdefined.h .
	deryaSE --define GRIDRES=${GRIDRES} --define COUNTRY=\"${COUNTRY}\" --define SHORTCOUNTRY=\"ita\" --define YEAR=\"2020\" --define REDUCE_FACTOR=0.1 CoV.lse
	./runderyaSE 1

part:
	deryaSE --partitionmap --filein ${COUNTRY}_${GRIDRES}km.asc --patches ${PROCS} --connect --seed 7
	mv wmap-final.ppm ${COUNTRY}_${GRIDRES}km-${PROCS}.ppm
	\rm wmap-*

clean:
	\rm -f deryaSE.log
	\rm -f DeryaSE-summary*
	\rm -f DeryaSE-map*
	\rm -f Test_*km*
	\rm -f runderyaSE

