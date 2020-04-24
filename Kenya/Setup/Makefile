PROCS := 012
COUNTRY := Kenya

map:
	deryaSE CoV.lse
	./runderyaSE 1

part:
	deryaSE --partitionmap --filein ${COUNTRY}_5km.asc --patches ${PROCS} --connect --seed 7
	mv wmap-final.ppm ${COUNTRY}_5km-${PROCS}.ppm
	\rm wmap-*

clean:
	\rm -f deryaSE.log
	\rm -f DeryaSE-summary*
	\rm -f DeryaSE-map*
	\rm -f Italy_5km*
	\rm -f runderyaSE

