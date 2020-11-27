CORES:=12
GRIDRES:=5
COUNTRY:=Kenya
OUTPUT  := ../Output
CWD := $(shell pwd)
OPTS = $(eval OPTS := $$(shell bash testTauleap.sh))$(OPTS)

all:
	make write
	make read
	make fits

prepare:
	\rm -rvf ${OUTPUT}
	\rm -f runderyaSE*

write:
	cd Sources && \
	\rm -f ./.tauleap && \
	g++ -c SplineInterpolator.cpp && \
	g++ -c Policy.cpp && \
	g++ -c Intervention.cpp && \
	python3 generate.py ${COUNTRY} Write > CoV.lse && \
	echo GRIDRES:=${GRIDRES}  > Makefile_make && \
	echo COUNTRY:=${COUNTRY} >> Makefile_make && \
	cat Makefile_make_base   >> Makefile_make && \
	make -f Makefile_make  deryaSEnorm && \
	mkdir -p ${OUTPUT} && \
	mkdir -p ${OUTPUT}/Base && \
	mv runderyaSE ${OUTPUT}/Base/runderyaSEwrite

runs:
	cd Sources && \
	\rm -f ./.tauleap && \
	g++ -c SplineInterpolator.cpp && \
	g++ -c Policy.cpp && \
	g++ -c Intervention.cpp && \
	python3 generate.py ${COUNTRY} Read > CoV.lse && \
	echo GRIDRES:=${GRIDRES}  > Makefile_make && \
	echo COUNTRY:=${COUNTRY} >> Makefile_make && \
	cat Makefile_make_base   >> Makefile_make && \
	make -f Makefile_make  deryaSEnorm && \
	mkdir -p ${OUTPUT} && \
	mv runderyaSE ${OUTPUT}/runderyaSEread && \
	cp scriptruns.sh ${OUTPUT}/ && \
	touch ${OUTPUT}/EXEC

fits:
	cd Sources && \
	\rm -f ./.tauleap && \
	g++ -c SplineInterpolator.cpp && \
	g++ -c Policy.cpp && \
	g++ -c Intervention.cpp && \
	python3 generate.py ${COUNTRY} Read > CoV.lse && \
	echo GRIDRES:=${GRIDRES}  > Makefile_make && \
	echo COUNTRY:=${COUNTRY} >> Makefile_make && \
	cat Makefile_make_base   >> Makefile_make && \
	make -f Makefile_make  deryaSEfits && \
	mkdir -p ${OUTPUT} && \
	mkdir -p ${OUTPUT}/Fits && \
	mv runderyaSE ${OUTPUT}/Fits/runderyaSEfit && \
	cp scriptfits.sh ${OUTPUT}/

sources:
	mkdir -p ${OUTPUT}/Sources
	cp CoV.lse ${OUTPUT}/Sources
	cp lse-userdefined.cpp ${OUTPUT}/Sources
	cp config.py ${OUTPUT}/Sources
	cp generate.py ${OUTPUT}/Sources
	cp lse-userdefined-base.cpp ${OUTPUT}/Sources
	cp ${COUNTRY}-setup.cpp ${OUTPUT}/Sources

clean:
	\rm -f ./.tauleap
	\rm -f CoV.lse
	\rm -f deryaSE.log
	\rm -f *.o
	\rm -rvf __pycache__

