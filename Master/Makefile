CORES:=12
GRIDRES:=5
COUNTRY:=Kenya
OUTDIR:=OutKenya3
OUTPUT  := ../Output/${OUTDIR}
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
	\rm -f ./.tauleap
	g++ -c SplineInterpolator.cpp
	g++ -c Policy.cpp
	python3 generate.py ${COUNTRY} Write > CoV.lse
	echo GRIDRES:=${GRIDRES}  > Makefile_make
	echo COUNTRY:=${COUNTRY} >> Makefile_make
	cat Makefile_make_base >> Makefile_make
	make -f Makefile_make  deryaSEnorm
	mkdir -p ${OUTPUT}
	mkdir -p ${OUTPUT}/Base
	mv runderyaSE runderyaSEwrite
	mv runderyaSEwrite ${OUTPUT}/Base

read:
	\rm -f ./.tauleap
	g++ -c SplineInterpolator.cpp
	g++ -c Policy.cpp
	python3 generate.py ${COUNTRY} Read > CoV.lse
	echo GRIDRES:=${GRIDRES}  > Makefile_make
	echo COUNTRY:=${COUNTRY} >> Makefile_make
	cat Makefile_make_base >> Makefile_make
	make -f Makefile_make  deryaSEnorm
	mkdir -p ${OUTPUT}
	mv runderyaSE ${OUTPUT}/runderyaSEread
	cp multirun.py ${OUTPUT}/
	touch ${OUTPUT}/EXEC

fits:
	\rm -f ./.tauleap
	g++ -c SplineInterpolator.cpp
	g++ -c Policy.cpp
	python3 generate.py ${COUNTRY} Read > CoV.lse
	echo GRIDRES:=${GRIDRES}  > Makefile_make
	echo COUNTRY:=${COUNTRY} >> Makefile_make
	cat Makefile_make_base >> Makefile_make
	make -f Makefile_make  deryaSEfits
	mkdir -p ${OUTPUT}
	mkdir -p ${OUTPUT}/Fits
	mv runderyaSE ${OUTPUT}/Fits/runderyaSEfit
	cp fitting.py ${OUTPUT}/

clean:
	\rm -f ./.tauleap
	\rm -f CoV.lse
	\rm -f deryaSE.log
	\rm -f *.o
	\rm -rvf __pycache__

