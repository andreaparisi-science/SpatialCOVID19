CORES:=12
GRIDRES:=20
COUNTRY:=Test

all:
	make counties

counties:
	deryaSE --extractadmins --mapfile ../Setup/${COUNTRY}_${GRIDRES}km.asc --namefield DEN_PROV --sortfield COD_PROV --export Shapefiles/Limiti2011_g/Prov2011_g/Prov2011_WGS.shp
	mv outputmap.asc  ${COUNTRY}_${GRIDRES}km_ids.asc
	mv outputmap.ppm  ${COUNTRY}_${GRIDRES}km_ids.ppm
	mkdir -p Counties
	mv *.dat Counties/

clean:
	\rm -f *.asc
	\rm -f *.ppm
	\rm -f Counties/*

