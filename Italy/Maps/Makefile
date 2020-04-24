all:
	make counties

counties:
	deryaSE --extractadmins --mapfile ../Setup/Italy_5km.asc --namefield DEN_PROV --sortfield COD_PROV --export Shapefiles/Limiti2011_g/Prov2011_g/Prov2011_WGS.shp
	mv outputmap.asc  Italy_5km_ids.asc
	mv outputmap.ppm  Italy_5km_ids.ppm
	mv *.dat Counties/

clean:
	\rm -f *.asc
	\rm -f *.ppm
	\rm -f Counties/*

