all:
	make macroareas
	make counties

macroareas:
	deryaSE --extractadmins --mapfile ../Setup/Kenya_5km.asc --namefield Name --sortfield id --shapeonly Shapefiles/MacroAreas.shp
	mv outputmap.asc  Kenya_5km_maonly.asc
	mv outputmap.ppm  Kenya_5km_maonly.ppm

counties:
	deryaSE --extractadmins --mapfile ../Setup/Kenya_5km.asc --namefield COUNTY --sortfield OBJECTID --export Shapefiles/County.shp
	mv outputmap.asc  Kenya_5km_ids.asc
	mv outputmap.ppm  Kenya_5km_ids.ppm
	mv *.dat Counties/

clean:
	\rm -f *.asc
	\rm -f *.ppm
	\rm -f Counties/*

