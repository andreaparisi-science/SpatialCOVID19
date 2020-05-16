COUNTRY=ken
YEAR=2018

for age in '0' '1' '5' '10' '15' '20' '25' '30' '35' '40' '45' '50' '55' '60' '65' '70' '75' '80'
do
	for sex in 'f' 'm'
	do
		wget -c ftp://ftp.worldpop.org.uk/GIS/AgeSex_structures/Global_2000_2020/${YEAR}/${COUNTRY^^}/${COUNTRY}_${sex}_${age}_${YEAR}.tif
	done
done

