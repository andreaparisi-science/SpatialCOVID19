COUNTRY=$1
GRIDRES=$2
NUMBER=0

deryaSE --extractmap ${COUNTRY}_${GRIDRES}km.asc 210 540 200 388 1 1 0 1 > PPP_5km.asc
mv PPP_${GRIDRES}km.asc ${COUNTRY}_${GRIDRES}km.asc
while [ ${NUMBER} -le 8 ] ; do \
	deryaSE --extractmap ${COUNTRY}_${GRIDRES}km_${NUMBER}.asc 210 540 200 388 1 1 0 1 > PPP_${GRIDRES}km_${NUMBER}.asc ; \
	mv PPP_${GRIDRES}km_${NUMBER}.asc ${COUNTRY}_${GRIDRES}km_${NUMBER}.asc ; \
	NUMBER=$(( ${NUMBER} + 1 )) ; \
done

