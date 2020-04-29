CORES=12
GRIDRES=5
COUNTRY=Italy
OUTDIR=OutItaly

PROCS=$(printf "%03d" ${CORES})

echo -e "\033[36;1mDownloading and installing the engine\033[0m"
cd Extern
bash install.sh $1

cd ..
echo -e "\033[36;1mDownloading maps\033[0m"
cd ${COUNTRY}
cd Gridded
bash install.sh

echo -e "\033[36;1mRescaling maps\033[0m"
cd ../Setup
echo CORES:=${CORES} > Makefile
echo GRIDRES:=${GRIDRES} >> Makefile
echo COUNTRY:=${COUNTRY} >> Makefile
echo OUTDIR:=Out${COUNTRY} >> Makefile
echo PROCS:=${PROCS} >> Makefile
cat Makefile_base >> Makefile
make
make part

echo -e "\033[36;1mWorking on maps\033[0m"
cd ../Maps
echo CORES:=${CORES} > Makefile
echo GRIDRES:=${GRIDRES} >> Makefile
echo COUNTRY:=${COUNTRY} >> Makefile
cat Makefile_base >> Makefile
make

echo -e "\033[36;1mSetting up mobility analysis\033[0m"
cd Mobility
echo -e "GENERAL:" > CoV.lse
echo -e "\tPartfile\t\"../../Setup/${COUNTRY}_${GRIDRES}km-${PROCS}.ppm\"" >> CoV.lse
echo -e "\tLocations\t50" >> CoV.lse
echo -e "\tReadinit\tWrite" >> CoV.lse
echo -e "\tStorefile\t\"storage\"" >> CoV.lse
echo -e "\tReduceFactor\t1" >> CoV.lse
echo -e "\tSimlength\t10" >> CoV.lse
echo -e "\tIndivPrefLocs\t0" >> CoV.lse
echo -e "\tIndivHasData\t0" >> CoV.lse
echo -e "\tAccessCycle\t1" >> CoV.lse
echo -e "\nPOPULATIONS:" >> CoV.lse
echo -e "\tDefault \"../../Setup/Kenya_5km.asc\"\tNone\tRadiationNorm\t0.05" >> CoV.lse
cat CoV_base.lse >> CoV.lse

echo CORES:=${CORES} > Makefile
echo GRIDRES:=${GRIDRES} >> Makefile
echo COUNTRY:=${COUNTRY} >> Makefile
cat Makefile_base >> Makefile

echo -e "\033[36;1mBuilding models\033[0m"
cd ../../../Master
echo CORES:=${CORES} > Makefile
echo GRIDRES:=${GRIDRES} >> Makefile
echo COUNTRY:=${COUNTRY} >> Makefile
echo OUTDIR:=${OUTDIR} >> Makefile
cat Makefile_base >> Makefile
make

