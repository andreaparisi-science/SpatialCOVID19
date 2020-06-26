CORES=12
GRIDRES=5
COUNTRY=Kenya

PROCS=$(printf "%03d" ${CORES})

echo -e "\033[36;1mDownloading and installing the engine\033[0m"
cd Extern
bash install.sh $1
if [ $? -eq 1 ]; then
	echo -e "\033[31;1mError in installing external engine.\033[0m"
	exit 1
fi

cd ..
echo -e "\033[36;1mDownloading maps\033[0m"
cd Data/${COUNTRY}
cd Gridded
#bash install.sh

echo -e "\033[36;1mRescaling maps\033[0m"
cd ../Setup
echo CORES:=${CORES} > Makefile
echo GRIDRES:=${GRIDRES} >> Makefile
echo COUNTRY:=${COUNTRY} >> Makefile
echo OUTDIR:=Out${COUNTRY} >> Makefile
echo PROCS:=${PROCS} >> Makefile
cat Makefile_base >> Makefile
#make
#make part

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

echo -e "\033[36;1mSetting up models\033[0m"
cd ../../../../Scenarios/
mkdir -p Default/Sources/
cp ../Main/makeScenario.sh .
cd Default/
cp ../../Main/Makefile_base .
echo CORES:=${CORES} > Makefile
echo GRIDRES:=${GRIDRES} >> Makefile
echo COUNTRY:=${COUNTRY} >> Makefile
cat Makefile_base >> Makefile

cd Sources/
cp ../../../Main/Makefile_make_base .
cp ../../../Main/config.py .
cp ../../../Main/generate.py .
cp ../../../Main/scriptruns.sh .
cp ../../../Main/scriptfits.sh .
cp ../../../Main/*-setup.cpp .
ln -sf ../../../Main/lse-userdefined-base.cpp
ln -sf ../../../Main/Intervention.cpp
ln -sf ../../../Main/Intervention.h
ln -sf ../../../Main/Policy.cpp
ln -sf ../../../Main/Policy.h
ln -sf ../../../Main/SplineInterpolator.cpp
ln -sf ../../../Main/SplineInterpolator.h

echo -e "\033[36;1mReady to build models.  Use make write, make read, make fits or make nofam.\033[0m"

