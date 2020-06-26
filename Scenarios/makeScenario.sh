if [ $# == 1 ]; then
	BASE=Default
else
	BASE=$2
fi


echo Generating Scenario [$1] from [${BASE}]
mkdir -p $1
cp -rvL ${BASE}/Sources/  $1/ &> /dev/null
cp -v ${BASE}/Make*  $1/ &> /dev/null
mkdir -p $1/Output/

