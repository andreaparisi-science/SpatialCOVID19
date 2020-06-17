import urllib.request
import sys
import csv
import os.path

countryCodes = {
	'Spain':	'ES',
	'Italy':	'IT',
	'Germany':	'GM',
	'Portugal':	'PT',
	'Switzerland':	'CH',
	'Canada':		'CA',
	'UK':		'GB'
}

country  = sys.argv[1]
code = countryCodes[country]

filename = country+".csv"

if os.path.isfile( filename ):
	print( "Data exists" )
else:
	url = "https://api.census.gov/data/timeseries/idb/1year?get=AGE,POP&FIPS=" + code + "&time=2020&SEX=0"
	response = urllib.request.urlopen( url )
	dataframe = response.read().decode('UTF-8')
	dataframe = eval(dataframe)
	header=dataframe[0]
	data=dataframe[1:]
	# Format is Age | Pop | Country | Year | Sex

	with open( filename, 'w', newline='' ) as csvfile:
		writer = csv.writer(csvfile, delimiter=' ',
				quotechar='"', quoting=csv.QUOTE_MINIMAL)
		writer.writerows(data)
	print ( "Data retrieved" )


