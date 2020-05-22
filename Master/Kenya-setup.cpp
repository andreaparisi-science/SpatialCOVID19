#include "../Data/Kenya/CSVReader.h"
#include "../Data/Kenya/CSVReader.cpp"

// Relevant files loaded at execution time
static  std::string  contactMatrixFile = "../../../Data/Kenya/Contacts/KenyaContactMatrix";
static  std::string  ageGroupEsriFile  = "../../../Data/Kenya/Setup/Kenya_%dkm_%d.dat";
static  std::string  identifiersFile   = "../../../Data/Kenya/Maps/Kenya_%dkm_ids.asc";
static  std::string  timeseriesFile    = "../../../Data/Kenya/Kenya_timeseries.dat";

enum  {EXTENT_NATIONAL = 1};
std::vector<Intervention>	interventions = {
	Intervention( POLICY_SOCIALDIST_PROB, EXTENT_NATIONAL,   1000000.0, 14, 0.75 ),
	Intervention( POLICY_SCHOOL_CLOSURE,  EXTENT_NATIONAL,   1000000.0,  0, 1.00 ),
	Intervention( POLICY_FAMILY_TRANSMIT, EXTENT_NATIONAL,   1000000.0,  0, 1.75 ),
	Intervention( POLICY_STAYATHOME_OTH,  EXTENT_NATIONAL,   1000000.0, 14, 0.30 ),
	Intervention( POLICY_STAYATHOME_AGE,  EXTENT_NATIONAL,   1000000.0, 14, 0.90 ),
	Intervention( POLICY_STAYATHOME_SCH,  EXTENT_NATIONAL,   1000000.0, 14, 0.80 ),
	Intervention( POLICY_STAYATHOME_FULL, EXTENT_NATIONAL,   1000000.0,  0, 1.00 )
};

//	POLICY_TRACING_PROB = 0,		// Contact tracing probability
//	POLICY_SOCIALDIST_PROB, 		// Generalized reduction of social interactions
//	POLICY_TRAVELREDUCTION, 		// Reduction of travel intensity
//	POLICY_TRAVELRED_ADMIN, 		// Reduction of travel intensity
//	POLICY_STAYATHOME_AGE, 			// Compliance of stay at home for oldest (non-working)
//	POLICY_STAYATHOME_OTH, 			// Compliance of stay at home for working individuals
//	POLICY_STAYATHOME_SCH, 			// Compliance of stay at home for school-aged individuals
//	POLICY_FAMILY_TRANSMIT, 		// Increae in family transmission
//	POLICY_STAYATHOME_FULL, 		// Whether stay-at-home for younger and older means avoiding all social contacts (ex. no shopping at all)
//	POLICY_SCHOOL_CLOSURE, 			// Fraction of schools closed (generalized)


std::vector< std::vector<double> >  firstInfections = {
	{0.0, 0.0, -1.0, 36.8172, -1.2864}//,  // First absolute
/*	{3.0, 3.0, 36.8172, -1.2864},  // 5th of March (x3)
	{3.0, 3.0, 36.8172, -1.2864},  // 5th of March
	{3.0, 3.0, 36.8172, -1.2864},  // 5th of March
	{7.0, 3.0, 36.8172, -1.2864},  // 9th of March
	{12.0, 3.0, 36.8172, -1.2864},  // 17th of March
	{17.0, 3.0,  36.8172, -1.2864}, // 22nd of March (x8)
	{17.0, 3.0,  36.8172, -1.2864}, // 22nd of March
	{17.0, 3.0,  36.8172, -1.2864}, // 22nd of March
	{17.0, 3.0,  36.8172, -1.2864}, // 22nd of March
	{17.0, 3.0,  36.8172, -1.2864}, // 22nd of March
	{17.0, 3.0,  36.8172, -1.2864}, // 22nd of March
	{17.0, 3.0,  36.8172, -1.2864}, // 22nd of March
	{17.0, 3.0,  36.8172, -1.2864}  // 22nd of March
*/};
// Note Day 0, the 13th is day 8


std::vector< std::vector<int> >    idsMap;
void  loadIdentifiers( const std::string datafile )  {
	char ch_filename[256];
	std::string  filename;
	sprintf( ch_filename, datafile.c_str(), GRIDRES );
	filename = std::string( ch_filename );
	EsriReader  reader(filename);
	MapData data = reader.readHeader();
	idsMap.resize( data.ncols );
	for (int jj = 0; jj < data.ncols; jj++)  idsMap[jj].resize( data.nrows, 0 );
	reader.readFile( 1.0, idsMap );
}



// Implements measures taking into account their spatial extent
void	evalLocalParameters()  {
	static int  prev_index = -1;
	static int  prev_day = -1;
	static bool  needsUpdating = true;



	int day = static_cast<int>(simStatus.getTime()+0.0000001);
	if (day != prev_day)  {
		prev_day = day;
		prev_index = -1;

//  	THIS TRIGGERS APPLICATION OF POLICY AFTER CROSSING 100 CASES PER DAY
		if (needsUpdating && totalCases > 100)  {
			for (int jj = 0; jj < interventions.size(); jj++)  {
				if (interventions[jj].getEndTime() >= 10000)  {
					interventions[jj].setActivationTime( simStatus.getTime() );
				}
			}
			needsUpdating = false;
		}
//
	}

	int index = 1;
	// To be applied only if required, otherwise values are already correct
	if (prev_index != index)  {
		if (index != 0)  {
			TRACING_PROB[0]    = TRACING_PROB[index]*params.tracing;
			SOCIALDIST_PROB[0] = SOCIALDIST_PROB[index];
			TRAVELREDUCTION[0] = TRAVELREDUCTION[index];
			TRAVELRED_ADMIN[0] = TRAVELRED_ADMIN[index];
			STAYATHOME_AGE[0]  = STAYATHOME_AGE[index];
			STAYATHOME_OTH[0]  = STAYATHOME_OTH[index];
			STAYATHOME_SCH[0]  = STAYATHOME_SCH[index];
			FAMILY_TRANSMIT[0] = FAMILY_TRANSMIT[index];
			STAYATHOME_FULL[0] = STAYATHOME_FULL[index];
			SCHOOL_CLOSURE[0]  = SCHOOL_CLOSURE[index];
		}

#ifndef  MODEL_FAMILY
		params.home		= FAMILY_TRANSMIT[0];
#endif
		params.work		= 1.0;
		params.school	= 1.0-SCHOOL_CLOSURE[0];
		params.other	= 1.0-SOCIALDIST_PROB[0];
		params.mobility = 1.0-TRAVELREDUCTION[0];
		updateContactMatrix();
		prev_index = index;
	}
}






inline  int getMobilityDuration(double dist)  {
	RandomGenerator *RNG = simStatus.getRandomGenerator();
//	int kk = static_cast<int>(210*dist/1000.0);
//	if (kk > 210)  kk = 210;
//	return  extractFromDistribution( rnd, durationDistr[kk], -1, 0 );
	return static_cast<int>( 1+RNG->get() * 14 );
}



/*
void  loadFirstInfections()  {
	std::ifstream  instream("county-level-data-byId.dat");
	CSVReader csv( ',', '"' );
	csv.removeQuotes( true );
	csv.trimSpaces( true );
	char line[65536];
	std::vector<std::string>  vec;

	// Build county popdistr;
	std::vector< std::vector<double> >  countyPopDistr();  // double allows using extractFromDistribution
	int maxCountyId = 0;
	for (int ii = 0; ii < simStatus.getMaxX(); ii++)  {
		for (int jj = 0; jj < simStatus.getMaxY(); jj++)  {
			if (popMap[xx][yy] > 0)  {
				id = idsMap[ii][jj];
				if (maxCountyId < id)  {
					maxCountyId = id;
					countyPopDistr.resize( maxCountyId+1 );
					if (countyPopDistr[id].size() == 0)  {
						countyPopDistr[id].push_back( popMap[xx][yy] );
					} else {
						countyPopDistr[id].push_back( popMap[xx][yy] + countyPopDistr[id].back() );
				}
			}
		}
	}


	while (instream.good())  {
		double time;
		int target, agegroup, importedSymp, countyId;
		instream.getline( line, 65536 );
		vec = csv.getline( line );
		time = atof(vec[1]);
		agegroup = atol(vec[2]);
		importedSymp = vec[3];
		countyId = vec[12];
		instream.peek();

		target = static_cast<int>( RNG->get()*countyPopDistr[countyId].back() );
		double rnd = RNG->get();
		target = extractFromDistribution( rnd, countyPopDistr[countyId], countyPopDistr[countyId].size(), 1 );
		int count = 0;
		for (int ii = 0; ii < simStatus.getMaxX(); ii++)  {
			for (int jj = 0; jj < simStatus.getMaxY(); jj++)  {
				if (popMap[xx][yy] > 0)  {
					id = idsMap[ii][jj];
					if (id == countyId)  {
						if (count == target)  {
							lon = simStatus.getLongitude( ii );
							lat = simStatus.getLatitude( jj );
							for (int kk = 0; kk < importedSymp; kk++)  {
								firstInfection.push_back( {time, 2.0, agegroup, lon, lat} );
							}
						}
					}
				}
			}
		}		
	}
}
*/


void  initCountrySpecific()  {
	loadIdentifiers( identifiersFile );
	//loadFirstInfections();
}




std::vector<int>  lockdownCountyList = {
	41, 47, 44, 46
}; // 41-Nairobi, 47-Mombasa, 44-Kilifi, 46-Kwale
bool  checkLockdown(int x0, int y0)  {
//	if (policy.size() > 1 && simStatus.getTime() > policyTime[1])  {
//		if (std::find(lockdownCountyList.begin(), lockdownCountyList.end(), idsMap[x0][y0]) != lockdownCountyList.end())  {
//			return true;
//		}
//	}
	return false;
}


// FITTING  PROTOTYPE FOR GENERALIZATION
//enum {PARAM_T0 = 0x01, PARAM_R0 = 0x02, PARAM_GAMMA = 0x04, PARAM_TRACING = 0x08};
//enum {DATA_CASES, DATA_SYMPT, DATA_ASYMPT, DATA_DEATHS};
std::vector< int >  inputTable = {DATA_DUMMY, DATA_CASES, DATA_DEATHS};
std::vector< int >  paramTable = {PARAM_T0, PARAM_R0, PARAM_GAMMA};
std::vector< int >  distsTable = {DATA_CASES};

