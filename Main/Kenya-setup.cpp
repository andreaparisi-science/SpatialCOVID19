//#include "../Data/Kenya/CSVReader.h"
//#include "../Data/Kenya/CSVReader.cpp"

// Relevant files loaded at execution time
static  std::string  contactMatrixFile = "../../../../Data/Kenya/Contacts/KenyaContactMatrix";
static  std::string  ageGroupEsriFile  = "../../../../Data/Kenya/Setup/Kenya_%dkm_%d.dat";
static  std::string  agePyramidFile    = "../../../../Data/Kenya/Setup/Kenya_%dkm_g%02u_stats.dat";
static  std::string  identifiersFile   = "../../../../Data/Kenya/Maps/Kenya_%dkm_ids.asc";
static  std::string  timeseriesFile    = "../../../../Data/Kenya/Private/Kenya_timeseries.dat";

static  std::string  importedCasesFile = "../../../../Data/Kenya/Private/imported_byId_byAge.dat";

enum  {EXTENT_NATIONAL = 1, EXTENT_COUNTY};
Intervention  travelBan;
std::vector<Intervention>	interventions = {
	Intervention( POLICY_REDUCE_INFLIGHT, EXTENT_NATIONAL,   0,  0, 1.00 ),
	Intervention( POLICY_TRACING_PROB,    EXTENT_NATIONAL,   0,  0, 0.00 )
/*	Intervention( POLICY_TRACING_PROB,    EXTENT_NATIONAL,   0,  0, 1.00 ),
	Intervention( POLICY_TRACING_PROB,    EXTENT_NATIONAL,   1000000,  10, 0.00 ),
	Intervention( POLICY_SOCIALDIST_PROB, EXTENT_NATIONAL,   2, 14, 0.30 ), // 15th first measures (check GOOGLE)
	Intervention( POLICY_SOCIALDIST_PROB, EXTENT_NATIONAL,  16, 21, 0.45 ), // 15th first measures (check GOOGLE)
	Intervention( POLICY_SOCIALDIST_PROB, EXTENT_NATIONAL,  44, 21, 0.35 ), // 15th first measures (check GOOGLE)

	Intervention( POLICY_STAYATHOME_OTH,  EXTENT_NATIONAL,   2, 14, 0.28 ), // Work from home (following GOOGLE and averaging last few days)
	Intervention( POLICY_STAYATHOME_OTH,  EXTENT_NATIONAL,  16, 14, 0.37 ), // Work from home (following GOOGLE and averaging last few days)
	Intervention( POLICY_STAYATHOME_OTH,  EXTENT_NATIONAL,  37, 28, 0.25 ), // Work from home (following GOOGLE and averaging last few days)

	Intervention( POLICY_SCHOOL_CLOSURE,  EXTENT_NATIONAL,   2,  0, 1.00 ), // 15th School closures
	Intervention( POLICY_FAMILY_TRANSMIT, EXTENT_NATIONAL,   2,  0, 1.75 ), //
	Intervention( POLICY_REDUCE_INFLIGHT, EXTENT_NATIONAL,   4,  8, 1.00 ), // 17th Only resident citizens entering (+ self quarantine ?), embargo on 25th
		travelBan = 
	Intervention( POLICY_TRAVELRED_ADMIN, EXTENT_COUNTY,	25,  0, 0.90 ),  //  6th April, travel ban between 4 counties (see below)
	//Intervention( POLICY_TRAVELREDUCTION, EXTENT_NATIONAL,   2, 14, 0.27-0.14 ), //  
	//Intervention( POLICY_TRAVELREDUCTION, EXTENT_NATIONAL,  16,  7, 0.35-0.19 ), //  
	//Intervention( POLICY_TRAVELREDUCTION, EXTENT_NATIONAL,  23, 14, 0.42-0.19 ), //  
	//Intervention( POLICY_TRAVELREDUCTION, EXTENT_NATIONAL,  44, 21, 0.33-0.13 )  //  
	Intervention( POLICY_TRAVELREDUCTION, EXTENT_NATIONAL,   2, 14, (0.27-0.14)/(1-0.14) ),
	Intervention( POLICY_TRAVELREDUCTION, EXTENT_NATIONAL,  16,  7, (0.35-0.19)/(1-0.19) ),
	Intervention( POLICY_TRAVELREDUCTION, EXTENT_NATIONAL,  23, 14, (0.42-0.19)/(1-0.19) ),
	Intervention( POLICY_TRAVELREDUCTION, EXTENT_NATIONAL,  44, 21, (0.33-0.13)/(1-0.13) )
*/};

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
//	POLICY_REDUCE_INFLIGHT,			// Stops external imports


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
int   maxId = 0;
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
	for (int ii = 0; ii < data.ncols; ii++)  {
		for (int jj = 0; jj < data.nrows; jj++)  {
			if (maxId < idsMap[ii][jj])  {
				maxId = idsMap[ii][jj];
			}
		}
	}
}



int  evalDaily()  {
	static bool  needsUpdating = true;

	if (needsUpdating && totalCases > 50)  {
		for (int jj = 0; jj < interventions.size(); jj++)  {
			if (interventions[jj].getEndTime() >= 10000 && simStatus.getTime() >= params.t0+2)  {
				interventions[jj].setActivationTime( simStatus.getTime()-params.t0 );
				interventions[jj].setDuration( 20 );
				needsUpdating = false;
			}
		}
	}
	if (totalCases > 0)  {
		params.tracing = 50.0/totalCases;
		if (params.tracing > 1)  params.tracing = 1.0;
	}
	return -1;
}



int  evalLocally()  {
	return EXTENT_NATIONAL;
}



inline  int getMobilityDuration(double dist)  {
	RandomGenerator *RNG = simStatus.getRandomGenerator();
//	int kk = static_cast<int>(210*dist/1000.0);
//	if (kk > 210)  kk = 210;
//	return  extractFromDistribution( rnd, durationDistr[kk], -1, 0 );
	return static_cast<int>( 1+RNG->get() * 14 );
}




std::vector< std::vector<double> >  importProbs;  
void  loadImportProbs()  {
	std::ifstream  instream( importedCasesFile );
	int nAgeGroups = groups[ groups.size()-1 ] + 1;
	int  ids, age;
	double dummy, prob;

	importProbs.resize( maxId+1 );
	for (int jj = 1; jj <= maxId; jj++)  {
		importProbs[jj].resize( nAgeGroups, 0.0 );
	}
	while (instream.good())  {
		instream >> ids;
		if (instream.eof())  break;
		instream >> age;
		instream >> dummy;
		instream >> dummy;
		instream >> prob;

//std::cout << "***** " << maxId << " " << ids << " "  << age << "\n" << std::flush;
		importProbs[ids][age] = prob;
//		instream.peek();
	}
}


bool  importedCase()  {
	RandomGenerator *RNG = simStatus.getRandomGenerator();
	int classId = simStatus.getCurrentIndividualClass();
	int age = ageNames[ classId ];
	int xx  = simStatus.getX() % simStatus.getMaxX();
	int yy  = simStatus.getY();
	int ids = idsMap[xx][yy];
	if (ids == 0)  simStatus.abort("IDS is ZERO");
	if (RNG->get() < importProbs[ids][age])  {
		return true;
	}
	return false;
}




void  initCountrySpecific()  {
	loadIdentifiers( identifiersFile );
	loadImportProbs();
	//loadFirstInfections();
}




std::vector<int>  lockdownCountyList = {
//:	41, 47, 44, 46
}; // 41-Nairobi, 47-Mombasa, 44-Kilifi, 46-Kwale
bool  checkLockdown(int x0, int y0)  {
	if (interventions.size() > 1 && simStatus.getTime() > travelBan.getActivationTime())  {
		if (std::find(lockdownCountyList.begin(), lockdownCountyList.end(), idsMap[x0][y0]) != lockdownCountyList.end())  {
			return true;
		}
	}
	return false;
}


// FITTING  PROTOTYPE FOR GENERALIZATION
//enum {PARAM_T0 = 0x01, PARAM_R0 = 0x02, PARAM_GAMMA = 0x04, PARAM_TRACING = 0x08};
//enum {DATA_CASES, DATA_SYMPT, DATA_ASYMPT, DATA_DEATHS};
std::vector< int >  inputTable = {DATA_CUMUL_ALL_CASES, DATA_CUMUL_DEATHS};
std::vector< int >  paramTable = {PARAM_T0, PARAM_BETA, PARAM_OMEGA};
std::vector< int >  distsTable = {DATA_DEATHS};

