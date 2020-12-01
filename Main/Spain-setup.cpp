// Relevant files loaded at execution time
static  std::string  contactMatrixFile = "../../../../Data/Spain/Contacts/SpainContactMatrix";
static  std::string  ageGroupEsriFile  = "../../../../Data/Spain/Setup/Spain_%dkm_%d.asc";
static  std::string  identifiersFile   = "../../../../Data/Spain/Maps/Spain_%dkm_reg.asc";
static  std::string  timeseriesFile    = "../../../../Data/Spain/Spain_timeseries.dat";
static  std::string  fileDeathsByAge   = "../../../../Data/Spain/Private/deathsByAge.dat";

static  std::string  importedCasesFile = "../../../../Data/Spain/Private/distrCases_byId.dat";

Intervention  nationwideLockdown;


// START DATEL 2020-03-05
enum  {EXTENT_LOCAL = 1, EXTENT_PROVINCE, EXTENT_NATIONAL};
std::vector<Intervention>	interventions = {
	//Intervention( POLICY_STAYATHOME_FULL, EXTENT_NATIONAL,  0,  0, 1.00 ),
	Intervention( POLICY_TRACING_PROB,    EXTENT_NATIONAL,  0,  0, 0.00 ),
	Intervention( POLICY_SOCIALDIST_PROB, EXTENT_NATIONAL,  0,  7, 0.64 ), 
	Intervention( POLICY_STAYATHOME_AGE,  EXTENT_NATIONAL,  0,  7, 0.64 ),
	Intervention( POLICY_STAYATHOME_SCH,  EXTENT_NATIONAL,  0,  7, 0.64 ),
	Intervention( POLICY_STAYATHOME_OTH,  EXTENT_NATIONAL,  0,  7, 0.39 ),

	Intervention( POLICY_REDUCE_INFLIGHT, EXTENT_NATIONAL,  5,  6, 0.14 ), // Frontiers closed on 16-March

		nationwideLockdown = // Effective 14-March
	Intervention( POLICY_SCHOOL_CLOSURE,  EXTENT_NATIONAL,  7,  0, 1.00 ),
	Intervention( POLICY_FAMILY_TRANSMIT, EXTENT_NATIONAL,  7,  0, 1.60 ),
	//Intervention( POLICY_REDUCE_INFLIGHT, EXTENT_NATIONAL,  9,  0, 1.00 ),

	Intervention( POLICY_STAYATHOME_FULL, EXTENT_NATIONAL,  7,  0, 1.00 ),
	Intervention( POLICY_SOCIALDIST_PROB, EXTENT_NATIONAL,  7,  7, 0.90 ), 
	Intervention( POLICY_STAYATHOME_AGE,  EXTENT_NATIONAL,  7,  7, 0.90 ),
	Intervention( POLICY_STAYATHOME_SCH,  EXTENT_NATIONAL,  7,  7, 0.90 ),
	Intervention( POLICY_STAYATHOME_OTH,  EXTENT_NATIONAL,  7,  7, 0.66 ),

	Intervention( POLICY_REDUCE_INFLIGHT, EXTENT_NATIONAL, 11,  7, 0.74 ),
	Intervention( POLICY_STAYATHOME_OTH,  EXTENT_NATIONAL, 14, 14, 0.72 ),
	Intervention( POLICY_REDUCE_INFLIGHT, EXTENT_NATIONAL, 18, 14, 0.92 )



//	Intervention( POLICY_SOCIALDIST_PROB, EXTENT_NATIONAL, 44,  7, 0.85 ),  //  Increases in social distancing and isolation follows data from google trends

//''	Intervention( POLICY_TRACING_PROB,    EXTENT_NATIONAL, 52, 10, 1.00 )  //  ~7 April - 10 days increses in testing with new test-isolate policy
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



// Time; Type (1.0 -> Exp, 2.0 -> Inf, 3.0 -> Asy; Lon; Lat;
std::vector< std::vector<double> >  firstInfections; //(1, {0.0, 0.0, 0.0, 9.705, 45.16});

std::vector< std::vector<int> >    idsMap;
std::vector< std::vector<double> > ageIdsPop;
std::vector<double>  idsPop;
int   maxId = 0;
void  loadIdentifiers( const std::string datafile )  {
	int nAgeGroups = groups[ groups.size()-1 ] + 1;
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

	ageIdsPop.resize( nAgeGroups, std::vector<double>(maxId+1, 0.0) );
	idsPop.resize( maxId+1, 0.0 );
	for (int ii = 0; ii < data.ncols; ii++)  {
		for (int jj = 0; jj < data.nrows; jj++)  {
			int ids = idsMap[ii][jj];
			for (int age = 0; age < nAgeGroups; age++)  {
				ageIdsPop[ age ][ ids ] += baseMap[ age ][ii][jj];
				idsPop[ ids ] += baseMap[ age ][ii][jj];
			}
		}
	}
}



int     evalDaily()  {
	if (totalCases > 0)  {
		params.tracing = 200.0/totalCases;
		if (params.tracing > 1)  params.tracing = 1.0;
	}
	return 0;
}



int     evalLocally()  {
	int index = EXTENT_NATIONAL;

	return  index;
}



std::vector< std::vector<double> >  importProbs;  
void  loadImportProbs()  {
	std::ifstream  instream( importedCasesFile );
	int nAgeGroups = groups[ groups.size()-1 ] + 1;
	int  ids, age;
	double dummy, prob;

	importProbs.resize( maxId+1 );
	params.iomega = 0.0;
	for (int jj = 1; jj <= maxId; jj++)  {
		importProbs[jj].resize( nAgeGroups, 0.0 );
	}
	while (instream.good())  {
		instream >> ids;
		if (instream.eof())  break;
		instream >> prob;

//std::cout << "***** " << maxId << " " << ids << " "  << age << "\n" << std::flush;
		for (int age = 0; age < nAgeGroups; age++)  {
			importProbs[ids][age] = prob / idsPop[ ids ];
			params.iomega = std::max( params.iomega, importProbs[ids][age] );
		}
//		instream.peek();
	}
	for (int ids = 1; ids <= maxId; ids++)  {
		for (int age = 0; age < nAgeGroups; age++)  {
			importProbs[ids][age] /= params.iomega;
		}
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
//std::cout << ids << " " << age << " " << maxId << " " << "\n" << std::flush;
	if (RNG->get() < importProbs[ids][age])  {
		return true;
	}
	return false;
}



int  ximportedCase( int nn )  {
	RandomGenerator *RNG = simStatus.getRandomGenerator();
	int classId = simStatus.getCurrentIndividualClass();
	int age = ageNames[ classId ];
	int xx  = simStatus.getX() % simStatus.getMaxX();
	int yy  = simStatus.getY();
	int ids = idsMap[xx][yy];
	if (ids == 0)  simStatus.abort("IDS is ZERO");
	return  RNG->binomial( nn, importProbs[ids][age] );
}



// Mobility: duration of trips.  Here we assume simple commuting and daily work. Parameter will be fitted
inline  double getMobilityDuration(double dist)  {
	return  1;
}



inline  void  initCountrySpecific()  {
	loadIdentifiers( identifiersFile );
	loadImportProbs();
}



// FITTING  PROTOTYPE FOR GENERALIZATION
std::vector< int >  inputTable = {DATA_CUMUL_ALL_CASES, DATA_CUMUL_DEATHS};
//std::vector< int >  paramTable = {PARAM_T0, PARAM_R0, PARAM_OMEGA};
std::vector< int >  paramTable = {PARAM_T0, PARAM_R0, PARAM_GAMMA, PARAM_OMEGA};
std::vector< int >  distsTable = {DATA_DEATHS};
std::vector< int >  printTable = {DATA_DEATHS};
std::vector< int >  contrTable = {POLICY_SOCIALDIST_PROB, POLICY_SCHOOL_CLOSURE, POLICY_TRACING_PROB, POLICY_REDUCE_INFLIGHT};






bool  checkLockdown(int x0, int y0)  {
	return  false;
}


