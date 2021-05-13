// Relevant files loaded at execution time
static  std::string  contactMatrixFile    = "../../../../Data/UK/Contacts/UKContactMatrix";
static  std::string  ageGroupEsriFile     = "../../../../Data/UK/Setup/UK_%dkm_%d.asc";
//static  std::string  identifiers_NHS_File = "../../../../Data/UK/Maps/UK_%dkm_NHS_ids.asc";
static  std::string  identifiers_LA_File  = "../../../../Data/UK/Maps/UK_%dkm_2019_ids.asc";
static  std::string  timeseriesFile       = "../../../../Data/UK/UK_timeseries.dat";
// TO BE MOFIDIED APPROPRIATELY
static  std::string  fileDeathsByAge      = "../../../../Data/Italy/Shared/deathsByAge.dat";

static  std::string  imported_LA_File     = "../../../../Data/UK/Shared/distrCases_byId.dat";
//static  std::string  imported_NHS_File    = "../../../../Data/UK/Shared/distrCases_byRegion_byAge.dat";

Intervention  nationwideLockdown;

// UK timeline
// 31st Jan: first detected case. Next two cases in early Feb local transmission
// 5th of March first death
//
// GVT and NON GVT ACTIONS
//
// 12 March (+41): PHE stops doing contact tracing (overwhelmed capacity)
//   Self isolation
// 13 March (+42): Premier league suspended. Some events suspended
// 16 March (+45): non essential travel and contacts "discouraged"
// 17 March (+46): non essential international travel discouraged
// 21 March (+50) (from): school close
// 21 March (+50) (from): closing all food estabilishments
// 26 March (+55): partial lockdown
// 25 March (+54): discourage non-essential journey
// 28 March (+57): regulations on in NOrthern Ireland (ask relationship betweeh NI and UK
//  3 April (+63): 73% reduction in travel mobility
// 10 May   (+100): Udaptes from stay-at-home to stay-safe ???



enum  {EXTENT_LOCAL = 1, EXTENT_PROVINCE, EXTENT_NATIONAL};
std::vector<Intervention>	interventions = {
	Intervention( POLICY_TRACING_PROB,    EXTENT_NATIONAL,   100000,  0, 1.00),
	//Intervention( POLICY_STAYATHOME_FULL, EXTENT_NATIONAL,   0,  0, 1.00 ), // Any stay at hom treated as full
	// Tracing suspended after 12 of March
	Intervention( POLICY_TRACING_PROB,    EXTENT_NATIONAL,   6,  0, 0.00 ),
	// Reduction of non essential flights to UK -> Pre-reduction from around the 2 of March. Maximum after 5 weeks
	Intervention( POLICY_REDUCE_INFLIGHT, EXTENT_NATIONAL,  -4, 35, 1.00 ),
	// 21 March
	Intervention( POLICY_SCHOOL_CLOSURE,  EXTENT_NATIONAL,  15,  0, 0.99 ), // 1% of school children attending
	Intervention( POLICY_FAMILY_TRANSMIT, EXTENT_NATIONAL,  15,  0, 2.00 ),
	// 10 March
	Intervention( POLICY_STAYATHOME_OTH,  EXTENT_NATIONAL,   4, 14, 0.62 ),  
	Intervention( POLICY_STAYATHOME_AGE,  EXTENT_NATIONAL,   4, 14, 0.77 ),
	Intervention( POLICY_STAYATHOME_SCH,  EXTENT_NATIONAL,   4, 14, 0.77 ),
	Intervention( POLICY_STAYATHOME_FULL, EXTENT_NATIONAL,   4,  0, 1.00 ),
// 24 March
	Intervention( POLICY_STAYATHOME_OTH,  EXTENT_NATIONAL,  18, 14, 0.69 ),  
	//  7 April
	Intervention( POLICY_STAYATHOME_OTH,  EXTENT_NATIONAL,  32, 14, 0.60 ),  
	Intervention( POLICY_SCHOOL_CLOSURE,  EXTENT_NATIONAL,  43,  0, 0.98 ),
	Intervention( POLICY_SCHOOL_CLOSURE,  EXTENT_NATIONAL,  78,  0, 0.94 )
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

std::vector< std::vector<int> >    idsMap_LA, idsMap_NHS, idsMap;
std::vector< std::vector<double> > ageIdsPop_LA, ageIdsPop_NHS;
std::vector<double>  idsPop_LA, idsPop_NHS;
int   maxId_LA = 0, maxId_NHS = 0, maxId = 0;
//void  loadIdentifiers( const std::string datafile_LA, const std::string datafile_NHS )  {
void  loadIdentifiers( const std::string datafile_LA )  {
	char ch_filename[256];
	std::string  filename;
	int ids, nAgeGroups = groups[ groups.size()-1 ] + 1;

	sprintf( ch_filename, datafile_LA.c_str(), GRIDRES );
	filename = std::string( ch_filename );
	EsriReader  reader_LA(filename);
	MapData data = reader_LA.readHeader();
	idsMap_LA.resize( data.ncols );
	for (int jj = 0; jj < data.ncols; jj++)  idsMap_LA[jj].resize( data.nrows, 0 );
	reader_LA.readFile( 1.0, idsMap_LA );
	for (int ii = 0; ii < data.ncols; ii++)  {
		for (int jj = 0; jj < data.nrows; jj++)  {
			if (maxId_LA < idsMap_LA[ii][jj])  {
				maxId_LA = idsMap_LA[ii][jj];
			}
		}
	}

	ageIdsPop_LA.resize( nAgeGroups, std::vector<double>(maxId_LA+1, 0.0) );
	idsPop_LA.resize( maxId_LA+1, 0.0 );
	for (int ii = 0; ii < data.ncols; ii++)  {
		for (int jj = 0; jj < data.nrows; jj++)  {
			ids = idsMap_LA[ii][jj];
			for (int age = 0; age < nAgeGroups; age++)  {
				ageIdsPop_LA[ age ][ ids ] += baseMap[ age ][ii][jj];
				idsPop_LA[ ids ] += baseMap[ age ][ii][jj];
			}
		}
	}
/*
	sprintf( ch_filename, datafile_NHS.c_str(), GRIDRES );
	filename = std::string( ch_filename );
	EsriReader  reader_NHS(filename);
	data = reader_NHS.readHeader();
	idsMap_NHS.resize( data.ncols );
	for (int jj = 0; jj < data.ncols; jj++)  idsMap_NHS[jj].resize( data.nrows, 0 );
	reader_NHS.readFile( 1.0, idsMap_NHS );
	for (int ii = 0; ii < data.ncols; ii++)  {
		for (int jj = 0; jj < data.nrows; jj++)  {
			if (maxId_NHS < idsMap_NHS[ii][jj])  {
				maxId_NHS = idsMap_NHS[ii][jj];
			}
		}
	}

	ageIdsPop_NHS.resize( nAgeGroups, std::vector<double>(maxId_LA+1, 0.0) );
	idsPop_NHS.resize( maxId_LA+1, 0.0 );
	for (int ii = 0; ii < data.ncols; ii++)  {
		for (int jj = 0; jj < data.nrows; jj++)  {
			ids = idsMap_NHS[ii][jj];
			for (int age = 0; age < nAgeGroups; age++)  {
				ageIdsPop_NHS[ age ][ ids ] += baseMap[ age ][ii][jj];
				idsPop_NHS[ ids ] += baseMap[ age ][ii][jj];
			}
		}
	}*/
	idsMap = idsMap_LA;
	maxId = maxId_LA;
}



int     evalDaily()  {
	static bool  needsUpdating = true;

	if (needsUpdating)  {
		for (int jj = 0; jj < interventions.size(); jj++)  {
			if (interventions[jj].getEndTime() >= 10000)  {
				interventions[jj].setActivationTime( simStatus.getTime()-params.t0 );
				interventions[jj].setDuration( 1 );
				needsUpdating = false;
			}
		}
	}

	if (totalCases > 0)  {
		params.tracing = 200.0/totalCases;
		if (params.tracing > 1)  params.tracing = 1.0;
	}
	return 0;
}


int     evalLocally()  {
	int index = EXTENT_NATIONAL;
	// Presume to apply nationwide values;
	int xx = simStatus.getX() % simStatus.getMaxX();
	int yy = simStatus.getY();
	return  index;
}



std::vector< std::vector<double> >  importProbs_byLA;
//std::vector< std::vector<double> >  importProbs_byNHS;
std::vector< std::vector<double> >  importProbs;
void  loadImportProbs()  {
	int nAgeGroups = groups[ groups.size()-1 ] + 1;
	int  ids, age;
	double dummy, prob;

	importProbs_byLA.resize( maxId_LA+1 );
	for (int jj = 1; jj <= maxId_LA; jj++)  {
		importProbs_byLA[jj].resize( nAgeGroups, 0.0 );
	}
//	importProbs_byNHS.resize( maxId_NHS+1 );
//	for (int jj = 1; jj <= maxId_NHS; jj++)  {
//		importProbs_byNHS[jj].resize( nAgeGroups, 0.0 );
//	}
	importProbs.resize( maxId+1 );
	for (int jj = 1; jj <= maxId; jj++)  {
		importProbs[jj].resize( nAgeGroups, 0.0 );
	}

	std::ifstream  instream_LA( imported_LA_File );
	while (instream_LA.good())  {
		instream_LA >> ids;
		if (instream_LA.eof())  break;
		instream_LA >> prob;

//std::cout << "***** " << maxId_2019 << " " << ids << " "  << age << "\n" << std::flush;
		for (age = 0; age < nAgeGroups; age++)  {
			importProbs_byLA[ids][age] = prob / idsPop_LA[ids];
		}
	}
	instream_LA.close();

/*	std::ifstream  instream_NHS( imported_NHS_File );
	while (instream_NHS.good())  {
		instream_NHS >> ids;
		instream_NHS >> age;
		if (instream_NHS.eof())  break;
		instream_NHS >> prob;

//std::cout << "***** " << maxId_2019 << " " << ids << " "  << age << "\n" << std::flush;
		importProbs_byNHS[ids][age] = prob;
	}
	instream_NHS.close();

	// Set data for region 8 (Scotland). We assume an age distribution based 
	for (int jj = 1; jj <= maxId_NHS; jj++)  {
		if (jj == 8) continue;
		for (int age = 0; age < nAgeGroups; age++)  {
			importProbs_byNHS[8][age] += importProbs_byNHS[jj][age] * ageIdsPop_NHS[age][jj];
		}
	}
	double extPop = 0;  for (int jj = 1; jj < 8; jj++)  extPop += idsPop_NHS[jj];
	extPop += idsPop_NHS[9];
	for (int age = 0; age < nAgeGroups; age++)  {
		importProbs_byNHS[8][age] /= ageIdsPop_NHS[age][8];
		importProbs_byNHS[8][age] *= (idsPop_NHS[8]/extPop);
	}
*/
	// This is to set up the initial condition
	params.iomega = 0;
	for (int xx = 0; xx < simStatus.getMaxX(); xx++)  {
		for (int yy = 0; yy < simStatus.getMaxY(); yy++)  {
			if (popMap[xx][yy] > 0)  {
				int ids_LA  = idsMap_LA[xx][yy];
//				int ids_NHS = idsMap_NHS[xx][yy];
//				importProbs[ids_LA][age] = importProbs_byLA[ids_LA][age] * importProbs_byNHS[ids_NHS][age];
//				params.iomega = std::max(params.iomega, importProbs_byLA[ids_LA][age] * importProbs_byNHS[ids_NHS][age]);
				params.iomega = std::max(params.iomega, importProbs_byLA[ids_LA][age]);
			}
		}
	}

	for (int ids = 1; ids <= maxId_LA; ids++)  {
		for (int age = 0; age < nAgeGroups; age++)  {
			importProbs_byLA[ids][age] /= params.iomega;
		}
	}
	
/*	for (int jj = 1; jj <= maxId_NHS; jj++)  {
		std::cout << "[[[" << jj << "]]] ";
		for (int age = 0; age < nAgeGroups; age++)  {
			std::cout << importProbs_byNHS[jj][age] << "\t";
		}
		std::cout << "\n";
	}*/
	importProbs = importProbs_byLA;
}




bool  importedCase()  {
	RandomGenerator *RNG = simStatus.getRandomGenerator();
	int classId = simStatus.getCurrentIndividualClass();
	int age = ageNames[ classId ];
	int xx  = simStatus.getX() % simStatus.getMaxX();
	int yy  = simStatus.getY();
	int ids_LA  = idsMap_LA[xx][yy];
//	int ids_NHS = idsMap_NHS[xx][yy];
	if (ids_LA == 0)   simStatus.abort("IDS_LA  is ZERO");
//	if (ids_NHS == 0)  simStatus.abort("IDS_NHS is ZERO");
//std::cout << ids << " " << age << " " << maxId_2019 << " " << "\n" << std::flush;
//	if (RNG->get() < importProbs_byLA[ids_LA][age] * importProbs_byNHS[ids_NHS][age])  {
	if (RNG->get() < importProbs_byLA[ids_LA][age])  {
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
	int ids = idsMap_LA[xx][yy];
	if (ids == 0)  simStatus.abort("IDS is ZERO");
	return  RNG->binomial( nn, importProbs_byLA[ids][age] );
}



// Mobility: duration of trips.  Here we assume simple commuting and daily work. Parameter will be fitted
inline  double getMobilityDuration(double dist)  {
	return  1;
}


inline  void  initCountrySpecific()  {
//	loadIdentifiers( identifiers_LA_File, identifiers_NHS_File );
	loadIdentifiers( identifiers_LA_File );
	loadImportProbs();
}


// FITTING  PROTOTYPE FOR GENERALIZATION
std::vector< int >  inputTable = {DATA_CUMUL_SYMPT, DATA_CUMUL_DEATHS};
//std::vector< int >  paramTable = {PARAM_T0, PARAM_R0, PARAM_OMEGA};
std::vector< int >  paramTable = {PARAM_T0, PARAM_R0, PARAM_GAMMA, PARAM_OMEGA};
std::vector< int >  distsTable = {DATA_DEATHS};
std::vector< int >  printTable = {DATA_DEATHS};
std::vector< int >  contrTable = {POLICY_SOCIALDIST_PROB, POLICY_SCHOOL_CLOSURE, POLICY_TRACING_PROB, POLICY_REDUCE_INFLIGHT};



bool  checkLockdown(int x0, int y0)  {
	return  true;
}


