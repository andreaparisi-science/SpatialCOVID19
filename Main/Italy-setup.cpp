// Relevant files loaded at execution time
static  std::string  contactMatrixFile = "../../../../Data/Italy/Contacts/ItalyContactMatrix";
static  std::string  ageGroupEsriFile  = "../../../../Data/Italy/Setup/Italy_%dkm_%d.asc";
static  std::string  identifiersFile   = "../../../../Data/Italy/Maps/Italy_%dkm_ids.asc";
static  std::string  timeseriesFile    = "../../../../Data/Italy/Italy_timeseries.dat";
static  std::string  fileDeathsByAge   = "../../../../Data/Italy/Private/deathsByAge.dat";

static  std::string  importedCasesFile = "../../../../Data/Italy/Private/distrCases_byId.dat";
// ITALY DATA
//
// 16th Feb first case (when the 1st detected appeared in hospital the first time)
// 22nd Feb (+8d/+8d) Red Zones (Lodi province) [no work, school closed, no social, no mobility]
// The municipalities are around these coordinates: the area spans roughly 4kmx8km (I hope I am correct) for 50k inhabitants,
// and covers a good part of the Lodi province.  It is roughly 2 grid cells (50 km^2)
//		45,2180  9.6949
//		45.1537  9.6959
// Vo' Euganeo: this is a single town of 3000 inhabitants.  We do not focus on municipalities thus we leave this out.
//  1st Mar (+6d/+14d)   Red Zones + Yellow zones = Lombardy, Veneto, Emilia-Romagna [school closed social reduced]
//      Red zones are the same as above, the yellow zones correspond to the following provinces
//				Lombardy:	16	Bergamo
//							17	Brescia
//							13	Como
//							19	Cremona
//							75	Lecce
//							98	Lodi
//							20	Mantova
//							15	Milano
//							108	Monza e Brianza
//							18	Pavia
//							14	Sondrio
//							12	Varese
//				Veneto:		25	Belluno
//							28	Padova
//							29	Rovigo
//							26	Treviso
//							27	Venezia
//							23	Verona
//							24	Vicenza
//				Emilia-Rom:	37	Bologna
//							34	Parma
//							36	Modena
//							35	Reggio Emilia
//							39	Ravenna
//							99	Rimini
//							38	Ferrara
//							40	Forli-Cesena
//							33	Piacenza
std::vector<int>  yellow_provinces = {16, 17, 13, 19, 75, 98, 20, 15, 108, 18, 14, 12, 25, 28, 29, 26, 27, 23, 24, 37, 34, 36, 35, 39, 99, 38, 40, 33};
//
//  4th Mar (+3d/+17d)   Nationwide school closures
//  9th Mar (+5d/+22d)   Nationwide shutdown 3 weeks 90%. This occurred in 2 phases, but the gap is just 1 day, 
//                       so we approximate this with a unique modification.  on the 8th the North was in lockdown and the centre-south was yellow area,
//                       On the 9th (executive 9th or 10th?) the full of Italy was in lockdown. In reality all of Italy was effectively entering
//                       lockdown due to people not moving around as before.  Industrial activities were blocked on the 22nd of Mar, but this is taken
//                       into account by the increased reduction from 9th to 22nd (2 weeks). Data shows that reduction effectively was already achieved
//						 in 10 days (Google)
//
Intervention  nationwideLockdown;

enum  {EXTENT_LOCAL = 1, EXTENT_PROVINCE, EXTENT_NATIONAL};
std::vector<Intervention>	interventions = {
//	Intervention( POLICY_REDUCE_INFLIGHT, EXTENT_NATIONAL,  0,  0, 1.00 )
	Intervention( POLICY_TRACING_PROB,    EXTENT_LOCAL,     0,  0, 1.00 ),  // 23rdFeb (lockdown local to some municipalities)
	Intervention( POLICY_TRACING_PROB,    EXTENT_PROVINCE,  0,  0, 1.00 ),  // 23rdFeb (lockdown local to some municipalities)
	Intervention( POLICY_TRACING_PROB,    EXTENT_NATIONAL,  0,  0, 1.00 ),  // 23rdFeb (lockdown local to some municipalities)
	Intervention( POLICY_SOCIALDIST_PROB, EXTENT_LOCAL,     7,  0, 0.80 ),
	Intervention( POLICY_TRAVELREDUCTION, EXTENT_LOCAL,     7,  0, 0.90 ),
	Intervention( POLICY_STAYATHOME_AGE,  EXTENT_LOCAL,     7,  0, 0.80 ),
	Intervention( POLICY_STAYATHOME_OTH,  EXTENT_LOCAL,     7,  0, 0.80 ),
	Intervention( POLICY_STAYATHOME_SCH,  EXTENT_LOCAL,     7,  0, 0.80 ),
	Intervention( POLICY_STAYATHOME_FULL, EXTENT_LOCAL,	7,  0, 1.00 ),
	Intervention( POLICY_FAMILY_TRANSMIT, EXTENT_LOCAL,     7,  0, 2.00 ),
	Intervention( POLICY_SCHOOL_CLOSURE,  EXTENT_LOCAL,     7,  0, 1.00 ),

	Intervention( POLICY_REDUCE_INFLIGHT, EXTENT_LOCAL,     8, 12, 1.00 ),
	Intervention( POLICY_REDUCE_INFLIGHT, EXTENT_PROVINCE,  8, 12, 1.00 ),
	Intervention( POLICY_REDUCE_INFLIGHT, EXTENT_NATIONAL,  8, 12, 1.00 ),
	Intervention( POLICY_SOCIALDIST_PROB, EXTENT_PROVINCE,  8,  0, 0.15 ),  //  1st Mar (restriction local to yellow provinces)
	Intervention( POLICY_STAYATHOME_OTH,  EXTENT_PROVINCE,  8,  0, 0.10 ),  // Work reduction
	Intervention( POLICY_STAYATHOME_AGE,  EXTENT_PROVINCE,  8,  0, 0.15 ),
	Intervention( POLICY_TRAVELREDUCTION, EXTENT_PROVINCE,  8,  0, 0.15/0.90 ),

	Intervention( POLICY_TRACING_PROB,    EXTENT_LOCAL,    11,  7, 0.00 ),  // 27th Feb testing on symptomatics only
	Intervention( POLICY_TRACING_PROB,    EXTENT_PROVINCE, 11,  7, 0.00 ),  // 27th Feb testing on symptomatics only
	Intervention( POLICY_TRACING_PROB,    EXTENT_NATIONAL, 11,  7, 0.00 ),  // 27th Feb testing on symptomatics only

	Intervention( POLICY_SOCIALDIST_PROB, EXTENT_NATIONAL, 14, 14, 0.85 ),  //  1st Mar (restriction local to yellow provinces)
	Intervention( POLICY_STAYATHOME_AGE,  EXTENT_NATIONAL, 14, 14, 0.85 ),
	Intervention( POLICY_STAYATHOME_SCH,  EXTENT_NATIONAL, 14, 14, 0.85 ),
	Intervention( POLICY_STAYATHOME_OTH,  EXTENT_NATIONAL, 14, 14, 0.65 ),
	Intervention( POLICY_FAMILY_TRANSMIT, EXTENT_PROVINCE, 14,  0, 2.00 ),
	Intervention( POLICY_SCHOOL_CLOSURE,  EXTENT_PROVINCE, 14,  0, 1.00 ),

	Intervention( POLICY_SCHOOL_CLOSURE,  EXTENT_NATIONAL, 17,  0, 1.00 ),  //  4th Mar (nationwide school closures)
	Intervention( POLICY_FAMILY_TRANSMIT, EXTENT_NATIONAL, 17,  0, 2.00 ),

		nationwideLockdown = // assignment within assignment, for specific referencing elsewhere
	Intervention( POLICY_TRAVELRED_ADMIN, EXTENT_NATIONAL, 22,  0, 0.00 ),  // 95% 
	Intervention( POLICY_STAYATHOME_FULL, EXTENT_NATIONAL, 22,  0, 1.00 )
/*
//	Intervention( POLICY_SOCIALDIST_PROB, EXTENT_NATIONAL, 44,  7, 0.85 ),  //  Increases in social distancing and isolation follows data from google trends

//''	Intervention( POLICY_TRACING_PROB,    EXTENT_NATIONAL, 52, 10, 1.00 )  //  ~7 April - 10 days increses in testing with new test-isolate policy
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



// Time; Type (1.0 -> Exp, 2.0 -> Inf, 3.0 -> Asy; Lon; Lat;
std::vector< std::vector<double> >  firstInfections; //(1, {0.0, 0.0, 0.0, 9.705, 45.16});

std::vector< std::vector<int> >    idsMap;
std::vector< std::vector<double> > ageIdsPop;
std::vector<double>  idsPop;
int   maxId = 0;
void  loadIdentifiers( const std::string datafile )  {
	char ch_filename[256];
	int ids, nAgeGroups = groups[ groups.size()-1 ] + 1;
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
			ids = idsMap[ii][jj];
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
//		45,2180  9.6949
//		45.1537  9.6959
	// Presume to apply nationwide values;
	int xx = simStatus.getX() % simStatus.getMaxX();
	int yy = simStatus.getY();
	int provinceId = idsMap[xx][yy];
	if (std::find( yellow_provinces.begin(), yellow_provinces.end(), provinceId ) != yellow_provinces.end() && simStatus.getTime() < params.t0+nationwideLockdown.getActivationTime())  {
		index = EXTENT_PROVINCE;
	}
	if ((simStatus.getX() - simStatus.getXfromLongitude(9.6954)) == 0 && simStatus.getTime() < params.t0+nationwideLockdown.getActivationTime())  {
		if (((simStatus.getY() - simStatus.getYfromLatitude(45.2180)) == 0) || ((simStatus.getY() - simStatus.getYfromLatitude(45.1537)) == 0))  {
			index = EXTENT_LOCAL;
		}
	}
	return  index;
}



std::vector< std::vector<double> >  importProbs;  
void  loadImportProbs()  {
	std::ifstream  instream( importedCasesFile );
	int nAgeGroups = groups[ groups.size()-1 ] + 1;
	int  ids, age;
	double dummy, prob;
	DataBuffer  shareInit;	

	importProbs.resize( maxId+1 );
	for (int jj = 1; jj <= maxId; jj++)  {
		importProbs[jj].resize( nAgeGroups, 0.0 );
	}
	params.iomega = 0.0;
	while (instream.good())  {
		instream >> ids;
		if (instream.eof())  break;
		instream >> prob;

//std::cout << "***** " << maxId << " " << ids << " "  << age << "\n" << std::flush;
		for (int age = 0; age < nAgeGroups; age++)  {
//			importProbs[ids][age] = prob / ageIdsPop[age][ids];
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
std::vector< int >  inputTable = {DATA_SYMPT, DATA_DEATHS};
//std::vector< int >  paramTable = {PARAM_T0, PARAM_R0, PARAM_OMEGA};
std::vector< int >  paramTable = {PARAM_T0, PARAM_R0, PARAM_GAMMA, PARAM_OMEGA};
//std::vector< int >  paramTable = {PARAM_R0};
std::vector< int >  distsTable = {DATA_DEATHS};
std::vector< int >  printTable = {DATA_DEATHS};
std::vector< int >  contrTable = {POLICY_SOCIALDIST_PROB, POLICY_SCHOOL_CLOSURE, POLICY_TRACING_PROB, POLICY_REDUCE_INFLIGHT};



enum  {REG_UNDEF = 0, REG_ABRUZZO, REG_VALLE_AOSTA, REG_PUGLIA, REG_BASILICATA,
	REG_CALABRIA, REG_CAMPANIA, REG_EMILIA_ROMAGNA, REG_FRIULI_VENEZIA_GIULIA,
	REG_LAZIO, REG_LIGURIA, REG_LOMBARDIA, REG_MARCHE, REG_MOLISE,
	REG_PIEMONTE, REG_SARDEGNA, REG_SICILIA, REG_TRENTINO, REG_TOSCANA,
	REG_UMBRIA, REG_VENETO};

int  provinceToRegion( int province )  {
	int  region = REG_UNDEF;

	switch (province)  {
		case 12: case 13: case 14: case 15: case 16: case 17: case 18: case 19: case 20: case 97: case 98: case 108:
			region = REG_LOMBARDIA;
			break;
		case 56: case 57: case 58: case 59: case 60:
			region = REG_LAZIO;
			break;
		case 61: case 62: case 63: case 64: case 65:
			region = REG_CAMPANIA;
			break;
		case 81: case 82: case 83: case 84: case 85: case 86: case 87: case 88: case 89:  
			region = REG_SICILIA;
			break;
		case 23: case 24: case 25: case 26: case 27: case 28: case 29: 
			region = REG_VENETO;
			break;
		case 33: case 34: case 35: case 36: case 37: case 38: case 39: case 40: case 99:
			region = REG_EMILIA_ROMAGNA;
			break;
		case  1: case  2: case  3: case  4: case  5: case  6: case 96: case 103:
			region = REG_PIEMONTE;
			break;
		case 71: case 72: case 73: case 74: case 75: case 110:
			region = REG_PUGLIA;
			break;
		case 45: case 46: case 47: case 48: case 49: case 50: case 51: case 52: case 53: case 100:
			region = REG_TOSCANA;
			break;
		case 78: case 79: case 80: case 101: case 102: 
			region = REG_CALABRIA;
			break;
		case 90: case 91: case 92: case 95: case 104: case 105: case 106: case 107:
			region = REG_SARDEGNA;
			break;
		case  8: case  9: case 10: case 11:
			region = REG_LIGURIA;
			break;
		case 41: case 42: case 43: case 44: case 109:
			region = REG_MARCHE;
			break;
		case 66: case 67: case 68: case  69:
			region = REG_ABRUZZO;
			break;
		case 30: case 31: case 32: case 93: 
			region = REG_FRIULI_VENEZIA_GIULIA;
			break;
		case 21: case 22: 
			region = REG_TRENTINO;
			break;
		case 54: case 55:
			region = REG_UMBRIA;
			break;
		case 76: case 77: 
			region = REG_BASILICATA;
			break;
		case 70: case 94: 
			region = REG_MOLISE;
			break;
		case 7:
			region = REG_VALLE_AOSTA;
			break;
		default:
			simStatus.abort( "Unrecognized province.");
	}
	return  region;
}



bool  checkLockdown(int x0, int y0)  {
	return  true;
}


