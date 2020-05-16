// Relevant files loaded at execution time
static  std::string  contactMatrixFile = "../../../Test/Contacts/KenyaContactMatrix";
static  std::string  ageGroupEsriFile  = "../../../Test/Setup/Test_%dkm_%d.dat";
static  std::string  identifiersFile   = "../../../Test/Maps/Kenya_%dkm_ids.asc";
static  std::string  timeseriesFile    = "../../../Test/Test_timeseries.dat";

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

// Parameters for each of the above policies. Rows indicate distinct policy packages, columns indicate distinct policy actions
std::vector< std::array<double, 10> >  policyParams = {  // 4 policies, each with 8 params/actions (POLICY_TYPE_LAST == 8 and POLICY_LAST == 4).
//	{1.00, 0.45, 0.10, 0.0, 0.0, 0.30,  0.0, 1.75, 0.0, 1.00},
//	{1.00, 0.45, 0.10, 0.9, 0.0, 0.30,  0.0, 1.75, 0.0, 1.00}
};
std::vector<double>  policyTime = {
//	8, 32
};    // Days of application from day zero.
std::vector< std::array<double, 10> >  policyDuration = {  // Duration of actions. Again, rows indicate distinct packages, columns distinct actions
//	{ 0, 14, 14, 14, 14, 14, 14, 14, 14,  8},
//	{ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0}
};
std::vector<int>  policyApplication = {
	1, 1
}; // Type and extent of policy: 1-Global; 0- applied value

std::vector< std::vector<double> >  firstInfections = {
	{0.0, 0.0, 36.8172, -1.2864},  // First absolute
	{0.0, 2.0, 36.8172, -1.2864},  // 5th of March (x3)
	{0.0, 2.0, 36.8172, -1.2864},  // 5th of March
	{0.0, 2.0, 36.8172, -1.2864},  // 5th of March
	{0.0, 2.0, 36.8172, -1.2864},  // 9th of March
	{0.0, 2.0, 36.8172, -1.2864},  // 17th of March
	{0.0, 2.0,  36.8172, -1.2864}, // 22nd of March (x8)
	{0.0, 2.0,  36.8172, -1.2864}, // 22nd of March
	{0.0, 2.0,  36.8172, -1.2864}, // 22nd of March
	{0.0, 2.0,  36.8172, -1.2864}, // 22nd of March
	{0.0, 2.0,  36.8172, -1.2864}, // 22nd of March
	{0.0, 2.0,  36.8172, -1.2864}, // 22nd of March
	{0.0, 2.0,  36.8172, -1.2864}, // 22nd of March
	{0.0, 2.0,  36.8172, -1.2864},  // 22nd of March
	{0.0, 0.0, 36.8172, -1.2864},  // First absolute
	{0.0, 2.0, 36.8172, -1.2864},  // 5th of March (x3)
	{0.0, 2.0, 36.8172, -1.2864},  // 5th of March
	{0.0, 2.0, 36.8172, -1.2864},  // 5th of March
	{0.0, 2.0, 36.8172, -1.2864},  // 9th of March
	{0.0, 2.0, 36.8172, -1.2864},  // 17th of March
	{0.0, 2.0,  36.8172, -1.2864}, // 22nd of March (x8)
	{0.0, 2.0,  36.8172, -1.2864}, // 22nd of March
	{0.0, 2.0,  36.8172, -1.2864}, // 22nd of March
	{0.0, 2.0,  36.8172, -1.2864}, // 22nd of March
	{0.0, 2.0,  36.8172, -1.2864}, // 22nd of March
	{0.0, 2.0,  36.8172, -1.2864}, // 22nd of March
	{0.0, 2.0,  36.8172, -1.2864}, // 22nd of March
	{0.0, 2.0,  36.8172, -1.2864}, // 22nd of March
	{0.0, 0.0, 36.8172, -1.2864},  // First absolute
	{0.0, 2.0, 36.8172, -1.2864},  // 5th of March (x3)
	{0.0, 2.0, 36.8172, -1.2864},  // 5th of March
	{0.0, 2.0, 36.8172, -1.2864},  // 5th of March
	{0.0, 2.0, 36.8172, -1.2864},  // 9th of March
	{0.0, 2.0, 36.8172, -1.2864},  // 17th of March
	{0.0, 2.0,  36.8172, -1.2864}, // 22nd of March (x8)
	{0.0, 2.0,  36.8172, -1.2864}, // 22nd of March
	{0.0, 2.0,  36.8172, -1.2864}, // 22nd of March
	{0.0, 2.0,  36.8172, -1.2864}, // 22nd of March
	{0.0, 2.0,  36.8172, -1.2864}, // 22nd of March
	{0.0, 2.0,  36.8172, -1.2864}, // 22nd of March
	{0.0, 2.0,  36.8172, -1.2864}, // 22nd of March
	{0.0, 2.0,  36.8172, -1.2864}, // 22nd of March
	{0.0, 2.0,  36.8172, -1.2864}, // 22nd of March (x8)
	{0.0, 2.0,  36.8172, -1.2864}, // 22nd of March
	{0.0, 2.0,  36.8172, -1.2864}, // 22nd of March
	{0.0, 2.0,  36.8172, -1.2864}, // 22nd of March
	{0.0, 2.0,  36.8172, -1.2864}, // 22nd of March
	{0.0, 2.0,  36.8172, -1.2864}, // 22nd of March
	{0.0, 2.0,  36.8172, -1.2864}, // 22nd of March
	{0.0, 2.0,  36.8172, -1.2864}  // 22nd of March
};
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
	int day = static_cast<int>(simStatus.getTime()+0.0000001);
	if (day != prev_day)  {
		prev_day = day;
		prev_index = -1;
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
		params.work		= 1.0-STAYATHOME_OTH[0];
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





void  initMobility()  {
	loadIdentifiers( identifiersFile );
}




std::vector<int>  lockdownCountyList = {
//	41, 47, 44, 46
}; // 41-Nairobi, 47-Mombasa, 44-Kilifi, 46-Kwale
bool  checkLockdown(int x0, int y0)  {
	RandomGenerator *RNG = simStatus.getRandomGenerator();
	if (simStatus.getTime() > policyTime[1])  {
		if (std::find(lockdownCountyList.begin(), lockdownCountyList.end(), idsMap[x0][y0]) != lockdownCountyList.end())  {
			return true;
		}
	}
	return false;
}



