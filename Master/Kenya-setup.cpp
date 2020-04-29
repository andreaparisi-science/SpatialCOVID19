// Relevant files loaded at execution time
static  std::string  contactMatrixFile = "../../../Kenya/Contacts/KenyaContactMatrix";
static  std::string  ageGroupEsriFile  = "../../../Kenya/Setup/Kenya_%dkm_%d.dat";
static  std::string  identifiersFile   = "../../../Kenya/Maps/Kenya_%dkm_ids.asc";
static  std::string  timeseriesFile    = "../../../Kenya/Kenya_timeseries.dat";

//	POLICY_SOCIALDIST_PROB = 0, 	// Generalized reduction of social interactions
//	POLICY_TRAVELREDUCTION, 		// Reduction of travel intensity
//	POLICY_STAYATHOME_AGE, 			// Compliance of stay at home for oldest (non-working)
//	POLICY_STAYATHOME_OTH, 			// Compliance of stay at home for working individuals
//	POLICY_STAYATHOME_SCH, 			// Compliance of stay at home for school-aged individuals
//	POLICY_FAMILY_TRANSMIT, 		// Increae in family transmission
//	POLICY_STAYATHOME_FULL, 		// Whether stay-at-home for younger and older means avoiding all social contacts (ex. no shopping at all)
//	POLICY_SCHOOL_CLOSURE, 			// Fraction of schools closed (generalized)

// Parameters for each of the above policies. Rows indicate distinct policy packages, columns indicate distinct policy actions
std::vector< std::array<double, 8> >  policyParams = {  // 4 policies, each with 8 params/actions (POLICY_TYPE_LAST == 8 and POLICY_LAST == 4).
};
std::vector<double>  policyTime = {    // Days of application from day zero.
};
std::vector< std::array<double, 8> >  policyDuration = {  // Duration of actions. Again, rows indicate distinct packages, columns distinct actions
};
std::vector<int>  policyApplication = {};  // Type and extent of policy: 1- initial outbreak areas(Lodi); 2- Regional (provinces); 3- Nationwide; 0- applied value

std::vector< std::vector<double> >  firstInfections(1, {0.0, 36.8172, -1.2864});


// Implements measures taking into account their spatial extent
void	evalLocalParameters()  {
	static int  prev_index = -1;
	static int  prev_day = -1;

	int day = static_cast<int>(simStatus.getTime()+0.0000001);
	if (day != prev_day)  {
		prev_day = day;
		prev_index = -1;
	}

	int index = 0;
	// To be applied only if required, otherwise values are already correct
	if (prev_index != index)  {
		if (index != 0)  {
			SOCIALDIST_PROB[0] = SOCIALDIST_PROB[index];
			TRAVELREDUCTION[0] = TRAVELREDUCTION[index];
			STAYATHOME_AGE[0]  = STAYATHOME_AGE[index];
			STAYATHOME_OTH[0]  = STAYATHOME_OTH[index];
			STAYATHOME_SCH[0]  = STAYATHOME_SCH[index];
			FAMILY_TRANSMIT[0] = FAMILY_TRANSMIT[index];
			STAYATHOME_FULL[0] = STAYATHOME_FULL[index];
			SCHOOL_CLOSURE[0]  = SCHOOL_CLOSURE[index];
		}

		params.work		= 1.0-STAYATHOME_OTH[0];
		params.school	= 1.0-SCHOOL_CLOSURE[0];
		params.other	= 1.0-SOCIALDIST_PROB[0];
		params.mobility = 1.0-TRAVELREDUCTION[0];
		updateContactMatrix();
		prev_index = index;
	}
}




static int MAXLEN = 210;
std::vector<double>  analyzeTripDistribution()  {
	std::vector<double> xx, yy;
	double fval, ymax, ymin, alpha, total = 0, weightot = 0;
	ifstream  handler("../../../Kenya/Maps/Mobility/TripDistribution.tsv");
	if (!handler.good())  {
		simStatus.abort( "Cannot find Trip distribution file." );
	}
	while (true)  {
		handler >> fval;
		if (handler.eof())  break;
		xx.push_back(fval);
		handler >> fval;
		yy.push_back(fval);
	}

	SplineInterpolator interp = SplineInterpolator(xx, yy, SPLINE_METHOD_NATURAL );
	xx.resize(MAXLEN+1);
	yy.clear();
	for (int kk = 0; kk <= MAXLEN; kk++)  {
		xx[kk] = kk;
	}
	yy = interp.interpolate(xx);
	// Rescaling to zero
	ymax = yy[0];
	ymin = yy[MAXLEN];
	alpha = ymax/(ymax-ymin);
	for (int kk = 0; kk <= MAXLEN; kk++)  {
		yy[kk] = alpha*(yy[kk]-ymin);
	}
	// Weighted and nonweighted sum
	for (int kk = 0; kk <= MAXLEN; kk++)  {
		total    += yy[kk];
		weightot += yy[kk]/(kk+1);
	}
	std::cout << "TRIP TOTALS: " << total << " " << weightot << "\n";
	return  yy;
}



std::vector< std::vector<double> >  makeDistr( double pp, std::vector<double>  tr )  {
	std::vector< std::vector<double> >  distr;
	distr.resize( 211 );
	for (int ii = 0; ii <= 210; ii++)  {
		distr[ii].resize( 211, 0.0 );
	}

	double  sum;
	for (int ii = 0; ii <= 210; ii++)  {
		sum = 0.0;
		for (int jj = 0; jj <= 210; jj++)  {
			distr[ii][jj] = pp*lse::Poisson(ii, jj) + (1-pp)*1.0/(211.0);
			sum += distr[ii][jj];
			distr[ii][jj] *= tr[jj];
		}
	}
	return distr;
}




// Mobility: duration of trips
static std::vector<double>  tripDistr;
static std::vector< std::vector<double> >  durationDistr;


inline  double getMobilityDuration(double rnd, double dist)  {
	RandomGenerator *RNG = simStatus.getRandomGenerator();
	rnd = RNG->get();
	int kk = static_cast<int>(210*dist/1000.0);
	if (kk > 210)  kk = 210;
	return  extractFromDistribution( rnd, durationDistr[kk], -1, 0 );
}


inline  void  initMobility()  {
	tripDistr = analyzeTripDistribution();
	durationDistr = makeDistr( 0.50, tripDistr );
}

