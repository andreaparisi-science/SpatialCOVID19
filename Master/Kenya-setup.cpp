// Relevant files loaded at execution time
static  std::string  contactMatrixFile = "../../Kenya/Contacts/KenyaContactMatrix";
static  std::string  ageGroupEsriFile  = "../../Kenya/Setup/Kenya_5km_%d.dat";
static  std::string  identifiersFile   = "../../Kenya/Maps/Kenya_5km_ids.asc";
static  std::string  timeseriesFile    = "../../Kenya/Kenya_timeseries.dat";

std::vector<int>  yellow_provinces = {16, 17, 13, 19, 75, 98, 20, 15, 108, 18, 14, 12, 25, 28, 29, 26, 27, 23, 24, 37, 34, 36, 35, 39, 99, 38, 40, 33};

enum  {POLICY_24thFeb = 0, POLICY_01stMar, POLICY_04thMar, POLICY_09thMar, POLICY_LAST};
enum  {	POLICY_SOCIALDIST_PROB = 0, 	// Generalized reduction of social interactions
		POLICY_TRAVELREDUCTION, 		// Reduction of travel intensity
		POLICY_STAYATHOME_AGE, 			// Compliance of stay at home for oldest (non-working)
		POLICY_STAYATHOME_OTH, 			// Compliance of stay at home for working individuals
		POLICY_STAYATHOME_SCH, 			// Compliance of stay at home for school-aged individuals
		POLICY_FAMILY_TRANSMIT, 		// Increae in family transmission
		POLICY_STAYATHOME_FULL, 		// Whether stay-at-home for younger and older means avoiding all social contacts (ex. no shopping at all)
		POLICY_SCHOOL_CLOSURE, 			// Fraction of schools closed (generalized)
		POLICY_TYPE_LAST
};
// Parameters for each of the above policies. Rows indicate distinct policy packages, columns indicate distinct policy actions
std::vector< std::array<double, 8> >  policyParams = {  // 4 policies, each with 8 params/actions (POLICY_TYPE_LAST == 8 and POLICY_LAST == 4).
	{0.80, 0.9, 0.8, 0.8, 0.8, 2.0, 0.0, 1.0},  // 24th Feb (lockdown local to some municipalities)
	{0.20, 0.0, 0.2, 0.2, 0.2, 2.0, 0.0, 1.0},  //  1st Mar (restriction local to yellow provinces)
	{0.20, 0.2, 0.2, 0.2, 0.2, 2.0, 0.0, 1.0},  //  4th Mar (nationwide school closures)
	{0.90, 0.0, 0.9, 0.9, 0.9, 2.0, 0.0, 1.0}   //  9th Mar (nationwide lockdown)
};
std::vector<double>  policyTime = {8, 14, 17, 22};  // Days of application from day zero.
std::vector< std::array<double, 8> >  policyDuration = {  // Duration of actions. Again, rows indicate distinct packages, columns distinct actions
	{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{7.0, 7.0, 7.0, 7.0, 7.0, 0.0, 0.0, 0.0},
	{10., 10., 10., 10., 0.0, 0.0, 10., 0.0}
};
std::vector<int>  policyApplication = {1, 2, 3, 3};  // Type and extent of policy: 1- initial outbreak areas(Lodi); 2- Regional (provinces); 3- Nationwide; 0- applied value

std::vector<double>  SOCIALDIST_PROB(3+1, 0.0); // [1..3] Current values at the different extent levels (above). [0] Current value in grid element
std::vector<double>  TRAVELREDUCTION(3+1, 0.0);
std::vector<double>  STAYATHOME_AGE(3+1, 0.0);
std::vector<double>  STAYATHOME_OTH(3+1, 0.0);
std::vector<double>  STAYATHOME_SCH(3+1, 0.0);
std::vector<double>  FAMILY_TRANSMIT(3+1, 1.0);
std::vector<double>  STAYATHOME_FULL(3+1, 0.0);
std::vector<double>  SCHOOL_CLOSURE(3+1, 0.0);

std::vector< std::vector<double> >  firstInfections(1, {0.0, 36.8172, -1.2864});


// R0 values for different transmission levels (calculated for the above values. They are *almost* independent on scaling of above values)
static std::vector<double>  R0_tau = {16.7546, 8.53905, 4.45365, 2.1001, 1.02682};

// Implements measures taking into account their spatial extent
void	evalLocalParameters()  {
	static int  prev_index = -1;
	static int  prev_day = -1;

	int day = static_cast<int>(simStatus.getTime()+0.0000001);
	if (day != prev_day)  {
		prev_day = day;
		prev_index = -1;
	}

	int index = 3;
//		45,2180  9.6949
//		45.1537  9.6959
	// Presume to apply nationwide values;
	int xx = simStatus.getX() % simStatus.getMaxX();
	int yy = simStatus.getY();
	int provinceId = idsMap[xx][yy];
	if (std::find( yellow_provinces.begin(), yellow_provinces.end(), provinceId ) != yellow_provinces.end())  {
		index = 2;
	}
	if ((simStatus.getX() - simStatus.getXfromLongitude(9.6954)) == 0)  {
		if (((simStatus.getY() - simStatus.getYfromLatitude(45.2180)) == 0) || ((simStatus.getY() - simStatus.getYfromLatitude(45.1537)) == 0))  {
			index = 1;
		}
	}
	// To be applied only if required, otherwise values are already correct
	if (prev_index != index)  {
		SOCIALDIST_PROB[0] = SOCIALDIST_PROB[index];
		TRAVELREDUCTION[0] = TRAVELREDUCTION[index];
		STAYATHOME_AGE[0]  = STAYATHOME_AGE[index];
		STAYATHOME_OTH[0]  = STAYATHOME_OTH[index];
		STAYATHOME_SCH[0]  = STAYATHOME_SCH[index];
		FAMILY_TRANSMIT[0] = FAMILY_TRANSMIT[index];
		STAYATHOME_FULL[0] = STAYATHOME_FULL[index];
		SCHOOL_CLOSURE[0]  = SCHOOL_CLOSURE[index];

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
	ifstream  handler("../../Kenya/Mobility/TripDistribution.tsv");
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


inline  double getMobilityDuration(double rnd, int kk)  {
	return  extractFromDistribution( rnd, durationDistr[kk], -1, 0 );
}


inline  void  initMobility()  {
	tripDistr = analyzeTripDistribution();
	durationDistr = makeDistr( 0.50, tripDistr );
}
