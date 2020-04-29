// Relevant files loaded at execution time
static  std::string  contactMatrixFile = "../../../Italy/Contacts/ItalyContactMatrix";
static  std::string  ageGroupEsriFile  = "../../../Italy/Setup/Italy_%dkm_%d.dat";
static  std::string  identifiersFile   = "../../../Italy/Maps/Italy_%dkm_ids.asc";
static  std::string  timeseriesFile    = "../../../Italy/Italy_timeseries.dat";

// ITALY DATA
//
// 16th Feb first case (when the 1st detected appeared in hospital the first time)
// 24th Feb (+8d/+8d) Red Zones (Lodi province) [no work, school closed, no social, no mobility]
// The municipalities are around these coordinates: the area spans roughly 4kmx8km (I hope I am correct), which is 2 grid cells (50 km^2)
//		45,2180  9.6949
//		45.1537  9.6959
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
//  POLICY IDENTIFIERS
enum  {POLICY_24thFeb = 0, POLICY_01stMar, POLICY_04thMar, POLICY_09thMar, POLICY_LAST};
//	POLICY_SOCIALDIST_PROB = 0, 	// Generalized reduction of social interactions
//	POLICY_TRAVELREDUCTION, 		// Reduction of travel intensity
//	POLICY_STAYATHOME_AGE, 			// Compliance of stay at home for oldest (non-working)
//	POLICY_STAYATHOME_OTH, 			// Compliance of stay at home for working individuals
//	POLICY_STAYATHOME_SCH, 			// Compliance of stay at home for school-aged individuals
//	POLICY_FAMILY_TRANSMIT, 		// Increae in family transmission
//	POLICY_STAYATHOME_FULL, 		// Whether stay-at-home for younger and older means avoiding all social contacts (ex. no shopping at all)
//	POLICY_SCHOOL_CLOSURE, 			// Fraction of schools closed (generalized)

// Parameters for each of the above policies. Rows indicate distinct policy packages, columns indicate distinct policy actions
std::vector< std::array<double, 8> >  policyParams = {  // 4 packages, each with 8 params/actions (POLICY_TYPE_LAST == 8 and POLICY_LAST == 4).
	{0.80, 0.9, 0.8, 0.8, 0.8, 2.0, 0.0, 1.0},  // 24th Feb (lockdown local to some municipalities)
	{0.25, 0.0, 0.2, 0.2, 0.2, 2.0, 0.0, 1.0},  //  1st Mar (restriction local to yellow provinces)
	{0.20, 0.2, 0.2, 0.2, 0.2, 2.0, 0.0, 1.0},  //  4th Mar (nationwide school closures)
	{0.90, 0.0, 0.8, 0.8, 0.8, 2.0, 0.0, 1.0}   //  9th Mar (nationwide lockdown)
};
std::vector<double>  policyTime = {8, 14, 17, 22};  // Days of application from day zero.
std::vector< std::array<double, 8> >  policyDuration = {  // Duration of actions. Again, rows indicate distinct packages, columns distinct actions
	{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{5.0, 5.0, 5.0, 5.0, 5.0, 0.0, 0.0, 0.0},
	{7.0, 7.0, 7.0, 7.0, 7.0, 0.0, 0.0, 0.0}
};
std::vector<int>  policyApplication = {1, 2, 3, 3};  // Type and extent of policy: 1- initial outbreak areas(Lodi); 2- Regional (provinces); 3- Nationwide; 0- applied value

std::vector< std::vector<double> >  firstInfections(1, {0.0, 9.705, 45.16});


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

	int index = 3;
//		45,2180  9.6949
//		45.1537  9.6959
	// Presume to apply nationwide values;
	int xx = simStatus.getX() % simStatus.getMaxX();
	int yy = simStatus.getY();
	int provinceId = idsMap[xx][yy];
	if (std::find( yellow_provinces.begin(), yellow_provinces.end(), provinceId ) != yellow_provinces.end() && simStatus.getTime() < params.t0+policyTime[3])  {
		index = 2;
	}
	if ((simStatus.getX() - simStatus.getXfromLongitude(9.6954)) == 0 && simStatus.getTime() < params.t0+policyTime[3])  {
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

#ifndef  MODEL_FAMILY
		params.home     = FAMILY_TRANSMIT[0];
#endif
		params.work	= 1.0-STAYATHOME_OTH[0];
		params.school	= 1.0-SCHOOL_CLOSURE[0];
		params.other	= 1.0-SOCIALDIST_PROB[0];
		params.mobility = 1.0-TRAVELREDUCTION[0];
		updateContactMatrix();
		prev_index = index;
	}

}



// Mobility: duration of trips.  Here we assume simple commuting and daily work. Parameter will be fitted
inline  double getMobilityDuration(double rnd, int kk)  {
	return  1;
}


inline  void  initMobility()  {
	loadIdentifiers( identifiersFile );
}
