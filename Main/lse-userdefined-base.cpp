#include <fstream>
#include <assert.h>
#include <queue>
#include <algorithm>
#include <set>
#include <sys/stat.h>  // mkdir (Linux)
#include <memory>  // smart pointers
#include "armadillo"

#include "SplineInterpolator.h"
#include "Policy.h"
#include "Intervention.h"

#ifdef FITTING
#warning( "Macro FITTING was defined." )
#endif

#ifdef MODEL_FAMILY
#warning( "Macro MODEL_FAMILY was defined." )
#endif

#ifdef MODEL_SEIR
#warning( "Macro MODEL_SEIR was defined." )
#endif

#ifdef TAULEAP
#warning( "Macro TAULEAP was defined." )
#endif


// Levels of transmissing levels for asymptomatics as for 
// Brand S. etal, 2020.
enum  {DET_0_00, DET_0_10, DET_0_25, DET_0_50, DET_1_00};

// Used in fitting procedure
static constexpr int TT = 1;   // Number of post-convergence generations

// If this is NOT a fitting run
int    constexpr DET_RATE = DET_1_00;	// Transmission rate for asymptomatics
double constexpr SYMPTOMATIC_MULTIPLIER = 0.9; // Transmission rate for oldest class

// If attack rate is close to natural
static double constexpr  BETA_FAMILY_MULTIPLIER = 1.0;
static double constexpr  R0_FAMILY_MULTIPLIER   = 1.17212;


int constexpr ISOLATE_GAPDAYS = 3;   // Gap days from probable case detection to isolation of traced individuals
               
#ifdef  TAULEAP
static bool tracingOn = false;   // Flags if tracing is being carried on (this should be always on, except for fast runs with no tracing)
#else
static bool tracingOn = true;   // Flags if tracing is being carried on (this should be always on, except for fast runs with no tracing)
#endif
static int  totalCases = 0;     // Process cumul number of cases
static int  gtotalCases = 0;    // Global cumul number of cases (for output)
static int  totalAsyCases = 0;  // Process cumul number of asymptomatic cases
static int  gtotalAsyCases = 0; // Global cumul number of asymptomatic cases (for output)
static int  mainUniqueId = 0;   // Unique Id for infectious (for contact tracing). It will get a unique starting value per process
static int  TOTCLASSES = 0; // Number of classes (to be set by ageNames())

static arma::mat  chineseKK;  // Stores the Chinese matrix for updating eigenvalues based on suscepitbility and transmissibility

// Vector storing policies (to be filled in country specific c++ setup files)
std::vector<double>  TRACING_PROB;
std::vector<double>  SOCIALDIST_PROB;
std::vector<double>  TRAVELREDUCTION;
std::vector<double>  TRAVELRED_ADMIN;
std::vector<double>  STAYATHOME_AGE;
std::vector<double>  STAYATHOME_OTH;
std::vector<double>  STAYATHOME_SCH;
std::vector<double>  FAMILY_TRANSMIT;
std::vector<double>  STAYATHOME_FULL;
std::vector<double>  SCHOOL_CLOSURE;
std::vector<double>  REDUCE_INFLIGHT;

static  std::set<int>  contactEvents;
static  std::map<int, int>  ageNames;  // from class to age group
static  std::map<int, int>  infectNames;
static  std::map<int, int>  susceptNames;
static  std::array<int,3>   occupancy;
static  std::array<int, 3>  goccupancy;
static  std::vector<int>    cumulDeaths;
static  std::vector<int>    gcumulDeaths;
static  std::vector<int>    cumulCases;
static  std::vector<int>    gcumulCases;
static  std::vector<double> gdistrDeaths;
static  double gdistrDeathsSum;
enum {OCC_HOS = 0, OCC_ICU, OCC_HOM};


// Age groups for: 1. (ages) tiff files from WorldPop database; 2. (groups) coarsed age groups to be used in this sim
std::vector<int> ages   = {0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80};
//std::vector<int> groups = {0, 0, 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16};
std::vector<int> groups = {0, 0, 0,  1,  1,  2,  2,  3,  3,  4,  4, 5, 5, 6, 6, 7, 7, 8};

std::array<int,2>  stayAtHomeLimit;

// Fitting macros
enum {PARAM_T0 = 1, PARAM_R0, PARAM_GAMMA, PARAM_TRACING, PARAM_TAU, PARAM_ZMAX, PARAM_OMEGA, PARAM_BETA, PARAM_SUSCEPTIBILITY, PARAM_SUSCEPTIBILITY_FULL, 
					PARAM_TRANSMISSIBILITY, PARAM_TRANSMISSIBILITY_FULL};
enum {DATA_ALL_CASES, DATA_SYMPT, DATA_ASYMPT, DATA_DEATHS, DATA_DUMMY, DATA_CUMUL_ALL_CASES, DATA_CUMUL_SYMPT, DATA_CUMUL_ASYMPT, DATA_CUMUL_DEATHS};

std::vector< std::vector< std::vector<double> > > baseMap;
std::vector< std::vector<double> > popMap;
std::vector<double>  sizes;

int firstinfected;
//bool  firstimported;
//bool  globalFirstimported;
//int   countDaysDelay;
class Infos;

// Prototype of functions
void  enforceFamilyTransmission(int, int, Infos*);
void  stayAtHome();
void  buildAgeAssignments();
void  familyTransmission(int, int, Infos*);
void  familySetup(int, int);
void  updateContactMatrix();
void  initContactMatrix();
void  getSymptomaticRate(int, double);
void  getSusceptibility();
void  enforceIsolation(int, int, Infos*);
void  shareContacts();
void  recovery(Infos*);
void  enforcePolicies();
void  shareOutputData();
void  doFitting(int, PolicyQueue&);
void  initChineseContactMatrix();
void  updateChineseEigenvalue();

#define  Test  -1
#define  Italy  1
#define  Kenya  2
#define  Spain  3
#define  Germany 4
#define  Switzerland 5
#define  China  6
#define  UK 7

#ifndef  COUNTRY
#error   COUNTRY is not defined.  Use '--define COUNTRY=XXX'  where  XXX is an implemented country.
#error   Currently implemented countries are: KENYA and ITALY.
#else
#	if COUNTRY == Italy
#		include "Italy-setup.cpp"
#	elif COUNTRY == Kenya
#		include "Kenya-setup.cpp"
#	elif COUNTRY == Test
#		include "Test-setup.cpp"
#	elif COUNTRY == UK
#		include "UK-setup.cpp"
#	elif COUNTRY == Spain
#		include "Spain-setup.cpp"
#	else
#		error Unrecognized name for COUNTRY.
#	endif
#endif



// This is the class storing individual data
class Infos : public IndivData {
public:
	Infos() { init(); }
	~Infos() {}
	void  init()  {
		isolate = 0;
		isolate_duration = 0;
		athome = 0;
		infamily = 0;
		length = 0;
		uniqueId = 0;
		athome_prev = 0.0;
		contacts.clear();
		contacts.shrink_to_fit();
		family.clear();
		family.shrink_to_fit();
	}
	void  packData( unsigned int& len, char* buffer )  {
		int sz, pos = 0;
		memcpy(buffer+pos, &isolate, sizeof(char)); pos += sizeof(char);
		memcpy(buffer+pos, &isolate_duration, sizeof(char)); pos += sizeof(char);
		memcpy(buffer+pos, &athome, sizeof(char)); pos += sizeof(char);
		memcpy(buffer+pos, &infamily, sizeof(char)); pos += sizeof(char);
		memcpy(buffer+pos, &length,  sizeof(int)); pos += sizeof(int);
		memcpy(buffer+pos, &uniqueId, sizeof(int)); pos += sizeof(int);
		memcpy(buffer+pos, &athome_prev, sizeof(double)); pos += sizeof(double);
		sz = contacts.size();
		memcpy(buffer+pos, &sz, sizeof(int)); pos += sizeof(int);
		for (int jj = 0; jj < sz; jj++)  {
			memcpy(buffer+pos, &contacts[jj], sizeof(int)); pos += sizeof(int);
		}
		sz = family.size();
		memcpy(buffer+pos, &sz, sizeof(int)); pos += sizeof(int);
		for (int jj = 0; jj < sz; jj++)  {
			memcpy(buffer+pos, &family[jj], sizeof(int)); pos += sizeof(int);
		}
		len = pos;
	}
	void  unpackData( unsigned int len, char* buffer )  {
		int sz, pos = 0;
		memcpy(&isolate, buffer+pos, sizeof(char)); pos += sizeof(char);
		memcpy(&isolate_duration, buffer+pos, sizeof(char)); pos += sizeof(char);
		memcpy(&athome, buffer+pos, sizeof(char)); pos += sizeof(char);
		memcpy(&infamily, buffer+pos, sizeof(char)); pos += sizeof(char);
		memcpy(&length, buffer+pos, sizeof(int)); pos += sizeof(int);
		memcpy(&uniqueId, buffer+pos, sizeof(int)); pos += sizeof(int);
		memcpy(&athome_prev, buffer+pos, sizeof(double)); pos += sizeof(double);
		memcpy(&sz, buffer+pos, sizeof(int)); pos += sizeof(int);
		contacts.assign(sz, 0);
		for (int jj = 0; jj < sz; jj++)  {
			memcpy(&contacts[jj], buffer+pos, sizeof(int)); pos += sizeof(int);
		}
		memcpy(&sz, buffer+pos, sizeof(int)); pos += sizeof(int);
		family.assign(sz, 0);
		for (int jj = 0; jj < sz; jj++)  {
			memcpy(&family[jj], buffer+pos, sizeof(int)); pos += sizeof(int);
		}
		assert(len == pos);
	}
	char isolate;
	char isolate_duration;
	char athome;
	char infamily;
	int length;
	int uniqueId;
	double  athome_prev;
	std::vector<int> contacts;
	std::vector<int> family;
private:
};


#ifndef  TAULEAP
Infos * indivDataFactory()  {
	Infos *data = nullptr;
	try {
		data = new Infos();
	} catch (std::bad_alloc &err)  {
		simStatus.abort("Bad Infos() allocation.");
	}
	return data;
}
void  indivDataDelete(IndivData *data)  {
	Infos* infos = reinterpret_cast<Infos*>( data );
	delete infos;
}
#endif



// This function reads in data from the age stratified maps,
// builds an additional global map (summing all age inhabitant
// per grid cell, this should be equivalent to the non age specific
// map generated by the procedures in the Setup directory).
//
#define  xstrstr(a) #a
#define  xstr(a)  xstrstr(a)

void  getAgeGroups()  {
	std::ifstream handler;
	char ch_filename[256];
	std::string  filename;
	int ncols;
	int nrows;
	std::string dummy;

	int nAgeGroups = groups[ groups.size()-1 ] + 1;
	cumulDeaths.resize(nAgeGroups);
	gcumulDeaths.resize(nAgeGroups);
	cumulCases.resize(nAgeGroups);
	gcumulCases.resize(nAgeGroups);

	ncols = simStatus.getMaxX();
	nrows = simStatus.getMaxY();
	std::cout << "Getting Age Group data... ";
	baseMap.resize( nAgeGroups );
	for (int aa = 0; aa < nAgeGroups; aa++)  {
		baseMap[aa].resize( ncols );
		for (int qq = 0; qq < ncols; qq++)  {
			baseMap[aa][qq].resize( nrows, 0 );
		}
	}

	popMap.resize( ncols );
	for (int kk = 0; kk < ncols; kk++)  {
		popMap[kk].resize( nrows, 0 );
	}

	for (int aa = 0; aa < nAgeGroups; aa++)  {
		sprintf( ch_filename, ageGroupEsriFile.c_str(), GRIDRES, aa );
		filename = std::string( ch_filename );
		std::cout << "Handling file [" << filename << "]\n" << std::flush;
		handler.open( filename );
		if (!handler.good())  {
			simStatus.abort( "File [" + filename + "] not found." );
		}
		// Discard header infos
		for (int qq = 0; qq < 12; qq++)  handler >> dummy;


		// Read contents
		for (int qq = nrows-1; qq >= 0; qq--)  {
			for (int kk = 0; kk < ncols; kk++)  {
				handler >> baseMap[aa][kk][qq];
			}
		}
		handler.close();
	}

	// Builds country level age group sizes
	// (not used elsewhere for the moment)
	// The same information is now generated by the
	// Setup procedures
	sizes.resize( nAgeGroups, 0 );
	for (int aa = 0; aa < nAgeGroups; aa++)  {
		for (int kk = 0; kk < ncols; kk++)  {
			for (int qq = 0; qq < nrows; qq++)  {
				popMap[kk][qq] += baseMap[aa][kk][qq];
				sizes[aa] += baseMap[aa][kk][qq];
			}
		}
	}

	std::cout << "done.\n";
}



void  resetCounters()  {
	totalCases = 0;     // Process cumul number of cases
	gtotalCases = 0;    // Global cumul number of cases (for output)
	totalAsyCases = 0;  // Process cumul number of asymptomatic cases
	gtotalAsyCases = 0; // Global cumul number of asymptomatic cases (for output)
	occupancy = goccupancy = {0, 0, 0};
	for (int jj = 0; jj < cumulDeaths.size(); jj++)  cumulDeaths[jj] = gcumulDeaths[jj] = 0;
	for (int jj = 0; jj < cumulCases.size(); jj++)   cumulCases[jj]  = gcumulCases[jj] = 0;
	firstinfected = 0;
}



void  extractFirstInfective()  {
	int  kk, maxKK, selKK, xx, yy, xstart, ystart, ids, age;
	int nAgeGroups = groups[ groups.size()-1 ] + 1;
	RandomGenerator *deryaRNG = simStatus.getRandomGenerator();

	// First: build cumulative by ids (and age)
	kk = 0;
	std::vector<double>  cumul((1+maxId)*nAgeGroups, 0);
	for (ids = 1; ids <= maxId; ids++)  {
		for (age = 0; age < nAgeGroups; age++)  {
			cumul[kk] = importProbs[ids][age];
			if (kk > 0)  cumul[kk] += cumul[kk-1];
			kk++;
		}
	}
	maxKK = kk;

	// Second: extract ids (and age)
	selKK = extractFromDistribution( deryaRNG->get(), cumul, maxKK, 1 );
	kk = 0;
	for (ids = 1; ids <= maxId; ids++)  {
		for (age = 0; age < nAgeGroups; age++)  {
			if (kk == selKK)  {
				goto extractedIds;
			}
			kk++;
		}
	}
	extractedIds:
		; // Nothing to do

	// Third: build cumulative of location|ids,age
	kk = 0;
	cumul.clear();
	cumul.resize( simStatus.getMaxX()*simStatus.getMaxY(), 0 );
	for (xx = 0; xx < simStatus.getMaxX(); xx++)  {
		for (yy = 0; yy < simStatus.getMaxY(); yy++)  {
			if (idsMap[xx][yy] == ids)  {
				cumul[kk] = baseMap[ age ][xx][yy];
				if (kk > 0)  cumul[kk] += cumul[kk-1];
				kk++;
			}
		}
	}
	maxKK = kk;

	// Fourth: extract location
	selKK = extractFromDistribution( deryaRNG->get(), cumul, maxKK, 1 );
	kk = 0;
	for (xx = 0; xx < simStatus.getMaxX(); xx++)  {
		for (yy = 0; yy < simStatus.getMaxY(); yy++)  {
			if (idsMap[xx][yy] == ids)  {
				if (kk == selKK)  {
					goto extractedLocation;
				}
				kk++;
			}
		}
	}
	extractedLocation:
		; // Nothing to do

	// Fifth: share among nodes
	int izero = 0;
	DataBuffer  shareInit;
	shareInit.setBuffer( sizeof(int)*3, 3 );
	shareInit.clear();
	if (simStatus.getProcessId() == 0)  {
		shareInit.pack( &xx, DATA_INT, DATA_ADD );
		shareInit.pack( &yy, DATA_INT, DATA_ADD );
		shareInit.pack( &age, DATA_INT, DATA_ADD );
	} else {
		shareInit.pack( &izero, DATA_INT, DATA_ADD );
		shareInit.pack( &izero, DATA_INT, DATA_ADD );
		shareInit.pack( &izero, DATA_INT, DATA_ADD );
	}
	shareInit.reduce();
	shareInit.unpack( &xstart );
	shareInit.unpack( &ystart );
	shareInit.unpack( &age );
//std::cout << "FIRST " << xstart << " " << ystart << " " << age << "\n" << std::flush;
	firstInfections.clear();
//std::cout << "CH 0\n" << std::flush;
	firstInfections.push_back( std::vector<double>( {0.0, 0.0, static_cast<double>(age), simStatus.getLongitude(xstart), simStatus.getLatitude(ystart)} ) );
//std::cout << "CH 1\n" << std::flush;
}




// Initialize simulation. Executed per grid cell, per individual
void  init()  {
	int nAgeGroups = groups[ groups.size()-1 ] + 1;
//	static std::vector<double>  probs;
	static std::vector<int>  occupy;
	static int xpre = -1, ypre = -1;
	int ncols = simStatus.getMaxX();
	int nrows = simStatus.getMaxY();
	int indivId = simStatus.getCurrentIndividualId();
	int classId = simStatus.getCurrentIndividualClass();
	int xx = simStatus.getX() % ncols;
	int yy = simStatus.getY();
	int nn = 0, count = 0;
	Infos *data;

#ifndef  TAULEAP
	data = reinterpret_cast<Infos*>( simStatus.getIndividualData( indivId, classId ) );
	if (data == nullptr)  {
		data = indivDataFactory();
		simStatus.setIndividualData( indivId, classId, data );
	}
	data->init();
#endif

	// If this is a new grid cell, build age distribution
	if (xx != xpre || yy != ypre)  {
		occupy.clear();
		occupy.resize( nAgeGroups, 0 );
		for (int jj = 0; jj < nAgeGroups; jj++)  {
			occupy[jj]  = baseMap[jj][xx][yy];
		}
		assert(popMap[xx][yy] > 0);
//		assert( probs[nAgeGroups-1] > 0.999998 && probs[nAgeGroups-1] <= 1.000002 );
		xpre = xx; ypre = yy;
		nn = 0;
		count = 0;
	}

	// Move individuals from the base class (Sus_0) to class required by age distribution
	while (nn < nAgeGroups)  {
		if (occupy[nn] > 0)  {
			simStatus.moveCurrentIndividualToClass( classmapper["Sus_" + std::to_string(nn)] );
			occupy[nn]--;
			count++;
			break;
		} else {
			nn++;
		}
	}
	if (nn == nAgeGroups)  {
		std::cout << "Problme is here: " << count << " " << simStatus.getLocalPopulationSize() << "\n" << std::flush;
		assert( count == simStatus.getLocalPopulationSize() );
	}
}



#ifndef  TAULEAP
bool  infect(int indivId, int classId, Infos* data)  {
	static std::vector<int>  infClasses;
	static bool first = true;
	int nAgeGroups = groups[ groups.size()-1 ] + 1;

	if (first)  {
		for (int jj = 0; jj < nAgeGroups; jj++)  {
			for (int kk = 0; kk < classcounter["Inf_"+std::to_string(jj)]; kk++)  {
				infClasses.push_back( classmapper["Inf_"+std::to_string(jj)] + kk );
			}
		}
		for (int jj = 0; jj < nAgeGroups; jj++)  {
			for (int kk = 0; kk < classcounter["Asy_"+std::to_string(jj)]; kk++)  {
				infClasses.push_back( classmapper["Asy_"+std::to_string(jj)] + kk );
			}
		}
		first = false;
	}

	//if (data->athome == 1)  return false;
	if (data->isolate == 1)  return false;

	int sz = infClasses.size();
	std::vector<double> pp(sz, 0.0);

	for (int kk = 0; kk < sz; kk++)  {
		pp[kk] = simStatus.getCountIndividuals( infClasses[kk] );
	}
	RandomGenerator *RNG = simStatus.getRandomGenerator();
	double rnd = RNG->get();
	int otherClass = infClasses[ extractFromDistribution( rnd, pp, -1, 0 ) ];
	int otherId    = RNG->get() * simStatus.getCountIndividuals(otherClass);
	Infos * other = reinterpret_cast<Infos*>( simStatus.getIndividualData( otherId, otherClass ) );
	if (other->isolate == 1)  {
		return false;
	}

	if (other->uniqueId == 0)  {
		other->uniqueId = ++mainUniqueId;
//std::cout << "ID0 " << other->uniqueId << " " << mainUniqueId << "\n";
	}
//std::cout << "ID1 " << other->uniqueId << " " << mainUniqueId << "\n";
	if (tracingOn && RNG->get() < TRACING_PROB[0])  {
		data->contacts.push_back(other->uniqueId);
	}
	return true;
}



bool  infectatwork()  {
	int indivId = simStatus.getCurrentIndividualId();
	int classId = simStatus.getCurrentIndividualClass();
	Infos * data = reinterpret_cast<Infos*>( simStatus.getIndividualData( indivId, classId ) );
	if (data->athome == 1)  return false;
	data->infamily = 1;
	return  infect( indivId, classId, data );
}



bool  infectother()  {
	int indivId = simStatus.getCurrentIndividualId();
	int classId = simStatus.getCurrentIndividualClass();
	Infos * data = reinterpret_cast<Infos*>( simStatus.getIndividualData( indivId, classId ) );
	if (STAYATHOME_FULL[0] == 1.0 && data->athome == 1)  return false;
	data->infamily = 1;
	return  infect( indivId, classId, data );
}



bool  infectnone()  {
	int indivId = simStatus.getCurrentIndividualId();
	int classId = simStatus.getCurrentIndividualClass();
	Infos * data = reinterpret_cast<Infos*>( simStatus.getIndividualData( indivId, classId ) );
	infect( indivId, classId, data );
	return false;
}



void  recovery(Infos *data)  {
	data->uniqueId = 0;
	data->infamily = 0;
	data->isolate  = 0;
	// Clear storage space to reduce memory footprint
	data->contacts.clear();
	data->contacts.shrink_to_fit();
	data->family.clear();
	data->family.shrink_to_fit();
}



// Return boolean because it is used directly in the SEIR model
// besides being used indirectly in the full model
bool  moveToCase()   {
	int indivId = simStatus.getCurrentIndividualId();
	int classId = simStatus.getCurrentIndividualClass();
	Infos * data = reinterpret_cast<Infos*>( simStatus.getIndividualData( indivId, classId ) );
	if (data->uniqueId > 0 && tracingOn)  {
		contactEvents.insert(data->uniqueId);
	}
	totalCases++;
	cumulCases[ ageNames[classId] ]++;
	// Clear storage space to reduce memory footprint
	recovery(data);
	return true;
}



bool moveToHos()  {
	occupancy[OCC_HOS]++;
	return true;
}



bool moveToIcu()  {
	occupancy[OCC_ICU]++;
	return true;
}



bool moveToHom()  {
	occupancy[OCC_HOM]++;
	moveToCase();
	return true;
}



bool  recHos()  {
	occupancy[OCC_HOS]--;
	return true;
}



bool  dedHos()  {
	int classId = simStatus.getCurrentIndividualClass();
	cumulDeaths[ ageNames[classId] ]++;
	occupancy[OCC_HOS]--;
	return true;
}



bool  recIcu()  {
	occupancy[OCC_ICU]--;
	return true;
}



bool  dedIcu()  {
	int classId = simStatus.getCurrentIndividualClass();
	cumulDeaths[ ageNames[classId] ]++;
	occupancy[OCC_ICU]--;
	return true;
}



bool  recHom()  {
	occupancy[OCC_HOM]--;
	return true;
}



bool  recoveryHidden()  {
	int indivId = simStatus.getCurrentIndividualId();
	int classId = simStatus.getCurrentIndividualClass();
	Infos * data = reinterpret_cast<Infos*>( simStatus.getIndividualData( indivId, classId ) );
	// Clear storage space to reduce memory footprint
	recovery(data);
	totalAsyCases++;
	return true;
}



bool  moveToInfect()  {
	return true;
}


#else  // if def TAULEAP

// Defined below.  See end of file

#endif



// Implements measures taking into account their spatial extent
void	evalLocalParameters()  {
	static int  prev_index = -1;
	static int  prev_day = -1;

	int day = static_cast<int>(simStatus.getTime()+0.0000001);
	if (day != prev_day)  {
		prev_day = day;
		prev_index = -1;
		evalDaily();
	}

	int index = evalLocally();
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
			REDUCE_INFLIGHT[0] = REDUCE_INFLIGHT[index];
		}

#ifndef  MODEL_FAMILY
		params.home     = FAMILY_TRANSMIT[0];
#endif
		if (params.stayQuadratic == 1.0)  {
			params.work		= 1.0-STAYATHOME_OTH[0];
		} else {
			params.work		= 1.0;
		}
		params.school	= 1.0-SCHOOL_CLOSURE[0];
		params.other	= 1.0-SOCIALDIST_PROB[0];
		params.mobility = 1.0-TRAVELREDUCTION[0];
		params.pomega	= 1.0-REDUCE_INFLIGHT[0];
		updateContactMatrix();
		prev_index = index;
	}
}




bool  firstinfect( double type, double case_lon, double case_lat, double age )  {
	int  MINAGEGROUP;
	RandomGenerator *gen = simStatus.getRandomGenerator();
	int nAgeGroups = groups[ groups.size()-1 ] + 1;
	std::vector<double>  counts;
	double lat = simStatus.getLatitude();
	double lon = simStatus.getLongitude();
	int pop, kk1 = -1, kk = -1, classId = -1;
	std::string prefix;

//std::cout << "FIRSTINFECT\n";
	auto it = find_if( ages.begin(), ages.end(), [](int jj)  {return static_cast<int>(params.minAgeFirstcase+0.01) < jj;} );
	//MINAGEGROUP = groups[it - ages.begin()-1];
	MINAGEGROUP = 0;

	if (type == 0.0)  {
		prefix = "Esp_";
	} else if (type == 1.0)  {
		prefix = "Inf_";
	} else if (type == 2.0)  {
		prefix = "Asy_";
	} else if (type == 3.0)  {
		prefix = "Esp_";
	} else if (type == 4.0)  {
		prefix = "Inf_";
	} else if (type == 5.0)  {
		prefix = "Asy_";
	} else {
		throw "Unknown type in firstinfect()";
	}

	if (haversine(lon, lat, case_lon, case_lat) < 0.5*GRIDRES)  {
//std::cout << "FIRST INFECTIVE FOUND\n";
		counts.clear();
		pop = simStatus.getLocalPopulationSize();
		for (int jj = MINAGEGROUP; jj < nAgeGroups; jj++)  {
			counts.push_back(
					static_cast<double>( simStatus.getCountIndividuals( classmapper["Sus_"+std::to_string(jj)] ) )
			);
//			std::cout << "&& " << counts[jj-MINAGEGROUP] << "\n";
		}
		double rnd = gen->get();
		kk1 = MINAGEGROUP+extractFromDistribution( rnd, counts, -1, 0 );
		if (age+0.1 >= 0)  {
			kk = static_cast<int>(age+0.1);
			if (kk > nAgeGroups)  simStatus.abort( "Age class out of bounds!" );
			classId = classmapper["Sus_"+std::to_string(kk)];
		}
		if (classId < 0 || simStatus.getCountIndividuals( classId ) == 0)  {
			kk = kk1;
			if (kk > nAgeGroups)  simStatus.abort( "Age class out of bounds!" );
			classId = classmapper["Sus_"+std::to_string(kk)];
		}
		if (classId >= 0)  {
//			std::cout << "&&&&&& [" << kk << "] " << simStatus.getCountIndividuals( classmapper["Sus_"+std::to_string(kk)] ) << "\n";
			std::cout << "CHECK " << baseMap[ kk ][ simStatus.getX() % simStatus.getMaxX() ][ simStatus.getY() % simStatus.getMaxY() ] << " " << simStatus.getCountIndividuals( classId ) << "\n" << std::flush;
			assert( simStatus.getCountIndividuals( classId ) > 0 );

#ifndef  TAULEAP
			Infos *data = reinterpret_cast<Infos*>( simStatus.getIndividualData( 0, classId ) );
			data->uniqueId = ++mainUniqueId;
			data->infamily = 1;
#endif

			simStatus.moveGenericIndividualFromToClass( classId, classmapper[prefix+std::to_string(kk)] );
//std::cout << "FIRSTINFECT MOVED\n";
//			if (type == 2.0)  simStatus.incrementEvents( eventmapper["Case"] );
//			else if (type == 3.0)  simStatus.incrementEvents( eventmapper["Hidden"] );
		} else {
			throw "Cannot find first infective";
		}
		return true;
	} else {
		return false;
	}
}



void  initPolicyParams()  {
	int maxval = 0;
	for (int jj = 0; jj < interventions.size(); jj++)  {
		if (maxval < interventions[jj].getExtent())  maxval = interventions[jj].getExtent();
	}
	TRACING_PROB.assign(maxval+1, 0.0);
	SOCIALDIST_PROB.assign(maxval+1, 0.0); // [1..3] Current values at the different extent levels (above). [0] Current value in grid element
	TRAVELREDUCTION.assign(maxval+1, 0.0);
	TRAVELRED_ADMIN.assign(maxval+1, 0.0);
	STAYATHOME_AGE.assign(maxval+1, 0.0);
	STAYATHOME_OTH.assign(maxval+1, 0.0);
	STAYATHOME_SCH.assign(maxval+1, 0.0);
	FAMILY_TRANSMIT.assign(maxval+1, 1.0);
	STAYATHOME_FULL.assign(maxval+1, 0.0);
	SCHOOL_CLOSURE.assign(maxval+1, 0.0);
	REDUCE_INFLIGHT.assign(maxval+1, 0.0);
}



void  accessCycle( int status )  {
	static bool switcher = true;
	static PolicyQueue queue;
	static DataBuffer  firstInfBuffer;
	double sum = 0;
	std::vector<int>::iterator it;

	switch (status)  {
		case CYCLE_INIT:
			firstInfBuffer.setBuffer( sizeof(int), 1 );
			mainUniqueId = 10000000*simStatus.getProcessId();
			simStatus.setDailyFractions( {8.0, 8.0+(params.weeklyWorkHours/7.0)} );
			getAgeGroups();
			buildAgeAssignments();
			//getSusceptibility();
			initContactMatrix();
			it = find_if( ages.begin(), ages.end(), [](int jj)  {return static_cast<int>(params.youngestLimit+0.01) < jj;} );
			stayAtHomeLimit[0] = groups[it - ages.begin()-1];
			it = find_if( ages.begin(), ages.end(), [](int jj)  {return static_cast<int>(params.eldestLimit+0.01) < jj;} );
			stayAtHomeLimit[1] = groups[it - ages.begin()-1];
std::cout << "STAY AT HOME " << stayAtHomeLimit[0] << " " << stayAtHomeLimit[1] << "\n";

#ifdef  MODEL_FAMILY
			params.home = 0;
			params.R0  *= params.betaMul;
#endif
			params.beta = params.R0 * params.gamma;

			updateContactMatrix();
			initCountrySpecific();
			params.somega = std::exp(params.somega);
			initChineseContactMatrix();
			updateChineseEigenvalue();
			break;

		case CYCLE_START:
			extractFirstInfective();
			initPolicyParams();
			resetCounters();
			queue.clear();
//			addAllPolicies( queue );
			for (auto it = interventions.begin(); it != interventions.end(); it++)  queue.addPolicy( &(*it) );
#ifndef FITTING
			queue.setStart( params.t0 );
#endif
			break;

		case CYCLE_RESET:
//if (simStatus.getProcessId() == PROC_MAIN)  {
//	std::cout << "[" << simStatus.getGlobalProcessId() << "] Policy count [" << queue.size() << "] - Schools day " << simStatus.getTime() << " -> " << SCHOOL_CLOSURE[0] << " " << SOCIALDIST_PROB[0] << "\n";
//}
//			std::cout << "@@@ " << params.eigen <<" " << params.R0 <<" " << 1.0/params.gamma << " " << 1.0 / params.sigma << "\n";
			if (tracingOn)  {
				shareContacts();
			}

			// ###################
			// # SOCIAL DISTANCING
			// ###################
			queue.applyPolicies( -1, simStatus.getTime() );
			//updateContactMatrix();
			break;

		case CYCLE_EVAL:
			evalLocalParameters();
			for (int jj = firstinfected; jj < firstInfections.size(); jj++)  {
				if (firstInfections[firstinfected][1] == 0.0 && simStatus.getTime() >= firstInfections[firstinfected][0])  {
					if (firstinfect( firstInfections[firstinfected][1], firstInfections[firstinfected][3], firstInfections[firstinfected][4], firstInfections[firstinfected][2] ))  {
						firstinfected++;
					}
				} else if (firstInfections[firstinfected][1] == 1.0 && simStatus.getTime() >= firstInfections[firstinfected][0])  {
					if (firstinfect( firstInfections[firstinfected][1], firstInfections[firstinfected][3], firstInfections[firstinfected][4], firstInfections[firstinfected][2] ))  {
						firstinfected++;
					}
				} else if (simStatus.getTime() >= params.t0+firstInfections[firstinfected][0])  {
					if (firstinfect( firstInfections[firstinfected][1], firstInfections[firstinfected][3], firstInfections[firstinfected][4], firstInfections[firstinfected][2] ))  {
						firstinfected++;
					}
				// We ran out of first infections?  Then we do not need to worry anymore
				} else if (simStatus.getTime() >= params.t0+firstInfections[firstInfections.size()-1][0]+1)  {
					firstinfected = firstInfections.size();
				}
			}
			enforcePolicies();
			stayAtHome();
			break;

		case CYCLE_PRE:
			firstInfBuffer.clear();
			firstInfBuffer.pack( &firstinfected, DATA_INT, DATA_MAX );
			firstInfBuffer.reduce();
			firstInfBuffer.unpack( &firstinfected );

			shareOutputData();
			contactEvents.clear();
			break;

		case CYCLE_POST:
#ifndef  FITTING
			if (simStatus.getTime() > 60 + params.t0)  {
				// gcumulXXXX have been set in printmap(), executed before CYCLE_POST
				if (gtotalCases+gtotalAsyCases < 10)  {
					simStatus.exit("This seems an early extinct outbreak");
				}
			}
#endif
			break;

		case CYCLE_LAST:
			break;

		case CYCLE_FINALIZE:
			break;
	}
#ifdef  FITTING
	doFitting(status, queue);
#endif
}



bool  handleIsolation( Infos *data, int classId )  {
	bool  retval;

	if (data->isolate > 1)  {
		retval = false;
	} else if (data->isolate == 1)  {
		retval = true;
	} else {
		retval = false;
	}
	return retval;
}



void  shareContacts()  {
	int  izero = 0;

	int sz = contactEvents.size();
	DataBuffer  buf;
	buf.setBuffer( sizeof(int)*(simStatus.getNumberOfProcesses()+1), (simStatus.getNumberOfProcesses()+1) );
	buf.clear();
	for (int jj = 0; jj < simStatus.getNumberOfProcesses(); jj++)  {
		if (jj == simStatus.getProcessId())  {
			buf.pack( &sz, DATA_INT, DATA_ADD );
		} else {
			buf.pack( &izero, DATA_INT, DATA_ADD );
		}
	}
	buf.reduce();
	int fullsz = 0;
	std::vector<int> sizes( simStatus.getNumberOfProcesses(), 0 );
	for (int jj = 0; jj < simStatus.getNumberOfProcesses(); jj++)  {
		buf.unpack( &sizes[jj] );
		fullsz += sizes[jj];
	}


//std::cout << "FULLSIZE IS: " << fullsz << "\n";
//for (auto &ele: contactEvents)  {
//	std::cout << "##### " << ele << " @@@@@\n";
//}

	buf.setBuffer( sizeof(int)*fullsz, fullsz );
	buf.clear();
	for (int jj = 0; jj < simStatus.getNumberOfProcesses(); jj++)  {
		if (jj == simStatus.getProcessId())  {
			assert( sizes[jj] == contactEvents.size() );
			for (auto it = contactEvents.begin(); it != contactEvents.end(); it++)  {
//		for (int kk = 0,; kk < sizes[jj]; kk++)  {
				int tval = *it;
				buf.pack( &tval, DATA_INT, DATA_ADD );
			}
		} else {
			for (int kk = 0; kk < sizes[jj]; kk++)  {
				buf.pack( &izero, DATA_INT, DATA_ADD );
			}
 		}
	}
	buf.reduce();
	//tmpContactEvents.resize( fullsz );
	contactEvents.clear();
	//tmpContactEvents.clear();
	int  tval;
	for (int jj = 0; jj < fullsz; jj++)  {
		//buf.unpack( &tmpContactEvents[jj] );
		buf.unpack( &tval );
		//tmpContactEvents.insert(tval);
		contactEvents.insert(tval);
//		std::cout << "TMPcontactEVENTS " << tmpContactEvents[jj] << "\n";
	}
	//contactEvents.swap( tmpContactEvents );
//	std::sort( contactEvents.begin(), contactEvents.end() );
//	auto it = std::unique( contactEvents.begin(), contactEvents.end() );
//	contactEvents.resize( std::distance(contactEvents.begin(), it) );

//	for (auto &ele: contactEvents)  {
//		std::cout << "@@@@@ " << ele << " @@@@@\n";
//	}
}



void  stayAtHome()  {
	// Skip if not needed (speed-up)
	if (STAYATHOME_SCH[0] == 0.0 && STAYATHOME_AGE[0] == 0.0 && STAYATHOME_OTH[0] == 0.0)  return;

	RandomGenerator *RNG = simStatus.getRandomGenerator();
	double *rate_athome, prob;
	for (int kk = 0; kk < TOTCLASSES; kk++)  {
		int size = simStatus.getCountIndividuals(kk);
		int age  = ageNames[kk];
		for (int jj = 0; jj < size; jj++)  {
			Infos * data = reinterpret_cast<Infos*>(simStatus.getIndividualData( jj, kk ));
			auto info = simStatus.getIndividualInfo( jj, kk );
			if (info.placeOfContact == 0)  {
				if (age <= stayAtHomeLimit[0])  {
					rate_athome = &STAYATHOME_SCH[0];
				} else if (age >= stayAtHomeLimit[1])  {
					rate_athome = &STAYATHOME_AGE[0];
				} else {
					rate_athome = &STAYATHOME_OTH[0];
				}
				prob = *rate_athome - data->athome_prev;
				if (prob > 0 && data->athome == 0)  {
					if (RNG->get() < prob)  {
						data->athome = 1;
					}
				} else if (prob < 0 && data->athome == 1)  {
					if (RNG->get() < -prob)  {
						data->athome = 0;
					}
				}
				data->athome_prev = *rate_athome;
			}
		}
	}
}


void  enforceIsolation(int jj, int kk, Infos *data)  {
	RandomGenerator *RNG = simStatus.getRandomGenerator();
	if (data->isolate > 1)  {
		data->isolate--;
	} else if (data->isolate == 1)  {
		data->isolate_duration--;
		if (data->isolate_duration == 0)  {
			data->isolate = 0;
		}
	} else if (data->isolate == 0)  {
		for (int qq = 0; qq < data->contacts.size(); qq++)  {
//					if (contactEvents.find( data->contacts[qq] ) != contactEvents.end())  {
			if (std::find(contactEvents.begin(), contactEvents.end(), data->contacts[qq]) != contactEvents.end())  {
				// If found
				//if (RNG->get() < TRACING_PROB[0])  {
					data->isolate = 1+ISOLATE_GAPDAYS;
					data->isolate_duration = 14;
					data->contacts.clear();
					data->contacts.shrink_to_fit();
//std::cout<<"ISOLATING !\n";
					break;
				//}
			}
		}
	}
}



void  enforceFamilyTransmission(int jj, int kk, Infos *data)  {
	if (data->infamily == 1)  {
		auto info = simStatus.getIndividualInfo( jj, kk );
		if (info.placeOfContact == 0)  {
			if (data->family.size() == 0)  {
				familySetup(jj, kk);
			}
			familyTransmission( jj, kk, data );
		}
	}
}



void  enforcePolicies()  {
	for (int kk = 0; kk < TOTCLASSES; kk++)  {
		for (int jj = 0; jj < simStatus.getCountIndividuals(kk); jj++)  {
			Infos * data = reinterpret_cast<Infos*>(simStatus.getIndividualData( jj, kk ));
			if (tracingOn)  {
				enforceIsolation(jj, kk, data);
			}
#ifdef  MODEL_FAMILY
			enforceFamilyTransmission(jj, kk, data);
#endif
		}
	}
}


bool  handleMobility(int x0, int y0, int x1, int y1, double prob)  {
	RandomGenerator *RNG = simStatus.getRandomGenerator();
	bool  retval = false;
	double  rnd, dist = haversine( x0, y0, x1, y1 );
	int  kk;

	x0 = x0 % simStatus.getMaxX();
	x1 = x1 % simStatus.getMaxX();

	if (prob > 0)  {
		if (TRAVELRED_ADMIN[0] > 0.0)  {
			if (checkLockdown(x0, y0) || checkLockdown(x1,y1))  {
				if (RNG->get() < TRAVELRED_ADMIN[0])  {
					return false;
				}
			}
		}
		if (RNG->get() >= params.eta)  return false;
		if (RNG->get() >= params.mobility)  return false;

#ifdef  TAULEAP
		return true;
#else
		int indivId = simStatus.getCurrentIndividualId();
		int classId = simStatus.getCurrentIndividualClass();
		Infos *data = reinterpret_cast<Infos*>( simStatus.getIndividualData( indivId, classId ) );
		if (data->athome == 1)  {
			return false;
		}
		if (data->length == 0)  {
			data->length = getMobilityDuration(dist);
			if (data->length == 0)  return false;
		}
		retval = !handleIsolation(data, classId);
		//contactEvents.clear();
#endif

	} else {

#ifdef  TAULEAP
		return true;
#else
		int indivId = simStatus.getCurrentIndividualId();
		int classId = simStatus.getCurrentIndividualClass();
		Infos *data = reinterpret_cast<Infos*>( simStatus.getIndividualData( indivId, classId ) );
		if (data->isolate == 0)  {
			data->length--;
			if (data->length == 0)  {
				retval = true;
			} else if (data->length > 0) {
				retval = false;
			} else {
				std::cout << data->length << "\n" << std::flush;
				throw "UNEXPECTED ERROR!!!";
			}
		}
#endif

	}
	return retval;
}



void  buildAgeAssignments()  {
	int nAgeGroups = groups[ groups.size()-1 ] + 1;
#ifdef  MODEL_SEIR
	std::vector<std::string>  names = {"Sus_", "Esp_", "Inf_", "Asy_", "Rec_", "Asyrec_"};
#else
	std::vector<std::string>  names = {"Sus_", "Esp_", "Inf_", "Asy_", "Wai_", "Hos_", "Icu_", "Hom_", "Rec_", "Ded_", "Asyrec_"};
#endif
	for (int jj = 0; jj < names.size(); jj++)  {
		for (int kk = 0; kk < nAgeGroups; kk++)  {
			int nm = classmapper[ names[jj] + std::to_string(kk) ];
			int sz = classcounter[ names[jj] + std::to_string(kk) ];
			for (int qq = 0; qq < sz; qq++)  {
				ageNames[nm+qq] = kk;
				TOTCLASSES++;
			}
		}
	}
	std::cout << "Total number of classes is [" << TOTCLASSES << "]\n";

	names = {"Sus_"};
	for (int jj = 0; jj < names.size(); jj++)  {
		for (int kk = 0; kk < nAgeGroups; kk++)  {
			int nm = classmapper[ names[jj] + std::to_string(kk) ];
			int sz = classcounter[ names[jj] + std::to_string(kk) ];
			for (int qq = 0; qq < sz; qq++)  {
				susceptNames[nm+qq] = kk;
			}
		}
	}

//	names = {"Inf_", "Asy_", "Hom_"};
	names = {"Inf_", "Asy_"};
	for (int jj = 0; jj < names.size(); jj++)  {
		for (int kk = 0; kk < nAgeGroups; kk++)  {
			int nm = classmapper[ names[jj] + std::to_string(kk) ];
			int sz = classcounter[ names[jj] + std::to_string(kk) ];
			for (int qq = 0; qq < sz; qq++)  {
				infectNames[nm+qq] = kk;
			}
		}
	}
}


bool  isSusceptible( int classid )  {
	try {
		int xx = susceptNames.at(classid);
	} catch (std::out_of_range &e)  {
		return  false;
	}
	return true;
}


bool  isInfective( int classid )  {
	try {
		int xx = infectNames.at(classid);
	} catch (std::out_of_range &e)  {
		return  false;
	}
	return true;
}


void familySetup(int indivId, int classId)  {
	int nAgeGroups = groups[ groups.size()-1 ] + 1;
	double  rnd;
	RandomGenerator *RNG = simStatus.getRandomGenerator();
	std::vector< std::array<int,3> > familymembers;
	int sz = 0;
	while (sz == 0)  {
		sz = RNG->poisson( params.avgFamilySize );
	}
	std::vector<double> pp;
	int count = 1, notfound = 0;
	auto  info = simStatus.getIndividualInfo( indivId, classId );
	assert( info.placeOfContact == 0);

	familymembers.push_back( {static_cast<int>(info.home.rank), indivId, classId} );
	for (int ii = 0; ii < nAgeGroups; ii++)  {
		pp.push_back( localKK_home[ ii ][ ageNames[classId] ] );
		if (ii > 0)  pp[ii] += pp[ii-1];
	}
	while  (count < sz)  {
		rnd = RNG->get();
		int otherClass = extractFromDistribution( rnd, pp, -1, 1 );
		otherClass = classmapper[ "Sus_" + std::to_string(otherClass) ];
		int otherSz = simStatus.getCountIndividuals( otherClass );
		if (otherSz > 0)  {
			int otherId = RNG->get() * otherSz;
			Infos *data = reinterpret_cast<Infos*>( simStatus.getIndividualData( otherId, otherClass ) );
			auto  otherinfo = simStatus.getIndividualInfo( otherId, otherClass );
			if (otherinfo.placeOfContact == 0 && data->isolate == 0 && data->infamily == 0)  {
				data->infamily = 1;
				count++;
				familymembers.push_back( {static_cast<int>(otherinfo.home.rank), otherId, otherClass} );
			} else {
				notfound++;
			}
		} else {
			notfound++;
		}
		if (notfound > 20)  {
			//simStatus.warning( "familySetup() ran out of individuals! Size is [" + std::to_string( sz ) + "]" );
			break;
		}
	}
	for (int kk = 0; kk < familymembers.size(); kk++)  {
		Infos *data = reinterpret_cast<Infos*>( simStatus.getIndividualData( familymembers[kk][1], familymembers[kk][2] ) );
		for (int qq = 0; qq < familymembers.size(); qq++)  {
			//if (qq == kk)  continue;
			data->family.push_back( familymembers[qq][0] );
		}
	}
}



void  familyTransmission(int indivId, int classId, Infos* data)  {
	IndividualInfo infosthis, infosother;
	RandomGenerator *RNG = simStatus.getRandomGenerator();
	// The factor 1.230 is used to ensure an attack rate of 15.8%
	// Note: we divide by R0 multiplier because in the formula R0 is the well-mixed value
	std::vector<double> zz = {params.zz_0, params.zz_1, params.zz_2, params.zz_3, params.zz_4, params.zz_5, params.zz_6, params.zz_7, params.zz_8};
	std::vector<double> tr = {params.tr_0, params.tr_1, params.tr_2, params.tr_3, params.tr_4, params.tr_5, params.tr_6, params.tr_7, params.tr_8};

	double  prob = params.R0 * params.gamma * FAMILY_TRANSMIT[0] * params.familyAttackMul / params.eigen;
	infosthis = simStatus.getIndividualInfo( indivId, classId );
	int infAge;
	int susAge;
	
	if (infosthis.placeOfContact == 0 && isInfective(classId))  {
		infAge = ageNames[ classId ];
		if (infAge == 1.0 || RNG->get() < tr[infAge])  {
			for (int kk = 0; kk < data->family.size(); kk++)  {
				int rank = data->family[kk];
				if (rank == infosthis.home.rank)  {
					continue;
				}

				try {
					infosother = simStatus.getIndividualInfoByRank(rank);
				} catch (std::exception &e)  {
	//				std::cout << "Not at home\n" << std::flush;
					continue;
				}

	//			std::cout << "This guy is at home!\n" << std::fflush;
				if (isSusceptible(infosother.status))  {
					susAge = ageNames[ infosother.status ];
	//				prob *= localKK_home[ susAge ][ infAge ];
					prob *= zz[susAge];
					if (RNG->get() < 1.0-std::exp( -prob ))  {
						Infos *other = reinterpret_cast<Infos*>( simStatus.getIndividualData( infosother.id, infosother.status ) );
						if (data->uniqueId == 0)  {
							data->uniqueId = ++mainUniqueId;
						}
				
						if (tracingOn && RNG->get() < TRACING_PROB[0])  {
							other->contacts.push_back(data->uniqueId);
						}

						int newClass = ageNames[ infosother.status ];
						newClass = classmapper[ "Esp_" + std::to_string(newClass) ];

	//					std::cout << "MOVING: " << indivId << "/" << simStatus.cell->indivhome->size() << " " << classId << " " << newClass << "\n" << std::flush;
						simStatus.moveSpecificIndividualToClass( infosother.id, infosother.status, newClass );
					}
				}
			}
		}
	}


}




void  shareOutputData()  {
	int nAgeGroups = groups[ groups.size()-1 ] + 1;

	DataBuffer  data;
	data.setBuffer( sizeof(int)*(5+cumulDeaths.size()+cumulCases.size()), (5+cumulDeaths.size()+cumulCases.size()) );	
	data.clear();
	data.pack( &totalCases, DATA_INT, DATA_ADD );
	data.pack( &totalAsyCases, DATA_INT, DATA_ADD );
	for (int kk = 0; kk < 3; kk++)  {
		data.pack( &occupancy[kk], DATA_INT, DATA_ADD );
	}
	for (int kk = 0; kk < nAgeGroups; kk++)  {
		data.pack( &cumulDeaths[kk], DATA_INT, DATA_ADD );
	}
	for (int kk = 0; kk < nAgeGroups; kk++)  {
		data.pack( &cumulCases[kk], DATA_INT, DATA_ADD );
	}
	data.reduce();
	data.unpack( &gtotalCases );
	data.unpack( &gtotalAsyCases );
	for (int kk = 0; kk < 3; kk++)  {
		data.unpack( &goccupancy[kk] );
	}
	for (int kk = 0; kk < nAgeGroups; kk++)  {
		data.unpack( &gcumulDeaths[kk] );
	}
	for (int kk = 0; kk < nAgeGroups; kk++)  {
		data.unpack( &gcumulCases[kk] );
	}
}




void  printextra( FILE *handler, int fd )  {
	int nAgeGroups = groups[ groups.size()-1 ] + 1;

	
	fprintf( handler, "\t%d %d\t%d %d %d\t", gtotalCases, gtotalAsyCases, goccupancy[OCC_HOS], goccupancy[OCC_ICU], goccupancy[OCC_HOM] );
	int sum = 0;
	for (int kk = 0; kk < nAgeGroups; kk++)  {
		sum += gcumulDeaths[kk];
		fprintf( handler, "%d ", gcumulDeaths[kk] );
	}
	fprintf( handler, "%d\t", sum );

	for (int kk = 0; kk < contrTable.size(); kk++)  {
		switch (contrTable[kk])  {
			case POLICY_TRACING_PROB:
				fprintf( handler, " %g", TRACING_PROB[0]*params.tracing );
				break;
			case POLICY_SOCIALDIST_PROB:
				fprintf( handler, " %g", SOCIALDIST_PROB[0] );
				break;
			case POLICY_TRAVELREDUCTION:
				fprintf( handler, " %g", TRAVELREDUCTION[0] );
				break;
			case POLICY_TRAVELRED_ADMIN:
				fprintf( handler, " %g", TRAVELRED_ADMIN[0] );
				break;
			case POLICY_STAYATHOME_AGE:
				fprintf( handler, " %g", STAYATHOME_AGE[0] );
				break;
			case POLICY_STAYATHOME_OTH:
				fprintf( handler, " %g", STAYATHOME_OTH[0] );
				break;
			case POLICY_STAYATHOME_SCH:
				fprintf( handler, " %g", STAYATHOME_SCH[0] );
				break;
			case POLICY_FAMILY_TRANSMIT:
				fprintf( handler, " %g", FAMILY_TRANSMIT[0] );
				break;
			case POLICY_STAYATHOME_FULL:
				fprintf( handler, " %g", STAYATHOME_FULL[0] );
				break;
			case POLICY_SCHOOL_CLOSURE:
				fprintf( handler, " %g", SCHOOL_CLOSURE[0] );
				break;
			case POLICY_REDUCE_INFLIGHT:
				fprintf( handler, " %g", REDUCE_INFLIGHT[0] );
				break;
		}
	}
	sum = 0;
	fprintf( handler, "\t" );
	for (int kk = 0; kk < nAgeGroups; kk++)  {
		sum += gcumulCases[kk];
		fprintf( handler, "%d ", gcumulCases[kk] );
	}
	fprintf( handler, "%d\t", sum );
}



void  printmapextra( FILE *handler )  {
	int nAgeGroups = groups[ groups.size()-1 ] + 1;
	int hosp = 0, icus = 0;

#ifndef MODEL_SEIR
	for (int jj = 0; jj < nAgeGroups; jj++)  {
		hosp += simStatus.getCountIndividuals( classmapper[ "Hos_" + std::to_string(jj) ] );
		icus += simStatus.getCountIndividuals( classmapper[ "Icu_" + std::to_string(jj) ] );
	}
	fprintf( handler, "\t%d %d", hosp, icus );
#endif
}



bool  importedBase()  {
	bool status = importedCase();
//	if (status)  {
//		firstimported = true;
//	}
	return status;
}


int constexpr STOREBUFFER_SZ = 5*1000*1000;

void  doFitting( int status, PolicyQueue &queue )  {
#ifndef BOLFI
	static std::unique_ptr<FittingABC> fitting;
#else
	static std::unique_ptr<FittingBOLFI5> fitting;
#endif
	static Particle    trialParticle;
	static RandomGenerator *deryaRNG;
	static std::vector<std::string> paramNames;
	static std::vector< std::array<int,4> >  timeseries;
	static int cases, deaths, sympt, asympt;
	static int prev_cases, prev_deaths, prev_sympt, prev_asympt;
	static unsigned int tsduration;
	static DataBuffer  chiSqShare, evtData, paramShare;
	static std::ifstream  in;
	static std::ofstream  outputFile, outputFileAvg, outputFitFile;
	static std::string    dirPrefix = "Generation-";
	static bool  hasrestarted = false, isgood;
	static std::array<std::vector<int>, 3>  dailyCount, cumulCount;
	double fzero = 0.0;
	static double chiSquared;
	static int    countEvents, countCases, countAsyCases;
	double dummy;
	static double ctime;
	static std::string outRunName    = "successfulRun-" + to_string( simStatus.getGlobalProcessId() / simStatus.getNumberOfProcesses() ) + ".dat";
	static std::string outDeathsName = "deathsCount-" + to_string( simStatus.getGlobalProcessId() / simStatus.getNumberOfProcesses() ) + ".dat";
	static int check_zero;
	static double tmpnum, tmpden;
	static int deltaParams;
	static bool doFitSusceptibility = false;
	static double officialImportationRate;
	static std::array<int,4> totals = {0};
	static int timer = 0, targetBaseline = 0;
	static bool isfirst = true;

	int nAgeGroups = groups[ groups.size()-1 ] + 1;

	switch( status )  {
		case CYCLE_INIT:
			in.open( timeseriesFile );
			if (in.fail())  {
				simStatus.abort( "Could not open file [" + timeseriesFile + "]" );
			}
			prev_cases = prev_deaths = prev_sympt = prev_asympt = 0;
			do {
				cases = deaths = sympt = asympt = 0;
				for (int jj = 0; jj < inputTable.size(); jj++)  {
					switch( inputTable[jj] )  {
						case DATA_DUMMY:
							in >> dummy;
							break;
						case DATA_ALL_CASES:
							in >> cases;
							break;
						case DATA_CUMUL_ALL_CASES:
							in >> cases;
							cases -= prev_cases;
							break;
						case DATA_DEATHS:
							in >> deaths;
							break;
						case DATA_CUMUL_DEATHS:
							in >> deaths;
							deaths -= prev_deaths;
							break;
						case DATA_SYMPT:
							in >> sympt;
							break;
						case DATA_CUMUL_SYMPT:
							in >> sympt;
							sympt -= prev_sympt;
							break;
						case DATA_ASYMPT:
							in >> asympt;
							break;
						case DATA_CUMUL_ASYMPT:
							in >> asympt;
							sympt -= prev_asympt;
							break;
						default:
							simStatus.abort("Unexpected data type in input file ["+std::to_string(paramTable[jj])+"]");
					}
				}
				timeseries.push_back({deaths, cases, sympt, asympt});
				totals[0] += deaths;
				totals[1] += cases;
				totals[2] += sympt;
				totals[3] += asympt;
				prev_cases  += cases;
				prev_deaths += deaths;
				prev_sympt  += sympt;
				prev_asympt += asympt;
			} while (in.good());
			in.close();
			tsduration = timeseries.size();
			std::cout << "Loaded data for " << tsduration << " days.\n" << std::flush;

			targetBaseline = 0;
			for (int jj = 0; jj < distsTable.size(); jj++)  {
				switch (distsTable[jj])  {
					case  DATA_DEATHS:
						targetBaseline += totals[0];
						break;
					case  DATA_ALL_CASES:
						targetBaseline += totals[1];
						break;
					case  DATA_SYMPT:
						targetBaseline += totals[2];
						break;
					case  DATA_ASYMPT:
						targetBaseline += totals[3];
						break;
					case  DATA_CUMUL_DEATHS:
					case  DATA_CUMUL_ALL_CASES:
					case  DATA_CUMUL_SYMPT:
					case  DATA_CUMUL_ASYMPT:
						simStatus.abort( "CUMULATIVE CASES NOT IMPLEMENTED." );
						break;
				}
			}

			if (params.fitMethod == 1.0)  {
#ifndef BOLFI
				fitting = make_unique<FittingABC>( "Adaptive" );
#endif
			} else if (params.fitMethod == 2.0)  {
#ifdef  BOLFI
				if (params.Restart == 1 || params.Restart == 0 || params.Restart == 2)  {
					fitting = make_unique<FittingBOLFI5>( "Efficient" );
				} else if (params.Restart == 3 || params.Restart == 4 || params.Restart == 5 || params.Restart == 6)  {
					fitting = make_unique<FittingBOLFI5>( "MCMC" );
				}
#endif
			}

			deryaRNG = simStatus.getRandomGenerator();
			fitting->setRNG( deryaRNG );
			if (params.fitMethod == 1)  {
#ifndef  BOLFI
				fitting->setTrials( 1 );  // Trials per particle
				fitting->setDistrs( 1 ); // Number of intermediate distributions after which we assume convergence
				fitting->setParticlesPerTrial( 2 ); // REPLICAS?
				fitting->setDiscardRatio( 0.10, params.fitRate ); // 0.40
				fitting->setThreshold( params.fitThreshold ); // 2,5 %
				fitting->setKernelAmplitude( 1.0 );
#endif
			} else {
				fitting->setThreshold( params.fitThreshold ); // 0.03
			}

			for (auto &elem: paramTable)  {
				switch (elem)  {
					case  PARAM_T0:
						fitting->addParameter( "t0",			0.0,	PARTYPE_UNIFORM,	{ 0., 60.} );
						break;

					case  PARAM_BETA:
						fitting->addParameter( "beta",		0.0,	PARTYPE_UNIFORM,	{std::log(1.0/4.0), std::log(15.0/4.0)} );
						break;

					case  PARAM_R0:
						fitting->addParameter( "R0",			0.0,	PARTYPE_UNIFORM,	{1.0, 8.0} );
						break;

					case  PARAM_GAMMA:
						fitting->addParameter( "1/gamma",	0.0,	PARTYPE_UNIFORM,	{0.2,  8.0} );
						break;

					case  PARAM_TRACING:
						fitting->addParameter( "tracing",	0.0,	PARTYPE_UNIFORM,	{ 0., 1.0} );
						break;

					case  PARAM_TAU:
						fitting->addParameter( "tau",		0.0,	PARTYPE_UNIFORM,	{ 0., 5.0} );
						break;

					case  PARAM_ZMAX:
						fitting->addParameter( "zmax",		0.0,	PARTYPE_UNIFORM,	{ 0.2, 1.0} );
						break;

					case  PARAM_OMEGA:
//						fitting->addParameter( "omega",		0.0,	PARTYPE_UNIFORM,	{ std::log(0.001), std::log(1.) } );
						fitting->addParameter( "omega",		0.0,	PARTYPE_UNIFORM,	{ std::log(0.0001*params.importOmega), std::log(100.*params.importOmega) } );
						break;

					case  PARAM_SUSCEPTIBILITY:
						fitting->addParameter( "z0",			0.0,	PARTYPE_UNIFORM,	{ 0., 1.} );
						fitting->addParameter( "z1",			0.0,	PARTYPE_UNIFORM,	{ 0., 1.} );
						fitting->addParameter( "z2",			0.0,	PARTYPE_UNIFORM,	{ 0., 1.} );
						fitting->addParameter( "z3",			0.0,	PARTYPE_UNIFORM,	{ 0., 1.} );
						//fitting->addParameter( "z4",			0.0,	PARTYPE_UNIFORM,	{ 0., 1.} );
						doFitSusceptibility = true;
						break;

					case  PARAM_SUSCEPTIBILITY_FULL:
						fitting->addParameter( "z0",			0.0,	PARTYPE_UNIFORM,	{ 0., 1.} );
						fitting->addParameter( "z1",			0.0,	PARTYPE_UNIFORM,	{ 0., 1.} );
						fitting->addParameter( "z2",			0.0,	PARTYPE_UNIFORM,	{ 0., 1.} );
						fitting->addParameter( "z3",			0.0,	PARTYPE_UNIFORM,	{ 0., 1.} );
						fitting->addParameter( "z4",			0.0,	PARTYPE_UNIFORM,	{ 0., 1.} );
						fitting->addParameter( "z5",			0.0,	PARTYPE_UNIFORM,	{ 0., 1.} );
						fitting->addParameter( "z6",			0.0,	PARTYPE_UNIFORM,	{ 0., 1.} );
						fitting->addParameter( "z7",			0.0,	PARTYPE_UNIFORM,	{ 0., 1.} );
						doFitSusceptibility = true;
						break;

					case  PARAM_TRANSMISSIBILITY:
						fitting->addParameter( "tr0",			0.0,	PARTYPE_UNIFORM,	{ 0., 1.} );
						fitting->addParameter( "tr1",			0.0,	PARTYPE_UNIFORM,	{ 0., 1.} );
						fitting->addParameter( "tr2",			0.0,	PARTYPE_UNIFORM,	{ 0., 1.} );
						fitting->addParameter( "tr3",			0.0,	PARTYPE_UNIFORM,	{ 0., 1.} );
						//fitting->addParameter( "z4",			0.0,	PARTYPE_UNIFORM,	{ 0., 1.} );
						doFitSusceptibility = true;
						break;

					case  PARAM_TRANSMISSIBILITY_FULL:
						fitting->addParameter( "tr0",			0.0,	PARTYPE_UNIFORM,	{ 0., 1.} );
						fitting->addParameter( "tr1",			0.0,	PARTYPE_UNIFORM,	{ 0., 1.} );
						fitting->addParameter( "tr2",			0.0,	PARTYPE_UNIFORM,	{ 0., 1.} );
						fitting->addParameter( "tr3",			0.0,	PARTYPE_UNIFORM,	{ 0., 1.} );
						fitting->addParameter( "tr4",			0.0,	PARTYPE_UNIFORM,	{ 0., 1.} );
						fitting->addParameter( "tr5",			0.0,	PARTYPE_UNIFORM,	{ 0., 1.} );
						fitting->addParameter( "tr6",			0.0,	PARTYPE_UNIFORM,	{ 0., 1.} );
						fitting->addParameter( "tr7",			0.0,	PARTYPE_UNIFORM,	{ 0., 1.} );
						doFitSusceptibility = true;
						break;
				}
			}
			paramNames = fitting->getNameOfParameters();

			if (doFitSusceptibility)  {
				double  dummy, fval;
				gdistrDeathsSum = 0.0;
				gdistrDeaths.clear();
				std::ifstream in( fileDeathsByAge );
				for (int jj = 0; jj < nAgeGroups; jj++)  {
					in >> dummy;
					in >> fval;
					gdistrDeaths.push_back(fval);
					gdistrDeathsSum += fval;
				}
				in.close();
//				for (int jj = 0; jj < nAgeGroups; jj++)  {
//					gdistrDeaths[jj] /= sum;
//				}
			}


			fitting->init( NN );

			//paramShare.setCommType( DATABUFFER_COMMGLOBAL );
			//chiSqShare.setCommType( DATABUFFER_COMMGLOBAL );
			//evtData.setCommType( DATABUFFER_COMMGLOBAL );
			paramShare.setBuffer( 6*sizeof(double), 6 );
			chiSqShare.setBuffer( 2*sizeof(double)+1*sizeof(int), 2+1 );
			evtData.setBuffer( 3*sizeof(int), 3 );

			if (simStatus.getGlobalProcessId() == 0 && (params.Restart == 1 || params.Restart == 3))  {
				std::system( ("\\rm -rv " + dirPrefix + "*").c_str() );
				// After evaluate, populationTimer has increased -> create dirs for storing future particles
				std::string  dirName = dirPrefix + std::to_string(fitting->getPopulationTimer());
				mkdir( dirName.c_str(), 0755 );
			}

			if (simStatus.getGlobalProcessId() == 0)  {
				if (params.Restart == 1 || params.Restart == 3)  {
					outputFileAvg.open( "outputAvg.dat", std::ios::out );
					outputFile.open( "output.dat", std::ios::out );
					outputFileAvg << "Timer";
					outputFile << "Timer\tSample";
					for (unsigned idx = 0; idx < paramNames.size(); idx++)  {
						outputFileAvg << "\t" << paramNames[idx] << "\t" << "Delta_" << paramNames[idx];
						outputFile << "\t" << paramNames[idx];
					}
					outputFileAvg << "\tError\n" << std::flush;
					outputFile << "\tWeight\n" << std::flush;
#ifdef  BOLFI
					if (params.fitMethod == 2.0)  {
						outputFitFile.open( "outputFit.dat", std::ios::out );
					}
#endif
				} else  {
					outputFileAvg.open( "outputAvg.dat", std::ios::app );
					outputFile.open( "output.dat", std::ios::app );
#ifdef  BOLFI
					if (params.fitMethod == 2.0)  {
						outputFitFile.open( "outputFit.dat", std::ios::app );
					}
#endif
				}
			}
			break;

		case CYCLE_START:
//std::cout << "[" << simStatus.getProcessId() << " START \n" << std::flush;
			if (simStatus.getTime() == 0.0 && params.Restart != 1 && params.Restart != 4)  {
				std::cout << "Loading data from storage.\n" << std::flush;
				MPI_Status  mpiStatus;

				int size, sz;
				std::unique_ptr<char[]> buffer = std::make_unique<char[]>(STOREBUFFER_SZ);
				std::filebuf fb;
				fb.open("storeSnapshot.dat", std::ios::in | std::ios::binary);
				std::istream os(&fb);
				sz = deryaRNG->loadState(nullptr);
				// THIS SHOULD BE USED AS A CHECK. MODIFY
				os.read( reinterpret_cast<char*>(&size), sizeof(int) );
				std::unique_ptr<char[]>  data = std::make_unique<char[]>(size);
				for (int jj = 0; jj < simStatus.getGlobalNumberOfProcesses(); jj++)  {
					os.read( data.get(), size );
					if (simStatus.getGlobalProcessId() == jj)  {
						// Do not overwrite RNG status if forecasting (we need different seeds)
//						if (params.Restart != 5)  {
							sz = deryaRNG->loadState( data.get() );
							if (sz != size)  {
								simStatus.abort("Error in reading stored RNG data.");
							}
//						}
					}
				}

				os.read( reinterpret_cast<char*>(&size), sizeof(int) );
				// THIS SHOULD BE USED AS CHECK. MODIFY
				os.read( buffer.get(), size );
				sz = fitting->loadState( buffer.get() );
				if (size != sz)  {
					simStatus.abort("Error in reading state data.");
				}
				if (size > STOREBUFFER_SZ)  {
					simStatus.abort("Size of stored data too large!");
				}
				os.read( reinterpret_cast<char*>(&timer), sizeof(int) );
				fb.close();

			} else if (simStatus.getTime() == 0.0 && ((params.Restart == 1 && !isfirst) || params.Restart == 4)) {
				std::cout << "Saving data to storage.\n" << std::flush;
				int  sz;

				MPI_Status  mpiStatus;
				sz = deryaRNG->saveState(nullptr);
				std::unique_ptr<char[]>  data = std::make_unique<char[]>(sz);
				if (simStatus.getGlobalProcessId() != 0)  {
					deryaRNG->saveState( data.get() );
					MPI_Send( data.get(), sz, MPI_BYTE, 0, 123456+simStatus.getGlobalProcessId(),  MPI_COMM_WORLD);
				}

				if (simStatus.getGlobalProcessId() == 0)  {
					std::unique_ptr<char[]> buffer = std::make_unique<char[]>(STOREBUFFER_SZ);
					std::filebuf fb;
					fb.open("storeSnapshot.dat", std::ios::out | std::ios::binary);
					std::ostream os(&fb);

					// Random generator status
					int size = deryaRNG->saveState(nullptr);
					os.write( reinterpret_cast<char*>(&size), sizeof(int) );
					deryaRNG->saveState( buffer.get() );
					os.write( buffer.get(), size );
					for (int jj = 1; jj < simStatus.getGlobalNumberOfProcesses(); jj++)  {
						MPI_Recv( data.get(), size, MPI_BYTE, jj, 123456+jj, MPI_COMM_WORLD, &mpiStatus);
						os.write( data.get(), size );
					}

					// Fitting status
					size = fitting->saveState( buffer.get() );
					if (size > STOREBUFFER_SZ)  exit(1);
					os.write( reinterpret_cast<char*>(&size), sizeof(int) );
					os.write( buffer.get(), size );
					os.write( reinterpret_cast<char*>(&timer), sizeof(int) );

					fb.close();
				}
			}
			isfirst = false;

//std::cout << "##### RESTART is [" << params.Restart << "]\n";
			if (params.Restart == 0)  {
				params.Restart = 1;
			}  else if (params.Restart == 2)  {
				fitting->reset();
				params.Restart = 1;
			}  else if (params.Restart == 3 || params.Restart == 5)  {
//				fitting->init( 1 );
//std::cout << "##### UPDATE\n";
				fitting->setUpdate( false );
				params.Restart = 4;
			}  else if (params.Restart == 6)  {
				fitting->setUpdate( false );
			}

			hasrestarted = false;

do {
			dailyCount[0].clear();
			dailyCount[1].clear();
			dailyCount[2].clear();
			check_zero = 0;
			params.R0 = 0.0;
			params.tracing = 1.0;
			params.pomega  = 1.0;
			trialParticle = fitting->start();
			deltaParams = 0;
			for (int jj = 0; jj < paramTable.size(); jj++)  {
				switch (paramTable[jj])  {
					case  PARAM_T0:
						params.t0      = trialParticle.parameters[jj+deltaParams];
						break;

					case  PARAM_BETA:
						params.beta    = std::exp( trialParticle.parameters[jj+deltaParams] );
						break;

					case  PARAM_R0:
						params.R0      = trialParticle.parameters[jj+deltaParams];
						break;

					case  PARAM_GAMMA:
						params.gamma   = 1.0/trialParticle.parameters[jj+deltaParams];
						break;

					case  PARAM_TRACING:
						params.tracing = trialParticle.parameters[jj+deltaParams];
						break;

					case  PARAM_TAU:
						params.tau = trialParticle.parameters[jj+deltaParams];
						break;

					case  PARAM_ZMAX:
						params.zmax = trialParticle.parameters[jj+deltaParams];
						break;

					case  PARAM_OMEGA:
						params.somega = std::exp(trialParticle.parameters[jj+deltaParams]);
						break;

					case  PARAM_SUSCEPTIBILITY:
						params.zz_0 = trialParticle.parameters[jj];
						params.zz_1 = trialParticle.parameters[jj+1];
						params.zz_2 = trialParticle.parameters[jj+2];
						params.zz_3 = params.zz_2;
						params.zz_4 = trialParticle.parameters[jj+3];
						params.zz_5 = params.zz_4;
						params.zz_6 = params.zz_4;
						//params.zz_7 = trialParticle.parameters[jj+4];
						//params.zz_8 = params.zz_5;
						deltaParams = 3;
						params.zz_norm = 9.0/(params.zz_0+params.zz_1+params.zz_2+params.zz_3+params.zz_4+params.zz_5+params.zz_6+params.zz_7+params.zz_8);
						updateChineseEigenvalue();
						break;

					case  PARAM_SUSCEPTIBILITY_FULL:
						params.zz_0 = trialParticle.parameters[jj];
						params.zz_1 = trialParticle.parameters[jj+1];
						params.zz_2 = trialParticle.parameters[jj+2];
						params.zz_3 = trialParticle.parameters[jj+3];
						params.zz_4 = trialParticle.parameters[jj+4];
						params.zz_5 = trialParticle.parameters[jj+5];
						params.zz_6 = trialParticle.parameters[jj+6];
						params.zz_7 = trialParticle.parameters[jj+7];
						//params.zz_7 = trialParticle.parameters[jj+4];
						//params.zz_8 = params.zz_5;
						deltaParams = 7;
						params.zz_norm = 9.0/(params.zz_0+params.zz_1+params.zz_2+params.zz_3+params.zz_4+params.zz_5+params.zz_6+params.zz_7+params.zz_8);
						updateChineseEigenvalue();
						break;

					case  PARAM_TRANSMISSIBILITY:
						params.tr_0 = trialParticle.parameters[jj];
						params.tr_1 = trialParticle.parameters[jj+1];
						params.tr_2 = trialParticle.parameters[jj+2];
						params.tr_3 = params.tr_2;
						params.tr_4 = trialParticle.parameters[jj+3];
						params.tr_5 = params.tr_4;
						params.tr_6 = params.tr_4;
						//params.zz_7 = trialParticle.parameters[jj+4];
						//params.zz_8 = params.zz_5;
						deltaParams = 3;
						params.tr_norm = 9.0/(params.tr_0+params.tr_1+params.tr_2+params.tr_3+params.tr_4+params.tr_5+params.tr_6+params.tr_7+params.tr_8);
						updateChineseEigenvalue();
						break;

					case  PARAM_TRANSMISSIBILITY_FULL:
						params.tr_0 = trialParticle.parameters[jj];
						params.tr_1 = trialParticle.parameters[jj+1];
						params.tr_2 = trialParticle.parameters[jj+2];
						params.tr_3 = trialParticle.parameters[jj+3];
						params.tr_4 = trialParticle.parameters[jj+4];
						params.tr_5 = trialParticle.parameters[jj+5];
						params.tr_6 = trialParticle.parameters[jj+6];
						params.tr_7 = trialParticle.parameters[jj+7];
						//params.zz_7 = trialParticle.parameters[jj+4];
						//params.zz_8 = params.zz_5;
						deltaParams = 7;
						params.tr_norm = 9.0/(params.tr_0+params.tr_1+params.tr_2+params.tr_3+params.tr_4+params.tr_5+params.tr_6+params.tr_7+params.tr_8);
						updateChineseEigenvalue();
						break;

				}
			}
			params.t0   = std::floor( params.t0 );

#ifdef  MODEL_FAMILY
			if (params.R0 == 0.0)  {
				params.beta *= params.betaMul;
				params.R0    = params.beta/params.gamma;
			} else {
				params.R0  *= params.betaMul;
				params.beta = params.R0 * params.gamma;
			}
#else
			if (params.R0 == 0.0)  {
				params.R0    = params.beta/params.gamma;
			} else {
				params.beta = params.R0 * params.gamma;
			}
#endif

			if (simStatus.getProcessId() == 0)  {
				std::cout << "[" << simStatus.__fullProcessId << "] Params: ";
				for (int jj = 0; jj < paramNames.size(); jj++)  {
					std::cout << paramNames[jj] << " " << trialParticle.parameters[jj] << ",";
				}
				std::cout << " -> " << params.t0 + timeseries.size() << "\n" << std::flush;
			}

} while (params.Restart == 6);

//std::cout << "CHECK PP: {" << params.zz_4 << ", " << params.zz_5 << ", " << params.zz_6 << ", " << params.zz_7 << ", " << params.zz_8 << "}\n" << std::flush;
//std::cout << "CHECK N0\n" << std::flush;
			fitting->acceptParticle(trialParticle);
//std::cout << "CHECK N1\n" << std::flush;



			chiSquared = 0.0;
			countEvents = 0;
			countCases  = 0;
			countAsyCases = 0;
			tmpnum = tmpden = 0.0;
//std::cout << "CHECK N1\n" << std::flush;
			simStatus.setSimulationLength( params.t0 + timeseries.size() );
//std::cout << "CHECK N2\n" << std::flush;
			queue.setStart( params.t0 );
//std::cout << "CHECK N3\n" << std::flush;
			break;

		case CYCLE_RESET:
//			if (!globalFirstimported)  {
//				countDaysDelay++;
//				evtData.clear();
//				evtData.pack( &firstimported, DATA_BOOL, DATA_OR );
//				evtData.reduce();
//				evtData.unpack( &globalFirstimported );
//				if (globalFirstimported && simStatus.getGlobalProcessId() == 0)  {
//					std::cout << "Delay is [" << countDaysDelay << "] days.\n";
//				}
//			}
			break;

		case CYCLE_EVAL:
			break;

		case CYCLE_PRE:
//			if (!globalFirstimported)  {
//				simStatus.preventOutputlineFrame(0);
//				simStatus.preventOutputlineFrame(1);
//				simStatus.preventOutputmapFrame();
//			} else  {
//				params.somega = officialImportationRate;
				if (simStatus.isOutputlineFrame(0))  {
#ifdef  MODEL_SEIR
					countCases  = simStatus.getCountEvents( 0, eventmapper["Case"] );
					evtData.clear();
					evtData.pack( &countCases,   DATA_INT, DATA_ADD );
					evtData.reduce();
					evtData.unpack( &countCases );
					dailyCount[0].push_back( 0 );
					dailyCount[1].push_back( countCases );
#else
					countEvents = simStatus.getCountEvents( 0, eventmapper["Deaths"] );
					countCases  = simStatus.getCountEvents( 0, eventmapper["CaseHos"] ) + simStatus.getCountEvents( 0, eventmapper["CaseIcu"] ) + simStatus.getCountEvents( 0, eventmapper["CaseHom"] );
					countAsyCases = simStatus.getCountEvents( 0, eventmapper["Hidden"] );
					evtData.clear();
					evtData.pack( &countEvents,  DATA_INT, DATA_ADD );
					evtData.pack( &countCases,   DATA_INT, DATA_ADD );
					evtData.pack( &countAsyCases,   DATA_INT, DATA_ADD );
					evtData.reduce();
					evtData.unpack( &countEvents );
					evtData.unpack( &countCases );
					evtData.unpack( &countAsyCases );

					dailyCount[0].push_back( countEvents );
					dailyCount[1].push_back( countCases );
					dailyCount[2].push_back( countAsyCases );
#endif
				}
//std::cout << "[" << simStatus.getProcessId() << "] START 1\n" << std::flush;
//			}
			break;

		case CYCLE_POST:
//			if (!globalFirstimported)  {
//				simStatus.setTime( 0 );
//			} else {
				ctime = simStatus.getTime();
				if (simStatus.isOutputlineFrame(0))  {
					if (ctime > params.t0)  {
						double numerator = 0.0, denominator = 0.0;
						double chiSquaredTmp = 0.0;
						int day = static_cast<int>(ctime-params.t0+0.001);

						if (day < timeseries.size())  {
							for (int jj = 0; jj < distsTable.size(); jj++)  {
								switch (distsTable[jj])  {
									case  DATA_DEATHS:
										tmpnum = countEvents;
										tmpden = countEvents;
										tmpnum += -timeseries[ day ][0];
										tmpden +=  timeseries[ day ][0];
										check_zero += countEvents;
										break;
									case  DATA_ALL_CASES:
										tmpnum = countCases + countAsyCases;
										tmpden = countCases + countAsyCases;
										tmpnum += -timeseries[ day ][1];
										tmpden +=  timeseries[ day ][1];
										check_zero += countCases + countAsyCases;
										break;
									case  DATA_SYMPT:
										tmpnum = countCases;
										tmpden = countCases;
										tmpnum += -timeseries[ day ][2];
										tmpden +=  timeseries[ day ][2];
										check_zero += countCases;
										break;
									case  DATA_ASYMPT:
										tmpnum = countAsyCases;
										tmpden = countAsyCases;
										tmpnum += -timeseries[ day ][3];
										tmpden +=  timeseries[ day ][3];
										check_zero += countAsyCases;
										break;
									case  DATA_CUMUL_DEATHS:
										tmpnum += countEvents;
										tmpden += countEvents;
										tmpnum += -timeseries[ day ][0];
										tmpden +=  timeseries[ day ][0];
										check_zero += countEvents;
										break;
									case  DATA_CUMUL_ALL_CASES:
										tmpnum += countCases + countAsyCases;
										tmpden += countCases + countAsyCases;
										tmpnum += -timeseries[ day ][1];
										tmpden +=  timeseries[ day ][1];
										check_zero += countCases + countAsyCases;
										break;
									case  DATA_CUMUL_SYMPT:
										tmpnum += countCases;
										tmpden += countCases;
										tmpnum += -timeseries[ day ][2];
										tmpden +=  timeseries[ day ][2];
										check_zero += countCases;
										break;
									case  DATA_CUMUL_ASYMPT:
										tmpnum += countAsyCases;
										tmpden += countAsyCases;
										tmpnum += -timeseries[ day ][3];
										tmpden +=  timeseries[ day ][3];
										check_zero += countAsyCases;
										break;

								}
								numerator   += tmpnum*tmpnum;
								denominator += tmpden*tmpden;

								chiSquaredTmp = denominator > 0 ? (numerator) / sqrt(denominator) : 0.0;
								chiSquared += chiSquaredTmp;
							}
						}
					}
					hasrestarted = fitting->step( chiSquared );
					if (hasrestarted)  {
						if (simStatus.getProcessId() == 0)  {
							std::cout << "Restarting in [Replica "+std::to_string( simStatus.getGlobalProcessId()/simStatus.getNumberOfProcesses() )+" ]\n" << std::flush;
						}
					}
					countEvents = 0;
					countCases  = 0;
					countAsyCases = 0;
				}
//			}
			break;

		case CYCLE_LAST:
//std::cout << "[" << simStatus.getGlobalProcessId() << "] - CYCLE_LAST start \n";
			if (doFitSusceptibility)   {
				double numerator = 0.0, denominator = 0.0;
				double chiSquaredTmp = 0.0;
				double sum = 0;
				for (int kk = 0; kk < nAgeGroups; kk++)  {
					sum += gcumulDeaths[kk];
				}
				for (int kk = 0; kk < nAgeGroups; kk++)  {
					gcumulDeaths[kk] /= sum;
				}
				for (int kk = 0; kk < nAgeGroups; kk++)  {
					tmpnum  = gcumulDeaths[kk]*gdistrDeathsSum - gdistrDeaths[kk];
					tmpden  = gcumulDeaths[kk]*gdistrDeathsSum + gdistrDeaths[kk];
					numerator   += tmpnum*tmpnum;
					denominator += tmpden*tmpden;
				}
				chiSquaredTmp = denominator > 0 ? (numerator) / sqrt(denominator) : 0.0;
				// Add and balance with other summary statistics
				chiSquared = chiSquared + chiSquaredTmp * static_cast<double>(targetBaseline) / static_cast<double>(gdistrDeathsSum);
			}


			if (params.allowZeroFits == 0.0 && check_zero == 0)  chiSquared = 2*fitting->getError();

			if (simStatus.getProcessId() == 0 && chiSquared < fitting->getError() && (params.Restart == 1 || params.Restart == 4))  {
				// Eval cumuls
				cumulCount = dailyCount;
				for (unsigned jj = 0; jj < 3; jj++)  {
					for (unsigned idx = 1; idx < cumulCount[0].size(); idx++)  {
						cumulCount[jj][idx] += cumulCount[jj][idx-1];
					}
				}

				std::string  dirName = dirPrefix + std::to_string(fitting->getPopulationTimer());
				std::ofstream  outRun = std::ofstream( dirName + "/" + outRunName, ios::app );
				for (unsigned idx = 0; idx < dailyCount[0].size(); idx++)  {
					outRun << (static_cast<int>(idx+1) - static_cast<int>(params.t0+0.1));
					for (int jj = 0; jj < printTable.size(); jj++)  {
						switch (printTable[jj])  {
							case  DATA_DEATHS:
								outRun << "\t" << dailyCount[0].at(idx);
								break;
							case  DATA_ALL_CASES:
								outRun << "\t" << dailyCount[1].at(idx) + dailyCount[2].at(idx);
								break;
							case  DATA_SYMPT:
								outRun << "\t" << dailyCount[1].at(idx);
								break;
							case  DATA_ASYMPT:
								outRun << "\t" << dailyCount[2].at(idx);
								break;
							case  DATA_CUMUL_DEATHS:
								outRun << "\t" << cumulCount[0].at(idx);
								break;
							case  DATA_CUMUL_ALL_CASES:
								outRun << "\t" << cumulCount[1].at(idx) + cumulCount[2].at(idx);
								break;
							case  DATA_CUMUL_SYMPT:
								outRun << "\t" << cumulCount[1].at(idx);
								break;
							case  DATA_CUMUL_ASYMPT:
								outRun << "\t" << cumulCount[2].at(idx);
								break;
						}
					}
					outRun << "\n";
				}
				outRun << "\n";
				outRun.close();

				std::ofstream  outDeaths = std::ofstream( dirName + "/" + outDeathsName, ios::app );
				for (unsigned idx = 0; idx < gdistrDeaths.size(); idx++)  {
					if (idx > 0)  outDeaths << " ";
					outDeaths << gdistrDeaths[idx];
				}
				outDeaths << "\n";
				outDeaths.close();
			}

//std::cout << "[" << simStatus.getGlobalProcessId() << "] - CYCLE_LAST CHECK 0 \n";
			isgood = fitting->evaluate( chiSquared, hasrestarted );
//std::cout << "[" << simStatus.getGlobalProcessId() << "] - CYCLE_LAST CHECK 1 \n";
			if (simStatus.getGlobalProcessId() == 0)  {
				int accepted, total;
				double  ratio = fitting->getCurrentAcceptanceRatio( &accepted, &total );
				std::cout << "Current ratio: [" << (accepted) << "/" << (total) << "]";
				std::cout << " -> [" << (100.0*(ratio)) << "%] ";
				ratio = fitting->getAcceptanceRatio( &accepted, &total );
				std::cout << "  (LAST WAS) [" << (accepted) << "/" << (total) << "]";
				std::cout << " -> [" << (100.0*(ratio)) << "%]\n" << std::flush;

			}

//std::cout << "[" << simStatus.getGlobalProcessId() << "] - CYCLE_LAST CHECK 2 \n";
			if (isgood && !hasrestarted && simStatus.getGlobalProcessId() == 0 && params.fitMethod == 2.0)  {
				std::vector< std::vector<double> >  validpars = fitting->getAcceptedParameters();
				for (unsigned int idx = 0; idx < validpars.size(); idx++)  {
					outputFile << fitting->getPopulationTimer();
					for (unsigned int kk = 0; kk < validpars[idx].size(); kk++)  {
						//if (kk > 0)  outputFile << "\t";
						outputFile << "\t" << validpars[idx][kk];
					}
					outputFile << "\n" << std::flush;
				}
			}

			if (hasrestarted)  {
//std::cout << "[" << simStatus.getGlobalProcessId() << "] - CYCLE_LAST CHECK 3 \n";
				int  timeIndex;

				if (simStatus.getGlobalProcessId() == 0)  {
					// After evaluate, populationTimer has increased -> create dirs for storing future particles
					std::string  dirName = dirPrefix + std::to_string(fitting->getPopulationTimer());
					mkdir( dirName.c_str(), 0755 );

					std::vector<double> avg = fitting->getAverages();
					std::vector<double> var = fitting->getVariances();
					outputFileAvg << avg[0];
					for (unsigned int idx = 1; idx < avg.size(); idx++)  {
						outputFileAvg << "\t" << avg[idx] << "\t" << var[idx];
					}
					outputFileAvg << "\t" << fitting->getError() << "\n" << std::flush;
					timeIndex = static_cast<int>(floor(avg[0] + 0.0001));

					std::vector< std::vector<double> >  validpars = fitting->getAcceptedParameters();
					for (unsigned int idx = 0; idx < validpars.size(); idx++)  {
						outputFile << fitting->getPopulationTimer();
						for (unsigned int kk = 0; kk < validpars[idx].size(); kk++)  {
							//if (kk > 0)  outputFile << "\t";
							outputFile << "\t" << validpars[idx][kk];
						}
						outputFile << "\n" << std::flush;
					}

#ifdef  BOLFI
					if (params.fitMethod == 2)  {
						std::vector<double>  fittingpars = fitting->getLikelihoodPoints();
						outputFitFile << fitting->getPopulationTimer();
						for (unsigned int idx = 0; idx < fittingpars.size(); idx++)  {
							//if (kk > 0)  outputFile << "\t";
							outputFitFile << "\t" << fittingpars[idx];
						}
						outputFitFile << "\n" << std::flush;
					}
#endif

				}

				if (simStatus.getGlobalProcessId() == 0)  {
					std::vector<double> avg = fitting->getAverages();
					std::vector<double> var = fitting->getVariances();
					std::cout << "VALS " << avg[0];
					for (unsigned int idx = 1; idx < avg.size(); idx++)  {
						std::cout << "\t" << avg[idx] << "\t" << var[idx];
					}
					std::cout << "\t" << fitting->getError() << "\n" << std::flush;
				}

			}

//std::cout << "CHECK A.1\n" << std::flush;
			if ((params.Restart == 1 || params.Restart == 4) && fitting->getPopulationSamplingTimer() <= TT)  {
//std::cout << "CHECK A.2\n" << std::flush;
				simStatus.setTime( 0 );
//std::cout << "CHECK A.3\n" << std::flush;
				simStatus.reinitDisease();
//std::cout << "CHECK A.4\n" << std::flush;
				if (fitting->getPopulationSamplingTimer() > 0)  {
					simStatus.setSimulationLength( params.forecastLen );
				}
//std::cout << "CHECK A.5\n" << std::flush;
				if (params.Restart == 4 && simStatus.getGlobalProcessId() == 0)  {
					char str[200];
					sprintf( str, "mkdir -p Run-%03d; cp -v * Run-%03d;", timer, timer );
					timer++;
					system( str );
				}
			}
//std::cout << "CHECK A.6\n" << std::flush;
			break;

		case CYCLE_FINALIZE:
			if (simStatus.getProcessId() == 0)  {
				outputFileAvg.close();
				outputFile.close();
			}
			fitting->finalize();
			//delete fitting;
			break;

	}
}






#ifdef  TAULEAP

bool  infectother()  {
	return true;
}

int  xinfectother(int nn)  {
	return nn;
}



bool  infectnone()  {
	return false;
}

int  xinfectnone(int nn)  {
	return 0;
}



bool  infectatwork()  {
	return true;
}

int  xinfectatwork(int nn)  {
	return nn;
}




bool  moveToInfect()  {
	return true;
}

int  xmoveToInfect(int nn)  {
	return nn;
}




bool  recoveryHidden()  {
	totalAsyCases++;
	return true;
}

int  xrecoveryHidden(int nn)  {
	totalAsyCases += nn;
	return nn;
}



// Return boolean because it is used directly in the SEIR model
// besides being used indirectly in the full model
bool  moveToCase()   {
	totalCases++;
	return true;
}

int  xmoveToCase(int nn)  {
	totalCases += nn;
	return nn;
}




bool moveToHos()  {
	occupancy[OCC_HOS]++;
	return true;
}

int  xmoveToHos(int nn)  {
	occupancy[OCC_HOS] += nn;
	return nn;
}




bool moveToIcu()  {
	occupancy[OCC_ICU]++;
	return true;
}

int  xmoveToIcu(int nn)  {
	occupancy[OCC_ICU] += nn;
	return nn;
}




bool moveToHom()  {
	occupancy[OCC_HOM]++;
	moveToCase();
	return true;
}

int xmoveToHom(int nn)  {
	occupancy[OCC_HOM] += nn;
	xmoveToCase(nn);
	return nn;
}




bool  recHos()  {
	occupancy[OCC_HOS]--;
	return true;
}

int  xrecHos(int nn)  {
	occupancy[OCC_HOS] -= nn;
	return nn;
}




bool  dedHos()  {
	//int classId = simStatus.getCurrentIndividualClass();
	//cumulDeaths[ ageNames[classId] ]++;
	occupancy[OCC_HOS]--;
	return true;
}

int  xdedHos(int nn)  {
	occupancy[OCC_HOS] -= nn;
	// Missing
	return nn;
}




bool  recIcu()  {
	occupancy[OCC_ICU]--;
	return true;
}

int  xrecIcu(int nn)  {
	occupancy[OCC_ICU] -= nn;
	return nn;
}




bool  dedIcu()  {
	//int classId = simStatus.getCurrentIndividualClass();
	//cumulDeaths[ ageNames[classId] ]++;
	occupancy[OCC_ICU]--;
	return true;
}

int  xdedIcu(int nn)  {
	//Missing
	occupancy[OCC_ICU] -= nn;
	return nn;
}




bool  recHom()  {
	occupancy[OCC_HOM]--;
	return true;
}

int  xrecHom(int nn)  {
	occupancy[OCC_HOM] -= nn;
	return nn;
}

#endif




std::vector< std::vector<double> >  reduceContactMatrix( std::vector< std::vector<double> >  &mm, std::vector<int>  &aP )  {
	std::vector< std::vector<double> >  KKreduced;
	KKreduced.resize(9, std::vector<double>(9, 0.0) );

	for (int xx=0; xx < 16; xx++)  {
		mm[xx][16] = mm[xx][15];
	}
	for (int xx=0; xx < 17; xx++)  {
		mm[16][xx] = mm[15][xx];
	}

	for (int yy=0; yy < 16; yy++)  {
		for (int xx=0; xx < 16; xx++)  {
			mm[xx][yy] *= aP[xx];
		}
	}
	for (int xx=0; xx < 8; xx++)  {
		for (int yy=0; yy < 8; yy++)  {
			KKreduced[xx][yy] += mm[2*xx][2*yy]+mm[2*xx+1][2*yy]+mm[2*xx][2*yy+1]+mm[2*xx+1][2*yy+1];
			KKreduced[xx][yy] /= aP[2*xx]+aP[2*xx+1];
		}
		KKreduced[xx][8] += mm[2*xx][16]+mm[2*xx+1][16];
		KKreduced[xx][8] /= aP[2*xx]+aP[2*xx+1];
	}
	for (int yy=0; yy < 8; yy++)  {
		KKreduced[8][yy] += mm[16][2*yy]+mm[16][2*yy+1];
		KKreduced[8][yy] /= aP[16];
	}
	return  KKreduced;
}




void  initChineseContactMatrix()  {
	int nAgeGroups = groups[ groups.size()-1 ] + 1;
	double fval;
	chineseKK = arma::zeros( nAgeGroups, nAgeGroups );
	std::vector< std::string >  types = {"home", "work", "school", "other"};
	for (int kk = 0; kk < types.size(); kk++)  {
		std::ifstream infile("../../../../Data/China/Contacts/ChinaContactMatrix_" + types[kk] + "_g09.csv");
		for (int ii = 0; ii < nAgeGroups; ii++)  {
			for (int jj = 0; jj < nAgeGroups; jj++)  {
				infile >> fval;
				chineseKK(ii, jj) += fval;
			}
		}
	}
}



void  updateChineseEigenvalue()  {
	int nAgeGroups = groups[ groups.size()-1 ] + 1;
	double eigval;
	arma::vec  zz(nAgeGroups), tr(nAgeGroups);
	zz(0) = params.zz_0;
	zz(1) = params.zz_1;
	zz(2) = params.zz_2;
	zz(3) = params.zz_3;
	zz(4) = params.zz_4;
	zz(5) = params.zz_5;
	zz(6) = params.zz_6;
	zz(7) = params.zz_7;
	zz(8) = params.zz_8;

	tr(0) = params.tr_0;
	tr(1) = params.tr_1;
	tr(2) = params.tr_2;
	tr(3) = params.tr_3;
	tr(4) = params.tr_4;
	tr(5) = params.tr_5;
	tr(6) = params.tr_6;
	tr(7) = params.tr_7;
	tr(8) = params.tr_8;
	arma::mat KK = chineseKK;
	KK = zz % KK.each_col();  // elementwise multiply each col
	KK = KK.each_row() % tr.t();
	arma::cx_vec eigvals = arma::eig_gen( KK );
	params.eigen = 0;
	for (int jj = 0; jj < eigvals.size(); jj++)  {
		eigval = std::sqrt( norm(eigvals(jj)) );
		params.eigen = max(params.eigen, eigval);
	}
	if (simStatus.getProcessId() == 0)  {
		std::cout << "China eigenvalue is [" << params.eigen << "]\n" << std::flush;
	}
}



