#include <fstream>
#include "GeoTiffReader.h"
#include <assert.h>
#include <queue>
#include <algorithm>
#include <set>
#include <sys/stat.h>  // mkdir (Linux)
#include <memory>  // smart pointers

#include "SplineInterpolator.h"
#include "Policy.h"

#ifdef FITTING
#warning( "Macro FITTING was defined." )
#endif

enum  {DET_0_00, DET_0_10, DET_0_25, DET_0_50, DET_1_00};

// Used in fitting procedure
static constexpr int NN = 16;  // Number of particles
static constexpr int TT = 1;   // Number of post-convergence generations

// If this is NOT a fitting run
int    constexpr DET_RATE = DET_1_00;	// Transmission rate for asymptomatics
double constexpr SYMPTOMATIC_MULTIPLIER = 1.0; // Transmission rate for oldest class

int constexpr ISOLATE_GAPDAYS = 4;   // Gap days from probable case detection to isolation of traced individuals
int constexpr TRACING_THRESHOLD = 100;  // When contact tracing stops
double constexpr  TRACING_PROB = 0.00;  // Contact tracing effectiveness 

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

static bool tracingOn = false;   // Flags if tracing is being carried on
static int  totalCases = 0;     // Process total number of cases
static int  cumulCases = 0;     // Process cumul number of cases
static int  gcumulCases = 0;    // Global cumul number of cases (for output)
static int  cumulAsyCases = 0;  // Process cumul number of asymptomatic cases
static int  gcumulAsyCases = 0; // Global cumul number of asymptomatic cases (for output)
static int  mainUniqueId = 0;   // Unique Id for infectious (for contact tracing). It will get a unique starting value per process
static double  timeOfAction = 100000000;   // When control measures get activated (modified during the run)
static int  TOTCLASSES = 0; // Number of classes (to be set by ageNames())

std::vector<double>  SOCIALDIST_PROB; // [1..3] Current values at the different extent levels (above). [0] Current value in grid element
std::vector<double>  TRAVELREDUCTION;
std::vector<double>  STAYATHOME_AGE;
std::vector<double>  STAYATHOME_OTH;
std::vector<double>  STAYATHOME_SCH;
std::vector<double>  FAMILY_TRANSMIT;
std::vector<double>  STAYATHOME_FULL;
std::vector<double>  SCHOOL_CLOSURE;

static  std::set<int>  contactEvents;
static  std::map<int, int>  ageNames;  // from class to age group
static  std::map<int, int>  infectNames;
static  std::map<int, int>  susceptNames;
static  std::array<int,3>   occupancy;
static  std::vector<int>    cumulDeaths;
static  std::array<int, 3>  goccupancy;
static  std::vector<int>    gcumulDeaths;
enum {OCC_HOS = 0, OCC_ICU, OCC_HOM};

// Age groups for: 1. (ages) tiff files from WorldPop database; 2. (groups) coarsed age groups to be used in this sim
std::vector<int> ages   = {0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80};
std::vector<int> groups = {0, 0, 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16};

// R0 values for different transmission levels (calculated for the above values. They are *almost* independent on scaling of above values)
//static std::vector<double>  R0_tau = {1.02682, 2.1001, 4.45365, 8.53905, 16.7546};


std::vector< std::vector< std::vector<double> > > baseMap;
std::vector< std::vector<double> > popMap;
std::vector<double>  sizes;

bool firstinfected;
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
void  doFitting(int);


#define  Kenya 1
#define  Italy 2

#ifndef  COUNTRY
#error   COUNTRY is not defined.  Use '--define COUNTRY=XXX'  where  XXX is an implemented country.
#error   Currently implemented countries are: KENYA and ITALY.
#else
#	if COUNTRY == Italy
#		include "Italy-setup.cpp"
#	elif COUNTRY == Kenya
#		include "Kenya-setup.cpp"
#	else
#		error Unrecognized name for COUNTRY.
#	endif
#endif

// Joe & Matt first approach.  NOT USED now
static  std::vector<double>  susceptibility = {
	0.00282026497661079,
	0.00364014459290738,
	0.00288739193982039,
	0.0120953741300004,
	0.0377834101532493,
	0.0579847625162934,
	0.068531081222694,
	0.0672540390061978,
	0.06158432860095,
	0.0965786717540175,
	0.118036726176092,
	0.156930489217356,
	0.284027163570426,
	0.393072078509953,
	0.411985317312746,
	0.641615696240083
};


// Sam's symptomatic fractions by age and transmissibility.  
// Rows are age groups, columns are transmission levels (0.0, 0.15, 0.25, 0.50, 1.00)
static std::vector< std::vector<double> >  symptomatic = {
	{0.0032162928721740547, 0.0019086591398639221, 0.001526200293108504, 0.0015050181239767404, 0.0014432192790456596},
	{0.0025163985012532524, 0.0013458086273749326, 0.0010464701838377496, 0.000962294586193556, 0.0009357685327922587},
	{0.00334013535150579, 0.001591552396296836, 0.0012034412294856776, 0.001122214491377564, 0.0010761342759786922},
	{0.003659635031500512, 0.0018821688155611128, 0.0014635316688939558, 0.0014231888053253985, 0.0013458086273749326},
	{0.03059388453983771, 0.018366923571798966, 0.014893246507446924, 0.014686542795053286, 0.013695268824959014},
	{0.03974987470281975, 0.021419249458555244, 0.01761274351004964, 0.01736829583592139, 0.016889531485125703},
	{0.08656450365699372, 0.04430308724092091, 0.03397095366723376, 0.031238410139814973, 0.03167807081661833},
	{0.06543702388380312, 0.03397095366723376, 0.026786815880391338, 0.02568689937146875, 0.024290277337257693},
	{0.06533371691360246, 0.03397095366723376, 0.029539946993910837, 0.028326981798033054, 0.028326981798033054},
	{0.12365117334923718, 0.06552219944626281, 0.04750977837915582, 0.04248392003306334, 0.040739451216136126},
	{0.1641805544108621, 0.10391828181800267, 0.08194176439475051, 0.07430476961279056, 0.07430476961279056},
	{0.12118804724498734, 0.12462334187198994, 0.12815601614527444, 0.12815601614527444, 0.13743204603117717},
	{0.1949095327760354, 0.17674391118172914, 0.18175403736365037, 0.17923146887955166, 0.17923146887955166},
	{0.26507498195578383, 0.27642552857253, 0.3091267376127518, 0.3178894947074111, 0.3315015659109177},
	{0.12004566675626463, 0.17674391118172914, 0.237033821675666, 0.2803160413878102, 0.2923192122329743},
	{0.6050993286002224, 0.5112690095796085, 0.4145743494560396, 0.3812269649128545, 0.375935907058367},
	{1.00, 1.00, 1.00, 1.00, 1.00}
};


// This is the class storing individual data
class Infos : public IndivData {
public:
	Infos() {isolate = 0; isolate_duration = 0; athome = 0; infamily = 0; length = 0; uniqueId = 0;};
	~Infos() {};
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
	double  athome_prev = 0.0;
	std::vector<int> contacts;
	std::vector<int> family;
private:
};


Infos * indivDataFactory()  {
	Infos *data = nullptr;
	try {
		data = new Infos();
	} catch (std::bad_alloc &err)  {
		simStatus.abort("Bad Infos() allocation.");
	}
	return data;
}



void  getAgeGroups()  {
	std::ifstream handler;
	char ch_filename[256];
	std::string  filename;
	int ncols;
	int nrows;

	int nAgeGroups = groups[ groups.size()-1 ] + 1;
	cumulDeaths.resize(nAgeGroups);
	gcumulDeaths.resize(nAgeGroups);

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
		for (int qq = nrows-1; qq >= 0; qq--)  {
//		for (int qq = 0; qq < nrows; qq++)  {
			for (int kk = 0; kk < ncols; kk++)  {
				handler >> baseMap[aa][kk][qq];
			}
		}
		handler.close();
	}
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
	totalCases = 0;     // Process total number of cases
	cumulCases = 0;     // Process cumul number of cases
	gcumulCases = 0;    // Global cumul number of cases (for output)
	cumulAsyCases = 0;  // Process cumul number of asymptomatic cases
	gcumulAsyCases = 0; // Global cumul number of asymptomatic cases (for output)
	occupancy = goccupancy = {0, 0, 0};
	for (int jj = 0; jj < cumulDeaths.size(); jj++)  cumulDeaths[jj] = gcumulDeaths[jj] = 0;
	firstinfected = false;
}


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
	//int yy = simStatus.getMaxY()-simStatus.getY()-1;
	int yy = simStatus.getY();
	int nn, count;
	Infos *data;
	//std::cout << "INIT " << xx << " " << yy << " " << xpre << " " << ypre << " " << ncols << "\n" << std::flush;
	//std::cout << "XXX\n" << std::flush;
	//std::cout << "POPS1 " << popMap[xx][yy] << "\n";
	//std::cout << "YYY\n" << std::flush;
	//std::cout << "POPS2 " << simStatus.getLocalPopulationSize() << "\n";
	//std::cout << "ZZZ\n" << std::flush;
	//assert( popMap[xx][yy] == simStatus.getLocalPopulationSize() );
	data = reinterpret_cast<Infos*>( simStatus.getIndividualData( indivId, classId ) );
	if (data != nullptr)  {
		delete data;
	}
	data = indivDataFactory();
	data->length = 0;
	simStatus.setIndividualData( indivId, classId, data );

	if (xx != xpre || yy != ypre)  {
/*		probs.clear();
		probs.resize( nAgeGroups, 0 );
//		std::cout << xx << " " << yy << " " << popMap.size() << " " << popMap[xx].size() << "\n" << std::flush;
//		std::cout << xx << " " << yy << " " << popMap[xx][yy] << "\n" << std::flush;
		for (int jj = 0; jj < nAgeGroups; jj++)  {
			probs[jj] = baseMap[jj][xx][yy]*1.0/(1.0*popMap[xx][yy]);
			if (jj > 0)  probs[jj] += probs[jj-1];
//			std::cout << "&& " << (probs[jj]-1)*1000 << "\n";
		}
*/
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
	while (nn < nAgeGroups)  {
		if (occupy[nn] > 0)  {
			simStatus.moveCurrentIndividualToClass( classmapper["Sus_" + std::to_string(nn)] );
//std::cout << "Moved (" << nn << ") to " << "Sus_" + std::to_string(nn) << "\n";
			occupy[nn]--;
			count++;
			break;
		} else {
			nn++;
		}
	}
	if (nn == nAgeGroups)  {
		assert( count == simStatus.getLocalPopulationSize() );
	}
}



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
	if (tracingOn && RNG->get() < TRACING_PROB)  {
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
	std::vector<int>().swap(data->contacts);
	std::vector<int>().swap(data->family);
}



void  moveToCase()   {
	int indivId = simStatus.getCurrentIndividualId();
	int classId = simStatus.getCurrentIndividualClass();
	Infos * data = reinterpret_cast<Infos*>( simStatus.getIndividualData( indivId, classId ) );
	if (data->uniqueId > 0 && tracingOn)  {
		contactEvents.insert(data->uniqueId);
	}
	totalCases++;
	cumulCases++;
	// Clear storage space to reduce memory footprint
	recovery(data);
}



bool moveToHos()  {
	occupancy[OCC_HOS]++;
	moveToCase();
	return true;
}



bool moveToIcu()  {
	occupancy[OCC_ICU]++;
	moveToCase();
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
	cumulAsyCases++;
	return true;
}



bool  moveToInfect()  {
	//totalCases++;
	return true;
}



void  firstinfect( double case_lon, double case_lat )  {
	int  constexpr  MINAGEGROUP = 7;
	RandomGenerator *gen = simStatus.getRandomGenerator();
	int nAgeGroups = groups[ groups.size()-1 ] + 1;
	std::vector<double>  counts;
	double lat = simStatus.getLatitude();
	double lon = simStatus.getLongitude();
	int pop, kk;
	// Vo' Euganeo
	if (haversine(lon, lat, case_lon, case_lat) < 1.0*GRIDRES)  {
		counts.clear();
		pop = simStatus.getLocalPopulationSize();
		for (int jj = MINAGEGROUP; jj < nAgeGroups; jj++)  {
			counts.push_back(
					static_cast<double>( simStatus.getCountIndividuals( classmapper["Sus_"+std::to_string(jj)] ) )
			);
//			std::cout << "&& " << counts[jj-MINAGEGROUP] << "\n";
		}
		double rnd = gen->get();
		kk = MINAGEGROUP+extractFromDistribution( rnd, counts, -1, 0 );
		if (kk >= 0)  {
			if (kk > nAgeGroups)  simStatus.abort( "Age class out of bounds!" );
//			std::cout << "&&&&&& [" << kk << "] " << simStatus.getCountIndividuals( classmapper["Sus_"+std::to_string(kk)] ) << "\n";
			int classId = classmapper["Sus_"+std::to_string(kk)];
			Infos *data = reinterpret_cast<Infos*>( simStatus.getIndividualData( 0, classId ) );
			data->uniqueId = ++mainUniqueId;

			simStatus.moveGenericIndividualFromToClass( classId, classmapper["Inf_"+std::to_string(kk)] );
			firstinfected = true;
			data->infamily = 1;
			//totalCases++;
		} else {
			throw "Cannot find first infective";
		}
	}
}



void  initPolicyParams()  {
	int maxval = 0;
	for (int jj = 0; jj < policyApplication.size(); jj++)  {
		if (maxval < policyApplication[jj])  maxval = policyApplication[jj];
	}
	SOCIALDIST_PROB.assign(maxval+1, 0.0); // [1..3] Current values at the different extent levels (above). [0] Current value in grid element
	TRAVELREDUCTION.assign(maxval+1, 0.0);
	STAYATHOME_AGE.assign(maxval+1, 0.0);
	STAYATHOME_OTH.assign(maxval+1, 0.0);
	STAYATHOME_SCH.assign(maxval+1, 0.0);
	FAMILY_TRANSMIT.assign(maxval+1, 1.0);
	STAYATHOME_FULL.assign(maxval+1, 0.0);
	SCHOOL_CLOSURE.assign(maxval+1, 0.0);
}



void  addAllPolicies( PolicyQueue &queue )  {
	for (int jj = 0; jj < policyParams.size(); jj++)  {
		for (int kk = 0; kk < policyParams[jj].size(); kk++)  {
			switch (kk)  {
				case  POLICY_SOCIALDIST_PROB:
					queue.addPolicy( Policy(policyApplication[jj], &SOCIALDIST_PROB[policyApplication[jj]], &policyTime[jj], policyDuration[jj][kk], policyParams[jj][kk]) );
					break;
				case  POLICY_TRAVELREDUCTION:
					queue.addPolicy( Policy(policyApplication[jj], &TRAVELREDUCTION[policyApplication[jj]], &policyTime[jj], policyDuration[jj][kk], policyParams[jj][kk]) );
					break;
				case  POLICY_STAYATHOME_AGE:
					queue.addPolicy( Policy(policyApplication[jj], &STAYATHOME_AGE[policyApplication[jj]], &policyTime[jj], policyDuration[jj][kk], policyParams[jj][kk]) );
					break;
				case  POLICY_STAYATHOME_OTH:
					queue.addPolicy( Policy(policyApplication[jj], &STAYATHOME_OTH[policyApplication[jj]], &policyTime[jj], policyDuration[jj][kk], policyParams[jj][kk]) );
					break;
				case  POLICY_STAYATHOME_SCH:
					queue.addPolicy( Policy(policyApplication[jj], &STAYATHOME_SCH[policyApplication[jj]], &policyTime[jj], policyDuration[jj][kk], policyParams[jj][kk]) );
					break;
				case  POLICY_FAMILY_TRANSMIT:
					queue.addPolicy( Policy(policyApplication[jj], &FAMILY_TRANSMIT[policyApplication[jj]], &policyTime[jj], policyDuration[jj][kk], policyParams[jj][kk]) );
					break;
				case  POLICY_STAYATHOME_FULL:
					queue.addPolicy( Policy(policyApplication[jj], &STAYATHOME_FULL[policyApplication[jj]], &policyTime[jj], policyDuration[jj][kk], policyParams[jj][kk]) );
					break;
				case  POLICY_SCHOOL_CLOSURE:
					queue.addPolicy( Policy(policyApplication[jj], &SCHOOL_CLOSURE[policyApplication[jj]], &policyTime[jj], policyDuration[jj][kk], policyParams[jj][kk]) );
					break;
			}
		}
	}
}



void  accessCycle( int status )  {
	static bool switcher = true;
	static bool activateStayAtHome = false;
	static PolicyQueue queue;
	double sum = 0;

	switch (status)  {
		case CYCLE_INIT:
			mainUniqueId = 10000000*simStatus.getProcessId();
			simStatus.setDailyFractions( {8.0, 8.0+(45.0/7.0)} );
			initMobility();
			getAgeGroups();
			buildAgeAssignments();
			//getSusceptibility();
			initContactMatrix();
//			std::cout << "CM " << params.KK_0_10 << " " << params.KK_10_0 << "\n";
#ifdef  MODEL_FAMILY
			params.home = 0;
#endif
			updateContactMatrix();

/*			if (SCHOOL_CLOSURE)  {
				queue.addPolicy( Policy( &params.school, &timeOfAction, 0.0 ) );
			}
			queue.addPolicy( Policy( &params.other,  &timeOfAction, 21, 1.0-SOCIALDIST_PROB ) );
			queue.addPolicy( Policy( &activateStayAtHome, &timeOfAction ) );
			queue.addPolicy( Policy( &params.work,   &timeOfAction, 1.0-STAYATHOME_OTH ) );
			policyQueue.setTime( params.t0 ); // days from first imported case
*/
			break;

		case CYCLE_START:
#ifndef FITTING
			getSymptomaticRate(DET_RATE, SYMPTOMATIC_MULTIPLIER);
			//params.eigen = R0_tau[DET_RATE];
//			std::cout << "@@@ " << params.eigen <<" " << params.R0 <<" " << 1.0/params.gamma << " " << 1.0 / params.sigma << "\n";
#endif
			initPolicyParams();
			resetCounters();
			addAllPolicies( queue );
			queue.setStart( params.t0 );
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
			updateContactMatrix();
			break;

		case CYCLE_EVAL:
			evalLocalParameters();
			if (!firstinfected && simStatus.getTime() >= firstInfections[0][0]) {
				firstinfect( firstInfections[0][1], firstInfections[0][2] );
			}
			enforcePolicies();
			//if (activateStayAtHome)  {
				stayAtHome();
//				activateStayAtHome = true;
			//}
			break;

		case CYCLE_PRE:
			shareOutputData();
			contactEvents.clear();
			break;

		case CYCLE_POST:
			totalCases = 0;
#ifndef  FITTING
			if (simStatus.getTime() > 60)  {
				// gcumulXXXX have been set in printmap(), executed before CYCLE_POST
				if (gcumulCases+gcumulAsyCases < 10)  {
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
	doFitting(status);
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
	std::set<int>  tmpContactEvents;

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
	buf.pack( &totalCases, DATA_INT, DATA_ADD );
	buf.reduce();
	int fullsz = 0;
	std::vector<int> sizes( simStatus.getNumberOfProcesses(), 0 );
	for (int jj = 0; jj < simStatus.getNumberOfProcesses(); jj++)  {
		buf.unpack( &sizes[jj] );
		fullsz += sizes[jj];
	}
	buf.unpack( &totalCases );
#ifndef FITTING
	std::cout << "TotalCases: " << totalCases << "\n";
#endif

	// EARLY EXIT!
	if (totalCases > TRACING_THRESHOLD)  {
		tracingOn = false;
		contactEvents.clear();
		return;
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
	tmpContactEvents.clear();
	int  tval;
	for (int jj = 0; jj < fullsz; jj++)  {
		//buf.unpack( &tmpContactEvents[jj] );
		buf.unpack( &tval );
		tmpContactEvents.insert(tval);
//		std::cout << "TMPcontactEVENTS " << tmpContactEvents[jj] << "\n";
	}
	contactEvents.swap( tmpContactEvents );
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
				if (age <= 3)  {
					rate_athome = &STAYATHOME_SCH[0];
				} else if (age >= 14)  {
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
				if (RNG->get() < TRACING_PROB)  {
					data->isolate = 1+ISOLATE_GAPDAYS;
					data->isolate_duration = 14;
					data->contacts.clear();
//std::cout<<"ISOLATING !\n";
					break;
				}
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
	if (prob > 0)  {
		if (RNG->get() >= params.eta)  return false;
		if (RNG->get() >= params.mobility)  return false;

		int indivId = simStatus.getCurrentIndividualId();
		int classId = simStatus.getCurrentIndividualClass();
		Infos *data = reinterpret_cast<Infos*>( simStatus.getIndividualData( indivId, classId ) );
		if (data->athome == 1)  {
			return false;
		}
		if (data->length == 0)  {
			data->length = getMobilityDuration(rnd, dist);
			if (data->length == 0)  return false;
		}
		retval = !handleIsolation(data, classId);
		contactEvents.clear();

	} else {
/*		int indivId = simStatus.getCurrentIndividualId();
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
*/
retval = true;
	}
	return retval;
}



void  buildAgeAssignments()  {
	int nAgeGroups = groups[ groups.size()-1 ] + 1;
	std::vector<std::string>  names = {"Sus_", "Esp_", "Inf_", "Asy_", "Hos_", "Icu_", "Hom_", "Rec_", "Ded_", "Asyrec_"};
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

	names = {"Inf_", "Asy_", "Hom_"};
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
		sz = RNG->poisson( 4.0 );
	}
	std::vector<double> pp;
	int count = 1, notfound = 0;
	auto  info = simStatus.getIndividualInfo( indivId, classId );
	assert( info.placeOfContact == 0);

	familymembers.push_back( {info.home.rank, indivId, classId} );
	for (int jj = 0; jj < nAgeGroups; jj++)  {
		pp.push_back( localKK_home[ ageNames[classId] ][ jj ] );
		if (jj > 0)  pp[jj] += pp[jj-1];
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
				familymembers.push_back( {otherinfo.home.rank, otherId, otherClass} );
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
	double  prob = params.R0 * params.gamma * FAMILY_TRANSMIT[0] / params.eigen;
//	double  prob = params.R0 * params.gamma / params.eigen;
	infosthis = simStatus.getIndividualInfo( indivId, classId );
	
	if (infosthis.placeOfContact == 0 && isInfective(classId))  {
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
				if (RNG->get() < 1.0-std::exp( -prob ))  {
					Infos *other = reinterpret_cast<Infos*>( simStatus.getIndividualData( infosother.id, infosother.status ) );
					if (data->uniqueId == 0)  {
						data->uniqueId = ++mainUniqueId;
					}
			
					if (tracingOn && RNG->get() < TRACING_PROB)  {
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


/*
	if (infosthis.placeOfContact == 0 && isSusceptible(classId))  {
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
			if (isInfective(infosother.status))  {
				if (RNG->get() < 1.0-std::exp( -prob ))  {
					Infos *other = reinterpret_cast<Infos*>( simStatus.getIndividualData( infosother.id, infosother.status ) );
					if (other->uniqueId == 0)  {
						other->uniqueId = ++mainUniqueId;
					}
			
					if (tracingOn && RNG->get() < TRACING_PROB)  {
						data->contacts.push_back(other->uniqueId);
					}

					int newClass = ageNames[ classId ];
					newClass = classmapper[ "Esp_" + std::to_string(newClass) ];

//					std::cout << "MOVING: " << indivId << "/" << simStatus.cell->indivhome->size() << " " << classId << " " << newClass << "\n" << std::flush;
					simStatus.moveSpecificIndividualToClass( indivId, classId, newClass );
				}
			}
		}
	}*/
}



void  getSymptomaticRate(int detect, double multiplier)  {
	params.zz_0 = symptomatic[0][detect]*multiplier;
	params.zz_1 = symptomatic[1][detect]*multiplier;
	params.zz_2 = symptomatic[2][detect]*multiplier;
	params.zz_3 = symptomatic[3][detect]*multiplier;
	params.zz_4 = symptomatic[4][detect]*multiplier;
	params.zz_5 = symptomatic[5][detect]*multiplier;
	params.zz_6 = symptomatic[6][detect]*multiplier;
	params.zz_7 = symptomatic[7][detect]*multiplier;
	params.zz_8 = symptomatic[8][detect]*multiplier;
	params.zz_9 = symptomatic[9][detect]*multiplier;
	params.zz_10 = symptomatic[10][detect]*multiplier;
	params.zz_11 = symptomatic[11][detect]*multiplier;
	params.zz_12 = symptomatic[12][detect]*multiplier;
	params.zz_13 = symptomatic[13][detect]*multiplier;
	params.zz_14 = symptomatic[14][detect]*multiplier;
	params.zz_15 = symptomatic[15][detect]*multiplier;
	params.zz_16 = symptomatic[16][detect]*multiplier;
//	std::cout << "ZZ " << params.zz_7 << " " << params.zz_15 << "\n";
}


void  getSusceptibility()  {
	params.zz_0 = susceptibility[0];
	params.zz_1 = susceptibility[1];
	params.zz_2 = susceptibility[2];
	params.zz_3 = susceptibility[3];
	params.zz_4 = susceptibility[4];
	params.zz_5 = susceptibility[5];
	params.zz_6 = susceptibility[6];
	params.zz_7 = susceptibility[7];
	params.zz_8 = susceptibility[8];
	params.zz_9 = susceptibility[9];
	params.zz_10 = susceptibility[10];
	params.zz_11 = susceptibility[11];
	params.zz_12 = susceptibility[12];
	params.zz_13 = susceptibility[13];
	params.zz_14 = susceptibility[14];
	params.zz_15 = susceptibility[15];
//	std::cout << "ZZ " << params.zz_7 << " " << params.zz_15 << "\n";
}



void  shareOutputData()  {
	int nAgeGroups = groups[ groups.size()-1 ] + 1;

	DataBuffer  data;
	data.setBuffer( sizeof(int)*(5+cumulDeaths.size()), (5+cumulDeaths.size()) );	
	data.clear();
	data.pack( &cumulCases, DATA_INT, DATA_ADD );
	data.pack( &cumulAsyCases, DATA_INT, DATA_ADD );
	for (int kk = 0; kk < 3; kk++)  {
		data.pack( &occupancy[kk], DATA_INT, DATA_ADD );
	}
	for (int kk = 0; kk < nAgeGroups; kk++)  {
		data.pack( &cumulDeaths[kk], DATA_INT, DATA_ADD );
	}
	data.reduce();
	data.unpack( &gcumulCases );
	data.unpack( &gcumulAsyCases );
	for (int kk = 0; kk < 3; kk++)  {
		data.unpack( &goccupancy[kk] );
	}
	for (int kk = 0; kk < nAgeGroups; kk++)  {
		data.unpack( &gcumulDeaths[kk] );
	}
}



void  printextra( FILE *handler, int fd )  {
	int nAgeGroups = groups[ groups.size()-1 ] + 1;

	
	fprintf( handler, "\t%d %d\t%d %d %d %g %g\t", gcumulCases, gcumulAsyCases, goccupancy[OCC_HOS], goccupancy[OCC_ICU], goccupancy[OCC_HOM], SOCIALDIST_PROB[0], SCHOOL_CLOSURE[0] );
	int sum = 0;
	for (int kk = 0; kk < nAgeGroups; kk++)  {
		sum += gcumulDeaths[kk];
		fprintf( handler, "%d ", gcumulDeaths[kk] );
	}
	fprintf( handler, "%d", sum );
}



void  printmapextra( FILE *handler )  {
	int nAgeGroups = groups[ groups.size()-1 ] + 1;
	int hosp = 0, icus = 0;

	for (int jj = 0; jj < nAgeGroups; jj++)  {
		hosp += simStatus.getCountIndividuals( classmapper[ "Hos_" + std::to_string(jj) ] );
		icus += simStatus.getCountIndividuals( classmapper[ "Icu_" + std::to_string(jj) ] );
	}
	fprintf( handler, "\t%d %d", hosp, icus );
}





int constexpr STOREBUFFER_SZ = 100*1000*1000;

void  doFitting( int status )  {
	enum {PARAM_R0 = 0, PARAM_TAU, PARAM_ZMAX, PARAM_T0};
	static FittingABC  fitting( "Adaptive" );
	static Particle    trialParticle;
	static RandomGenerator *deryaRNG;
	static std::vector<std::string> paramNames;
	static std::vector< std::array<int,2> >  timeseries;
	static int cases, deaths;
	static unsigned int tsduration;
	static DataBuffer  chiSqShare, evtData, paramShare;
	static std::ifstream  in;
	static std::ofstream  outputFile, outputFileAvg;
	static std::string    dirPrefix = "Generation-";
	static bool  hasrestarted, isgood;
	static std::array<std::vector<int>, 2>  dailyCount;
	double fzero = 0.0;
	static double chiSquared;
	static int    countEvents, countCases;
	static double ctime;
	static std::string outRunName = "successfulRun-" + to_string( simStatus.getProcessId() ) + ".dat";
	enum {DATA_CASES = 0, DATA_DEATHS};

	switch( status )  {
		case CYCLE_INIT:
			in.open( timeseriesFile );
			if (in.fail())  {
				simStatus.abort( "Could not open file [timeseries.dat]" );
			}
			do {
				in >> cases >> deaths;
				timeseries.push_back({cases, deaths});
			} while (in.good());
			in.close();
			tsduration = timeseries.size();
			std::cout << "Loaded data for " << tsduration << " days.\n";

			deryaRNG = simStatus.getRandomGenerator();
			fitting.setRNG( deryaRNG );
			fitting.setTrials( 1 );  // Trials per particle
			fitting.setDistrs( 1 ); // Number of intermediate distributions after which we assume convergence
			fitting.setParticlesPerTrial( 2 ); // REPLICAS?
			fitting.setDiscardRatio( 0.10, 0.30 );
			fitting.setThreshold( 2.5 );
			fitting.setKernelAmplitude( 1.0 );

			fitting.addParameter( "R0",		0.0,	PARTYPE_UNIFORM,	{ 1., 7.} );
			fitting.addParameter( "gamma",	0.0,	PARTYPE_UNIFORM,	{ 1., 8.} );
			fitting.addParameter( "tau",	0.0,	PARTYPE_UNIFORM,	{ 0., 5.} );
			fitting.addParameter( "zmax",	0.0,	PARTYPE_UNIFORM,	{ 0.2, 1.0} );
			fitting.addParameter( "t0",		0.0,	PARTYPE_UNIFORM,	{ 1., 60.} );
//			fitting.addParameter( "t0",		0.0,	PARTYPE_UNIFORM,	{ 1., 10.} );
			fitting.addParameter( "eta",	0.0,	PARTYPE_UNIFORM,	{ 0.4, 1.0} );
			fitting.addParameter( "theta",	0.0,	PARTYPE_UNIFORM,	{ 0.25, 4.0} );
			fitting.init( NN );
			paramNames = {"R0", "gamma", "tau", "zmax", "t0", "eta", "theta"};

			//paramShare.setCommType( DATABUFFER_COMMGLOBAL );
			//chiSqShare.setCommType( DATABUFFER_COMMGLOBAL );
			//evtData.setCommType( DATABUFFER_COMMGLOBAL );
			paramShare.setBuffer( 6*sizeof(double), 6 );
			chiSqShare.setBuffer( 2*sizeof(double)+1*sizeof(int), 2+1 );
			evtData.setBuffer( 2*sizeof(int), 2 );

			if (simStatus.getGlobalProcessId() == 0 && params.Restart == 0)  {
				std::system( ("\\rm -rv " + dirPrefix + "*").c_str() );
				// After evaluate, populationTimer has increased -> create dirs for storing future particles
				std::string  dirName = dirPrefix + std::to_string(fitting.getPopulationTimer());
				mkdir( dirName.c_str(), 0755 );
			}

			if (simStatus.getGlobalProcessId() == 0)  {
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
			}
			break;

		case CYCLE_START:
//			simStatus.setSimulationLength( nweeks * 7 + 1 );
			if (simStatus.getTime() == 0.0 && params.Restart)  {
				std::cout << "Loading data from storage.\n";
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
				for (int jj = 0; jj < simStatus.getNumberOfProcesses(); jj++)  {
					os.read( data.get(), size );
					if (simStatus.getProcessId() == jj)  {
						sz = deryaRNG->loadState( data.get() );
						if (sz != size)  {
							simStatus.abort("Error in reading stored RNG data.");
						}
					}
				}

				os.read( reinterpret_cast<char*>(&size), sizeof(int) );
				// THIS SHOULD BE USED AS CHECK. MODIFY
				os.read( buffer.get(), size );
				sz = fitting.loadState( buffer.get() );
				if (size != sz)  {
					simStatus.abort("Error in reading state data.");
				}
				if (size > STOREBUFFER_SZ)  {
					simStatus.abort("Size of stored data too large!");
				}
				fb.close();

//				firstTime = true;
//				lastTime = -1;
//				initialized = false;
//				simStatus.setTime( 0 );

//				simStatus.setTime(0);
//				simStatus.reinitDisease();

				if (params.Restart == 2)  fitting.reset();
				params.Restart = 0;

			} else if (simStatus.getTime() == 0.0 && hasrestarted) {
				std::cout << "Saving data to storage.\n";
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
					size = fitting.saveState( buffer.get() );
					if (size > STOREBUFFER_SZ)  exit(1);
					os.write( reinterpret_cast<char*>(&size), sizeof(int) );
					os.write( buffer.get(), size );

					fb.close();
				}
			}
			hasrestarted = false;

			dailyCount[0].clear();
			dailyCount[1].clear();
			trialParticle = fitting.start();
			params.R0 = trialParticle.parameters[0];
			params.gamma = 1.0/trialParticle.parameters[1];
			params.tau   = trialParticle.parameters[2];
			params.zmax  = trialParticle.parameters[3];
			params.t0    = trialParticle.parameters[4];
			params.eta   = trialParticle.parameters[5];
			params.theta = trialParticle.parameters[6];
			// Share parameters among nodes
/*			paramShare.clear();
			if (simStatus.getProcessId() == 0)  {
				std::cout << "Params: " << params.R0 << " " << params.gamma << " " << params.tau << " " << params.zmax << " " << params.t0 << " "  << params.eta << " -> " << params.t0 + timeseries.size() << "\n";
				paramShare.pack( &params.R0,    DATA_DOUBLE, DATA_ADD );
				paramShare.pack( &params.gamma, DATA_DOUBLE, DATA_ADD );
				paramShare.pack( &params.tau,   DATA_DOUBLE, DATA_ADD );
				paramShare.pack( &params.zmax,  DATA_DOUBLE, DATA_ADD );
				paramShare.pack( &params.t0,    DATA_DOUBLE, DATA_ADD );
				paramShare.pack( &params.eta,   DATA_DOUBLE, DATA_ADD );
			} else {
				paramShare.pack( &fzero, DATA_DOUBLE, DATA_ADD );
				paramShare.pack( &fzero, DATA_DOUBLE, DATA_ADD );
				paramShare.pack( &fzero, DATA_DOUBLE, DATA_ADD );
				paramShare.pack( &fzero, DATA_DOUBLE, DATA_ADD );
				paramShare.pack( &fzero, DATA_DOUBLE, DATA_ADD );
				paramShare.pack( &fzero, DATA_DOUBLE, DATA_ADD );
			}
			paramShare.reduce();
			paramShare.unpack( &params.R0 );
			paramShare.unpack( &params.gamma );
			paramShare.unpack( &params.tau );
			paramShare.unpack( &params.zmax );
			paramShare.unpack( &params.t0 );
			paramShare.unpack( &params.eta );
*/
			if (simStatus.getProcessId() == 0)  {
				std::cout << "[" << simStatus.__fullProcessId << "] Params: R0 " << params.R0 << ", gamma " << params.gamma << ", tau " << params.tau;
				std::cout << ", zmax " << params.zmax << ", t0 " << params.t0 << ", eta "  << params.eta << ", theta " << params.theta << " -> " << params.t0 + timeseries.size() << "\n";
			}
			fitting.acceptParticle(trialParticle);

			params.beta = params.R0 * params.gamma;
			getSymptomaticRate( static_cast<int>(params.tau), params.zmax );
//			params.eigen = R0_tau[static_cast<int>(params.tau)];

			chiSquared = 0.0;
			simStatus.setSimulationLength( params.t0 + timeseries.size() );
			countEvents = 0;
			countCases  = 0;
			break;

		case CYCLE_RESET:
			break;

		case CYCLE_EVAL:
			break;

		case CYCLE_PRE:
			if (simStatus.isOutputlineFrame(0))  {
				countEvents = simStatus.getCountEvents( 0, eventmapper["Deaths"] );
				countCases  = simStatus.getCountEvents( 0, eventmapper["CaseHos"] ) + simStatus.getCountEvents( 0, eventmapper["CaseIcu"] ) + simStatus.getCountEvents( 0, eventmapper["CaseHom"] );
				evtData.clear();
				evtData.pack( &countEvents,  DATA_INT, DATA_ADD );
				evtData.pack( &countCases,   DATA_INT, DATA_ADD );
				evtData.reduce();
				evtData.unpack( &countEvents );
				evtData.unpack( &countCases );

				dailyCount[0].push_back( countEvents );
				dailyCount[1].push_back( countCases );
			}
			break;

		case CYCLE_POST:
			ctime = simStatus.getTime();
			if (simStatus.isOutputlineFrame(0))  {
				if (ctime > params.t0)  {
					double numerator = 0.0, denominator = 0.0;
					double chiSquaredTmp = 0.0;
					int day = static_cast<int>(ctime-params.t0+0.001);

					if (day < timeseries.size())  {
						double tmpnum = countCases;
						double tmpden = countCases;
						tmpnum += -timeseries[ day ][0];
						tmpden +=  timeseries[ day ][0];
						numerator   += tmpnum*tmpnum;
						denominator += tmpden*tmpden;

						chiSquaredTmp = denominator > 0 ? (numerator) / sqrt(denominator) : 0.0;
						chiSquared += chiSquaredTmp;
						//hasrestarted = fitting.step( chiSquared );
						//if (hasrestarted)  {
						//	if (simStatus.getProcessId() == 0)  {
						//		std::cout << "Restarting in [All Processes]\n";
						//	}
						//}
					}
				}
				countEvents = 0;
				countCases  = 0;
			}
			break;

		case CYCLE_LAST:
//std::cout << "[" << simStatus.getGlobalProcessId() << "] - CYCLE_LAST start \n";
			if (chiSquared < fitting.getError())  {
				std::string  dirName = dirPrefix + std::to_string(fitting.getPopulationTimer());
				std::ofstream  outRun = std::ofstream( dirName + "/" + outRunName, ios::app );
				for (unsigned idx = 0; idx < dailyCount[1].size(); idx++)  {
					outRun << idx;
					outRun << "\t" << dailyCount[1].at(idx);
					outRun << "\n";
				}
				outRun << "\n";
				outRun.close();
			}
//std::cout << "[" << simStatus.getGlobalProcessId() << "] - CYCLE_LAST CHECK 0 \n";
			isgood = fitting.evaluate( chiSquared, hasrestarted );
//std::cout << "[" << simStatus.getGlobalProcessId() << "] - CYCLE_LAST CHECK 1 \n";
			if (simStatus.getGlobalProcessId() == 0)  {
				int accepted, total;
				double  ratio = fitting.getCurrentAcceptanceRatio( &accepted, &total );
				std::cout << "Current ratio: [" << (accepted) << "/" << (total) << "]";
				std::cout << " -> [" << (100.0*(ratio)) << "%] ";
				ratio = fitting.getAcceptanceRatio( &accepted, &total );
				std::cout << "  (LAST WAS) [" << (accepted) << "/" << (total) << "]";
				std::cout << " -> [" << (100.0*(ratio)) << "%]\n";

			}

//std::cout << "[" << simStatus.getGlobalProcessId() << "] - CYCLE_LAST CHECK 2 \n";
			if (hasrestarted)  {
//std::cout << "[" << simStatus.getGlobalProcessId() << "] - CYCLE_LAST CHECK 3 \n";
				int  timeIndex;

				if (simStatus.getGlobalProcessId() == 0)  {
					// After evaluate, populationTimer has increased -> create dirs for storing future particles
					std::string  dirName = dirPrefix + std::to_string(fitting.getPopulationTimer());
					mkdir( dirName.c_str(), 0755 );

					std::vector<double> avg = fitting.getAverages();
					std::vector<double> var = fitting.getVariances();
					outputFileAvg << avg[0];
					for (unsigned int idx = 1; idx < avg.size(); idx++)  {
						outputFileAvg << "\t" << avg[idx] << "\t" << var[idx];
					}
					outputFileAvg << "\t" << fitting.getError() << "\n" << std::flush;
					timeIndex = static_cast<int>(floor(avg[0] + 0.0001));

					std::vector< std::vector<double> >  validpars = fitting.getAcceptedParameters();
					for (unsigned int idx = 0; idx < validpars.size(); idx++)  {
						outputFile << fitting.getPopulationTimer();
						for (unsigned int kk = 0; kk < validpars[idx].size(); kk++)  {
							//if (kk > 0)  outputFile << "\t";
							outputFile << "\t" << validpars[idx][kk];
						}
						outputFile << "\n" << std::flush;
					}

				}

				if (simStatus.getGlobalProcessId() == 0)  {
					std::vector<double> avg = fitting.getAverages();
					std::vector<double> var = fitting.getVariances();
					std::cout << "VALS " << avg[0];
					for (unsigned int idx = 1; idx < avg.size(); idx++)  {
						std::cout << "\t" << avg[idx] << "\t" << var[idx];
					}
					std::cout << "\t" << fitting.getError() << "\n" << std::flush;
				}

			}

			if (fitting.getPopulationSamplingTimer() <= TT)  {
				simStatus.setTime( 0 );
				simStatus.reinitDisease();
			}
			break;

		case CYCLE_FINALIZE:
			if (simStatus.getProcessId() == 0)  {
				outputFileAvg.close();
				outputFile.close();
			}
			fitting.finalize();
			break;

	}
}




void  evaluateR0( int tau_lvl )  {
	int nAgeGroups = groups[ groups.size()-1 ] + 1;

	std::vector< std::vector<double> >  MM;
	MM.resize(3*nAgeGroups);;
	for (int jj = 0; jj < nAgeGroups; jj++)  {
		MM[jj].resize(3*nAgeGroups, 0.0);
	}

	// Build FF matrix:
	
}

