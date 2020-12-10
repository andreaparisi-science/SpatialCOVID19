#include <assert.h>
#include <algorithm>
#include <fstream>
#include <cmath>

#include "SimplexSearch.h"


void assignPreferredLocations();
void storeIndividualData( const char * );

enum {METHOD_CPC = 1, METHOD_DISTANCE = 2, METHOD_NORMDIST = 3};
constexpr int METHOD = METHOD_CPC;
constexpr double NR = 0.2;
constexpr double minLim = 90.0;
constexpr double delta  = 0.25;
constexpr double TIMEFRAME = 1.0;
constexpr double REL_ERROR_X = 1.e-3;
constexpr double REL_ERROR_Y = 1.e-4;
const std::string commutingFile = "Data/Commuting.dat";
const std::string country = "Spain";

class MobData {
public:
	int  key;
	unsigned int size;
	int nlocs;
	std::vector<int>  xx, yy;
	std::vector<double>  cumulProb;
};

class GroupData {
public:
	GroupData()  {populationSize = 0; centroid[0] = centroid[1] = 0.0; cellsCount = 0;}
	std::vector<double> commutingData;
	std::vector<double> commutingModel;
	double populationSize;
	double centroid[2];
	int    cellsCount;
};

void  loadAscFile(const std::string&, std::vector< std::vector<MobData> > &);
void	  loadStoreFile(const std::string&, std::vector< std::vector<MobData> > &);
void  loadCommutingData(const std::string&, std::vector<GroupData> &, int);
double  evaluateModel(double factor, int method);
std::map<int, int>  buildGroups( const std::string );

std::vector< std::vector<MobData> >  popMap;
std::map<int, int>  groups;
std::string strgroups;
std::string strgroups_East_England  = "34+33+106+105+101+99+102+108+107+103+110+100+104+109,135+136+131+137+128+130+134+132+133+129,32+55+56,61+65+64+63+62+31,181+182+178+177+180+176+179,224+227+228+222+225+223+226";
std::string strgroups_East_Midlands = "77+75+79+76+72+78+74+73+15,203+199+197+200+202+201+198+18,171+172+174+173+169+170+175,163+166+164+168+162+165+167+16,17,188+187+184+189+186+183+185";
std::string strgroups_Greater_London = "326+313+306+325+315+321+323+305+312+300+298+302+311+320+314+317+322+301+299+316+304+297+309+295+319+318+324+307+303+296+308+310+294";
std::string strgroups_NorthEast_England = "48,278+277+279+280+281,47+5+1,4+3+2";
std::string strgroups_NorthWest_England = "49+50+6+7,67+71+69+66+70+68,258+259+260+261+262+263+264+265+266+267,160+151+159+152+156+161+154+157+155+150+158+153+9+8,268+269+270+271+272";
std::string strgroups_SouthEast_England = "37+38+41+36+40+39,59+58+60+57+42,95+97+98+94+96+43,120+121+127+123+118+122+125+117+126+119+124+45+44,46,140+142+144+148+149+143+146+138+145+139+141+147+35,205+204+206+207+208,235+234+236+239+229+231+238+232+230+233+237,251+246+247+249+248+250+245";
std::string strgroups_SouthWest_England = "22+24+211+213+212+210+209,23,25+114+116+111+112+115+113,30+54,93+92+90+91+89+28+29+88,81+80+82+83+86+87+84+85+27+26,53+52";
std::string strgroups_West_Midlands = "19,51+20,214+215+216+217+218+219+220+221+21,240+241+242+243+244,282+283+284+285+286+287+288,252+253+254+255+256+257";
std::string strgroups_Yorkshire = "276+275+273+274,293+291+290+289+292,196+192+190+193+191+194+195+14,11+10,13+12";
std::string strgroups_Scotland = "354+355+327+334+330+337+342+350+333+341+346+353,357+352+328+329+339+358+336+331+351+356+344+345,347+348,335+338+349+332+340+343";
std::string strgroups_Wales = "359+360+361+362+365+367+366+380+373+375+374+376+370+369+368,377+378+372+371+363+364+379";
std::vector<GroupData> region;
MapData  data;
int nGroups = 0;
DataBuffer bufferData;

inline double  ToRad( double angle )  {
	return angle * M_PI / 180.0;
}


double  ToAngleX( double x0, MapData &data )  {
	return (x0-0.5)*data.cellsize+data.xOrig;
}



double  ToAngleY( double y0, MapData &data )  {
//	std::cout << "TOANGLE-Y: " << y0 << " " << data.cellsize << " " << data.yOrig << "\n";
	return (y0-0.5)*data.cellsize+data.yOrig;
}


double  haversine( double, double, double, double );
/*double  haversine( double lon1, double lat1, double lon2, double lat2 )  {
	static constexpr double RR = 6371; // Km
	double dLat2 = ToRad(lat2-lat1)/2;
	double dLon2 = ToRad(lon2-lon1)/2;
	lat1 = ToRad(lat1);
	lat2 = ToRad(lat2);
	double tmp1 = std::sin(dLat2) * std::sin(dLat2);
	double tmp2 = std::sin(dLon2) * std::sin(dLon2);
	double aa = tmp1 + tmp2 * std::cos(lat1) * std::cos(lat2);
	//c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
	double cc = 2 * std::asin(std::sqrt(aa));
	return  RR * cc;
}*/



int  insertLogHistogram(double xx, double yy, double x0, double delta, std::vector<double> &histo)  {
	int nn = static_cast<int>( (std::log(xx)-std::log(x0))/delta );
	if (nn >= histo.size())  {
		histo.resize(nn+1);
	}
	histo[nn] += yy;
	return nn;
}






void  loadAscFile(const std::string &filename, std::vector< std::vector<MobData> > &map )  {
	std::ifstream  handler;
	double fval;
	std::string  var;
	int  val;
	handler.open(filename);
	if (handler.eof() || handler.fail())  {
		std::cerr << "Could not read header of file [" << filename << "]\n";
		exit(1);
	}

	for (int loop = 0; loop < 6; loop++)  {
		handler >> var;
		std::cout << "Header: [" << var << "] ";
		if (var == "ncols")  {
			handler >> val;
			data.ncols = val;
			std::cout << "[" << val << "]\n";
		} else if (var == "nrows")  {
			handler >> val;
			data.nrows = val;
			std::cout << "[" << val << "]\n";
		} else if (var == "xllcorner")  {
			handler >> fval;
			data.xOrig = fval;
			std::cout << "[" << fval << "]\n";
		} else if (var == "yllcorner")  {
			handler >> fval;
			data.yOrig = fval;
			printf("[%g]\n", fval );
		} else if (var == "cellsize")  {
			handler >> fval;
			data.cellsize  = fval;
			printf("[%g]\n", fval );
		} else if (var == "NODATA_value")  {
			handler >> data.nodata;
			std::cout << "[" << data.nodata << "]\n";
		} else {
			std::cerr << "Unknown header [" << var << "]\n";
			exit(1);
		}
	}
	std::cout << std::flush;

	map.resize( data.ncols );
	for (int ii = 0; ii < data.ncols; ii++)  {
		map[ii].resize( data.nrows );
	}
	for (int jj = 0; jj < data.nrows; jj++)  {
		for (int ii = 0; ii < data.ncols; ii++)  {
			handler >> map[ii][jj].key;
//std::cout << ii << "\t" << jj << "\t" << map[ii][jj].key << "\n";
		}
	}
	handler.close();
}



void  loadStoreFile(const std::string& filename, std::vector< std::vector<MobData> > &populationMap)  {
	std::ifstream  handler;
	int xs, ys, xx, yy, x0, y0, nlocs, nelems;
	double  cumulProb, rf;

	handler.open(filename, std::ios::in & std::ios::binary);
	if (handler.eof() || handler.fail())  {
		std::cerr << "Unable to open file [" << filename << "].\n";
		exit(1);
	}
	handler.read( reinterpret_cast<char*>(&rf), sizeof(double) );
	handler.read( reinterpret_cast<char*>(&xs), sizeof(int) );
	handler.read( reinterpret_cast<char*>(&ys), sizeof(int) );
   	assert( xs == data.ncols && ys == data.nrows );

	while( handler.read( reinterpret_cast<char*>(&nelems), sizeof(int) ) )  {
		for (int idx = 0; idx < nelems; idx++)  {
			handler.read( reinterpret_cast<char*>(&x0), sizeof(int) );
//printf( "[%d] CHECK 3.1: x0 [%d] - y0 [%d] - sz [%d]\n", 0, x0, y0, populationMap[x0][y0].size ); fflush(stdout);
			handler.read( reinterpret_cast<char*>(&y0), sizeof(int) );
//printf( "[%d] CHECK 3.1: x0 [%d] - y0 [%d] - sz [%d]\n", 0, x0, y0, populationMap[x0][y0].size ); fflush(stdout);
			handler.read( reinterpret_cast<char*>(&populationMap[x0][y0].size), sizeof(unsigned int) );
//printf( "[%d] CHECK 3.1: x0 [%d] - y0 [%d] - sz [%d]\n", 0, x0, y0, populationMap[x0][y0].size ); fflush(stdout);
//			for (int nn = 0; nn < populationMap[x0][y0].size; nn++)  {
				handler.read( reinterpret_cast<char*>(&nlocs), sizeof(int) );
				populationMap[x0][y0].nlocs = nlocs;
//printf( "[%d] CHECK 3.1: x0 [%d] - y0 [%d] - sz [%d] - nlocs [%d]\n", 0, x0, y0, populationMap[x0][y0].size, nlocs ); fflush(stdout);
//printf( "[%d] CHECK 3.2\n", simStatus.__processId ); fflush(stdout);
				for (int loc = 0; loc < nlocs; loc++)  {
//printf( "[%d] CHECK 3.2.1\n", simStatus.__processId ); fflush(stdout);
					handler.read( reinterpret_cast<char*>(&xx), sizeof(int) );
					handler.read( reinterpret_cast<char*>(&yy), sizeof(int) );
					handler.read( reinterpret_cast<char*>(&cumulProb), sizeof(double) );
					populationMap[x0][y0].xx.push_back( xx );
					populationMap[x0][y0].yy.push_back( yy );
					populationMap[x0][y0].cumulProb.push_back( cumulProb );
//printf( "[%d] CHECK 3.2.2\n", simStatus.__processId ); fflush(stdout);
//printf( "[%d] CHECK 3.2.3\n", simStatus.__processId ); fflush(stdout);
				}
//			}
		}
	}
	handler.close();
}


// "1+2,3+4,5,6,7+9,8,10"
std::map<int,int>  buildGroups( const std::string strgroups )  {
	std::vector<std::string>  elems;
	std::vector<int>  list;
	std::map<int, int>  groups;
	{
		auto beg = strgroups.begin();
		auto it  = strgroups.begin();
		while (it != strgroups.end())  {
			it = std::find_if( beg, strgroups.end(), [](int ch) {return ch == ',';});
			elems.push_back( std::string( beg, it ) );
//			std::cout << elems.back() << "\n";
			beg = it+1;
		}
	}
	for (int ii = 0; ii < elems.size(); ii++)  {
		auto beg = elems[ii].begin();
		auto it  = elems[ii].begin();
		while (it != elems[ii].end())  {
			it = std::find_if( beg, elems[ii].end(), [](int ch) {return ch == '+';});
			groups[ std::stoi( std::string( beg, it ) ) ] = ii+1;
			beg = it+1;
		}
	}

//	for (auto&& elem: groups)  {
//		std::cout << "Group " << elem.first <<" " << elem.second << "\n";
//	}
	return groups;
}





void  loadCommutingData(const std::string& commutingFile, std::vector<GroupData> &region, int nGroups)  {
	std::ifstream  handler;
	handler.open(commutingFile, std::ios::in);
	if (handler.eof() || handler.fail())  {
		std::cerr << "Unable to open file [" << commutingFile << "].\n";
		exit(1);
	}
	for (int ii = 0; ii <= nGroups; ii++)  {
		for (int jj = 0; jj <= nGroups; jj++)  {
			region[ii].commutingData[jj] = 0;
		}
	}
	for (int ii = 1; ii <= nGroups; ii++)  {
		for (int jj = 1; jj <= nGroups; jj++)  {
			handler >> region[ii].commutingData[jj];
		}
	}
	handler.close();
}




void  resetModel()  {
	for (int x0 = 0; x0 < data.ncols; x0++)  {
		for (int y0 = 0; y0 < data.nrows; y0++)  {
			popMap[x0][y0].xx.clear();
			popMap[x0][y0].yy.clear();
			popMap[x0][y0].cumulProb.clear();
			popMap[x0][y0].nlocs = 0;
		}
	}

	for (int ii = 0; ii <= nGroups; ii++)  {
		for (int jj = 0; jj <= nGroups; jj++)  {
			region[ii].commutingModel[jj] = 0;
		}
		region[ii].populationSize = 0;
		region[ii].centroid[0] = 0;
		region[ii].centroid[1] = 0;
		region[ii].cellsCount = 0;
//		region[ii].populationSize = 0;
	}
}




void  updateModelData()  {
	int keyFrom, keyTo, xx, yy;
	double prob;

	for (int x0 = 0; x0 < data.ncols; x0++)  {
		for (int y0 = 0; y0 < data.nrows; y0++)  {
			if (popMap[x0][y0].key == 0)  continue;

			keyFrom = groups[ popMap[x0][y0].key ];
//if (keyFrom == 1)  {
//std::cout << "BEFORE: [" << (region[keyFrom].centroid[0]/region[keyFrom].cellsCount) << "]-[" << (region[keyFrom].centroid[1]/region[keyFrom].cellsCount) << "\n";
//std::cout << "  " << x0 << " " << y0 << " " << data.cellsize << " " << data.yOrig << " " << ToAngleY(y0, data) << "\n";
//}
			region[keyFrom].centroid[0] += ToAngleX(x0, data);
			region[keyFrom].centroid[1] += ToAngleY(y0, data);
			region[keyFrom].cellsCount++;
//if (keyFrom == 1)  {
//std::cout << "AFTER: [" << (region[keyFrom].centroid[0]/region[keyFrom].cellsCount) << "]-[" << (region[keyFrom].centroid[1]/region[keyFrom].cellsCount) << "\n";
//}
			for (int kk = 0; kk < popMap[x0][y0].nlocs; kk++)  {
				xx = popMap[x0][y0].xx[kk];
				yy = popMap[x0][y0].yy[kk];
				prob = popMap[x0][y0].cumulProb[kk];
				if (kk > 0)  prob -= popMap[x0][y0].cumulProb[kk-1];
				if (popMap[xx][yy].key == 0)  continue;
				keyTo   = groups[ popMap[xx][yy].key ];
				if ( keyFrom != keyTo )  {
					region[keyFrom].commutingModel[keyTo] += prob*popMap[x0][y0].size;
				}
//std::cout << "** " << keyFrom <<" " << keyTo << " " << popMap[x0][y0].size<< " " << prob << "\n";
			}
			region[keyFrom].populationSize += popMap[x0][y0].size;
		}
	}
	int sum = 0;
	for (int ii = 1; ii <= nGroups; ii++)  {
//		std::cout << "# " << ii << " " << region[ii].populationSize << "\n";
		sum += region[ii].populationSize;
		region[ii].centroid[0] /= region[ii].cellsCount;
		region[ii].centroid[1] /= region[ii].cellsCount;
//		std::cout << "@@ " << region[ii].centroid[0] << " " << region[ii].centroid[1] << "\n";
	}
//	std::cout << "## " << sum << "\n";
}



double  evaluateModel(double factor, int method)  {
	double  tmp, chiSq = 0.0;
	double  NCC = 0;
	double  NCM = 0;
	double  NCD = 0;
	for (int keyFrom = 1; keyFrom <= nGroups; keyFrom++)  {
		for (int keyTo = 1; keyTo <= nGroups; keyTo++)  {
			if (keyFrom == keyTo)  continue;

			if (method == METHOD_DISTANCE)  {
//std::cout << "]] " <<commutingModel[keyFrom][keyTo] << " " << commutingData[keyFrom][keyTo] << "\n";
				tmp    = region[keyFrom].commutingModel[keyTo] - (region[keyFrom].commutingData[keyTo]*region[keyFrom].populationSize*factor);
				//chiSq += tmp*tmp/region[keyFrom].populationSize;
				chiSq += tmp*tmp;
			} else if (method == METHOD_NORMDIST)  {
				tmp    = region[keyFrom].commutingModel[keyTo]*1.0/region[keyFrom].populationSize - (region[keyFrom].commutingData[keyTo]*factor);
				//chiSq += tmp*tmp/region[keyFrom].populationSize;
				chiSq += tmp*tmp;
			} else  {
				NCC += std::min( region[keyFrom].commutingModel[keyTo], region[keyFrom].commutingData[keyTo]*region[keyFrom].populationSize*factor );
				NCM += region[keyFrom].commutingModel[keyTo];
				NCD += region[keyFrom].commutingData[keyTo]*region[keyFrom].populationSize*factor;
//				std::cout << "[" << keyFrom << "][" << keyTo << "] " << NCC <<" " << NCM << " " << NCD << "\n";
			}
		}
		//std::cout << "\n";
	}
	if (method == METHOD_DISTANCE || method == METHOD_NORMDIST)  {
		return std::sqrt(chiSq);
	} else {
		return -2*NCC/(NCM+NCD);
	}
} 



void  buildHistograms(const arma::vec &nr, double delta)  {
	std::vector<double> dataHistogram, simsHistogram;
	double dist, dataLim = 0, simsLim = 0, maxLim;
	std::ofstream handler = std::ofstream( "pdf.dat" );
	for (int keyFrom = 1; keyFrom <= nGroups; keyFrom++)  {
		for (int keyTo = 1; keyTo <= nGroups; keyTo++)  {
			if (keyFrom == keyTo)  continue;
			handler << region[keyFrom].commutingData[keyTo]*region[keyFrom].populationSize*(1.0/TIMEFRAME) << "\t" << region[keyFrom].commutingModel[keyTo] << "\n";
			dist = haversine( region[keyFrom].centroid[0], region[keyFrom].centroid[1], region[keyTo].centroid[0], region[keyTo].centroid[1] );
			if (dist < minLim)  continue;
			dataLim = std::max( dataLim, static_cast<double>(insertLogHistogram( dist, region[keyFrom].commutingData[keyTo]*region[keyFrom].populationSize*(1.0/TIMEFRAME), minLim, delta, dataHistogram )) );
			simsLim = std::max( simsLim, static_cast<double>(insertLogHistogram( dist, region[keyFrom].commutingModel[keyTo], minLim, delta, simsHistogram )) );
		}
	}
	handler.close();
	maxLim = std::max(dataLim, simsLim);
	dataHistogram.resize(maxLim, 0);
	simsHistogram.resize(maxLim, 0);
//std::cout << maxLim/delta << "\n";

	std::ofstream fout("histoData.dat");
	for (int ii = 0; ii < maxLim; ii++)  {
		fout << minLim*std::exp(ii*delta)*0.5*(1+std::exp(delta)) << "\t" << dataHistogram[ii]/(minLim*std::exp(ii*delta)*(std::exp(delta)-1)) << "\n";
	}
	fout.close();


	fout.open("histoSims.dat");
	for (int ii = 0; ii < maxLim; ii++)  {
		fout << minLim*std::exp(ii*delta)*0.5*(1+std::exp(delta)) << "\t" << simsHistogram[ii]/(minLim*std::exp(ii*delta)*(std::exp(delta)-1)) << "\n";
	}
	fout.close();

	{
		static std::map<double,double>  distances;
		double sum = 0.0;
		for (int ii = 0; ii < maxLim; ii++)  {
			sum += (dataHistogram[ii] - simsHistogram[ii])*(dataHistogram[ii] - simsHistogram[ii]);
		}
		distances[nr(0)] = std::sqrt(sum);
		std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++\n";	
		std::cout << "Histogram distances are:\n";
		for (auto &el: distances)  {
			std::cout << "Nr = [" << el.first << "] -> [" << el.second << "]\n";
		}
		std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++\n";
	}
}



double  buildMobility(const arma::vec &nr)  {
	double eval, cpc;
	double dummy = 0.0;
	bool   caughtException = false;

	std::cout << "Checking parameters [" << nr(0);
	for (int jj = 1; jj < nr.size(); jj++)  {
		std::cout << ", " << nr(jj);
	}
	std::cout << "]\n";

	for (int jj = 0; jj < nr.size(); jj++)  {
		simStatus.setMobilityParameter(jj, nr(jj));
		//simStatus.setMobilityParameter(0, nr[0]);
		//simStatus.setMobilityParameter(1, nr[1]);
		//simStatus.setMobilityParameter(4, nr[2]);
	}
	try  {
		assignPreferredLocations();
	} catch(MobilityException &ee)  {
		caughtException = true;
	}
	bufferData.clear();
	bufferData.pack(&caughtException, DATA_BOOL, DATA_OR);
	bufferData.reduce();
	bufferData.unpack(&caughtException);
	if (caughtException)  {
		std::cout << "### Caught exception in process [" << simStatus.getProcessId() << "] ###\n";
		eval = 1.e300;
	} else {
		storeIndividualData( "storage" );
		if (simStatus.getProcessId() == 0)  {
			resetModel();
			loadStoreFile("storage", popMap);
			updateModelData();
			eval = evaluateModel( (1.0/TIMEFRAME), METHOD_DISTANCE );
			cpc  = evaluateModel( (1.0/TIMEFRAME), METHOD_CPC );
			std::cout << "CPC " << cpc << "\n";
			buildHistograms(nr, delta);
		}
		bufferData.clear();
		if (simStatus.getProcessId() == 0)  {
			bufferData.pack(&eval, DATA_DOUBLE, DATA_ADD);
			bufferData.pack(&cpc,  DATA_DOUBLE, DATA_ADD);
		} else {
			bufferData.pack(&dummy, DATA_DOUBLE, DATA_ADD);
			bufferData.pack(&dummy, DATA_DOUBLE, DATA_ADD);
		}
		bufferData.reduce();
		bufferData.unpack(&eval);
		bufferData.unpack(&cpc);
	}
	std::cout << "Evaluated value for parameters [" << nr(0);
	for (int jj = 1; jj < nr.size(); jj++)  {
		std::cout << ", " << nr(jj);
	}
	std::cout << "] is [" << eval << "]; ";
	std::cout << "CPC is [" << cpc << "]\n";
	if (METHOD == METHOD_DISTANCE || METHOD == METHOD_NORMDIST)
		return  eval;
	else
		return  cpc;
}



double  constrainParams(arma::vec &vv, bool &isConstrained)  {
	isConstrained = false;
	for (int jj = 0; jj < vv.size(); jj++)  {
//		if (vec[jj] < 0)  vec[jj] = -vec[jj];
		if (vv(jj) < 0)  {
			isConstrained = true;
			return  1.e300;
		//if (vec[jj] > 1)  vec[kk] = 2.0 - vec[jj];
		}
	}
	return 0.0;
}






arma::vec  __convertVectorToVec( const std::vector<double> &vv )  {
	arma::vec yy( vv.size() );
	for (int jj = 0; jj < vv.size(); jj++)  yy(jj) = vv[jj];
	return yy;
}



void  accessCycle( int status )  {
	static  double phi = 0.5*(1.0+std::sqrt(5.));
	double  optimalVal, optimalPar, sum;
	double  nrmin, nrmax, nr1, nr2, f1, f2, nropt, dd, err = 1.0;
	bool    skip1 = false, skip2 = false;
	std::ofstream  handler, fout;
	std::vector<double> start;
	//double  minLim = 20, dataLim = 0, simsLim = 0, maxLim = 0;
	double  dist;
	std::vector<double> dataHistogram, simsHistogram;
	SimplexSearch searcher(1);
	arma::vec  vstart, pt, scale;


	switch (status)  {
		case CYCLE_INIT:
			strgroups = strgroups_East_England + "," + strgroups_East_Midlands + "," + strgroups_Greater_London + "," + strgroups_NorthEast_England + "," + strgroups_NorthWest_England +
							"," + strgroups_SouthEast_England + "," + strgroups_SouthWest_England + "," + strgroups_West_Midlands + "," + strgroups_Yorkshire +
							"," + strgroups_Scotland + "," + strgroups_Wales;
			//analyzeTripDistribution();
			loadAscFile( std::string("../") + country + "_" + std::to_string(GRIDRES) + "km_ids.asc", popMap);
			loadStoreFile("storage", popMap);
			groups = buildGroups(strgroups);

			for (auto&& elem: groups)  {
				nGroups = std::max( nGroups, elem.second );
			}
			//nGroups++;

			region.resize(nGroups+1);
			for (int ii = 0; ii < nGroups+1; ii++)  {
				region[ii].commutingData.resize(nGroups+1, 0);
				region[ii].commutingModel.resize(nGroups+1, 0);
			}

			loadCommutingData( commutingFile, region, nGroups);
			bufferData.setBuffer( 2*sizeof(double), 2 );
			break;


		case CYCLE_START:
			nrmin = 0.0;
			nrmax = 1.0;
			nropt = 1.e300;
			start.resize(1);
			start[0] = simStatus.getMobilityParameter(0);
			vstart = __convertVectorToVec( start );
			//start[0] = simStatus.getMobilityParameter(0);
			//start[1] = simStatus.getMobilityParameter(2);
			//start[2] = simStatus.getMobilityParameter(3);
			//start[1] = simStatus.getMobilityParameter(4);
			scale = {1.0};
			std::tie<arma::vec, double>(pt, f1) = searcher.search(buildMobility, vstart, scale, {REL_ERROR_X, REL_ERROR_Y}, constrainParams, simStatus.getProcessId() == 0);
//				if (simStatus.getProcessId() == 0)  {
//					std::cout << "*** Current optimal for mobility is [" << start[0];
//					for (int jj = 1; jj < start.size(); jj++)  {
//						std::cout << ", " << start[jj];
//					}
//					std::cout << "]\n";
//					std::cout << "*** Current error is [" << err << "], higher than the target [" << 1.e-3 << "]\n";
//				}

//			buildMobility(searcher.getParamsMean());
//			buildHistograms(searcher.getParamsMean(), delta);
//			}

//			std::cout << "*** Optimal parameter for mobility: [" << nropt << "] with error [" << err << "]\n";
			simStatus.exit("Done.");
			break;
	
		case CYCLE_RESET:
			break;
		case CYCLE_EVAL:
			break;
		case CYCLE_PRE:
			simStatus.preventOutputmapFrame();
			simStatus.preventOutputlineFrame(0);
			break;
		case CYCLE_POST:
			break;
		case CYCLE_LAST:
			break;
		case CYCLE_FINALIZE:
			break;
	}
}


void  printmap( FILE *handler )  {
	int count;
	RandomGenerator *deryaRNG = simStatus.getRandomGenerator();
	count = deryaRNG->binomial( simStatus.getCountEvents( eventmapper[ "Case" ] ), 0.05 );
	fprintf( handler, " %d", count );
}





