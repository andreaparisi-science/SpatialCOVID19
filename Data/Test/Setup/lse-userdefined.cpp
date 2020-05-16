#include <fstream>
#include "GeoTiffReader.h"

void  accessCycle( int status )  {
	switch (status)  {
		case CYCLE_INIT:
			transformAgeGroupData();
			break;
		case CYCLE_START:
			simStatus.exit("Early exit");
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


std::vector< std::vector< std::vector<unsigned int> > >  baseMap, agegrMap;
std::vector< std::vector<unsigned int> >  tmpMap;
void  transformAgeGroupData()  {
	std::string  filename;
	char ch_filename[256];
	std::vector<int> ages   = {0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80};
	std::vector<int> groups = {0, 0, 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16};
	std::vector<char>  sex  = {'m','f'};
	int nAgeGroups = groups[ groups.size()-1 ] + 1;
	MapData  data;
	std::vector<int>  ageSizes;
	int aa, extrastep = 0;

	RandomGenerator *deryaRNG = simStatus.getRandomGenerator();

	for (int qq = 0; qq < ages.size(); qq++)  {
		for (int kk = 0; kk < sex.size(); kk++)  {
			sprintf( ch_filename, "../Gridded/ken_%c_%d_2018.tif", sex[kk], ages[qq] );
			filename = std::string( ch_filename );
			std::cout << "Handling file [" << filename << "]\n" << std::flush;
			GeoTiffReader reader( filename );
			data = reader.readHeader();
			if (qq == 0 && kk == 0)  {
				baseMap.resize( nAgeGroups );
				for (int aa = 0; aa < nAgeGroups; aa++)  {
					baseMap[aa].resize( data.ncols );
					for (int qq = 0; qq < data.ncols; qq++)  {
						baseMap[aa][qq].resize( data.nrows, 0 );
					}
				}
			}

			tmpMap.clear();
			tmpMap.resize( data.ncols );
			for (int ii = 0; ii < data.ncols; ii++)  {
				tmpMap[ii].resize( data.nrows, 0 );
			}
			reader.readFile( 1.0, tmpMap );

			aa = groups[qq];
			for (int ii = 0; ii < data.ncols; ii++)  {
				for (int jj = 0; jj < data.nrows; jj++)  {
					baseMap[aa][ii][jj] += deryaRNG->binomial( tmpMap[ii][jj], 47564296.0/51576655.0 );
				}
			}
		}
	}
	int delta = 10*GRIDRES;
	int xsize = data.ncols / delta + (data.ncols % delta == 0 ? 0 : 1);
	int ysize = data.nrows / delta + (data.nrows % delta == 0 ? 0 : 1);

	if (xsize == 1 || ysize == 1)  {
		xsize += 2;
		ysize += 2;
		data.xOrig -= delta*data.cellsize;
		data.yOrig -= delta*data.cellsize;
		extrastep = 1;
	}

	ageSizes.resize( nAgeGroups, 0 );
	agegrMap.resize( nAgeGroups );
	for (int aa = 0; aa < nAgeGroups; aa++)  {
		agegrMap[aa].resize( xsize );
		for (int xx = 0; xx < xsize; xx++)  {
			agegrMap[aa][xx].resize( ysize, 0 );
		}
		for (int xx = 0; xx < data.ncols; xx++)  {
			for (int yy = 0; yy < data.nrows; yy++)  {
				agegrMap[aa][xx/delta][yy/delta] += baseMap[aa][xx][yy];
				ageSizes[aa] += baseMap[aa][xx][yy];
			}
		}
	}

	// If new rows/cols are needed, shift everything
	if (extrastep)  {
		for (int aa = 0; aa < nAgeGroups; aa++)  {
			for (int ii = xsize-1; ii > 0; ii--)  {
				for (int jj = ysize-1; jj > 0; jj--)  {
					agegrMap[aa][ii][jj] = agegrMap[aa][ii-1][jj-1];
				}
			}
			for (int ii = 0; ii < xsize; ii++)  {
				agegrMap[aa][ii][0] = 0;
			}
			for (int jj = 0; jj < ysize; jj++)  {
				agegrMap[aa][0][jj] = 0;
			}
		}
	}

	tmpMap.clear();
	tmpMap.resize( xsize );
	for (int xx = 0; xx < xsize; xx++)  {
		tmpMap[xx].resize( ysize, 0 );
		for (int yy = 0; yy < ysize; yy++)  {
			for (int aa = 0; aa < nAgeGroups; aa++)  {
				tmpMap[xx][yy] += agegrMap[aa][xx][yy];
			}
		}
	}


	filename = std::string("Test_") + std::to_string(GRIDRES) + "km.asc";
	std::ofstream handler = std::ofstream(filename.c_str());
	handler << "ncols\t" << xsize << "\n";
	handler << "nrows\t" << ysize << "\n";
	handler << "xllcorner\t" << data.xOrig << "\n";
	handler << "yllcorner\t" << data.yOrig << "\n";
	handler << "cellsize\t" << (delta*data.cellsize) << "\n";
	handler << "NODATA_value\t" << data.nodata << "\n";
	for (int yy = ysize-1; yy >= 0; yy--)  {
		for (int xx = 0; xx < xsize; xx++)  {
			if (xx > 0)  handler << " ";
			handler << tmpMap[xx][yy];
		}
		handler << "\n";
	}
	handler.close();

	for (int aa = 0; aa < nAgeGroups; aa++)  {
		filename = std::string("Test_") + std::to_string(GRIDRES) + std::string("km_") + std::to_string(aa) + ".dat";
		std::ofstream handler = std::ofstream(filename.c_str());
		for (int yy = ysize-1; yy >= 0; yy--)  {
			for (int xx = 0; xx < xsize; xx++)  {
				if (xx > 0)  handler << " ";
				if (tmpMap[xx][yy] > 0)  {
					handler << agegrMap[aa][xx][yy];
				} else {
					handler << 0;
				}
			}
			handler << "\n";
		}
		handler.close();
	}

	for (int aa = 0; aa < nAgeGroups; aa++)  {
		filename = std::string("Test_") + std::to_string(GRIDRES) + std::string("km_") + std::to_string(aa) + ".asc";
		std::ofstream handler = std::ofstream(filename.c_str());
		handler << "ncols\t" << xsize << "\n";
		handler << "nrows\t" << ysize << "\n";
		handler << "xllcorner\t" << data.xOrig << "\n";
		handler << "yllcorner\t" << data.yOrig << "\n";
		handler << "cellsize\t" << (delta*data.cellsize) << "\n";
		handler << "NODATA_value\t" << data.nodata << "\n";
		for (int yy = ysize-1; yy >= 0; yy--)  {
			for (int xx = 0; xx < xsize; xx++)  {
				if (xx > 0)  handler << " ";
				if (tmpMap[xx][yy] > 0)  {
					handler << agegrMap[aa][xx][yy];
				} else {
					handler << 0;
				}
			}
			handler << "\n";
		}
		handler.close();
	}

	{
		filename = std::string("Test_") + std::to_string(GRIDRES) + "km_stats.dat";
		std::ofstream handler = std::ofstream(filename.c_str());
		for (int aa = 0; aa < nAgeGroups; aa++)  {
			handler << 5*aa << " " << 5*aa+4 << " " << ageSizes[aa] << "\n";
		}
		handler.close();
	}

}


