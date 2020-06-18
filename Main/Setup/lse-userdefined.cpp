#include <fstream>
#include "GeoTiffReader.h"

#ifdef  REDUCE_FACTOR
#	warning Defined REDUCE_FACTOR
#endif


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


std::vector< std::vector< std::vector<unsigned int> > >  agegrMap;
std::vector< std::vector<unsigned int> >  baseMap, tmpMap;
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
	int delta, xsize, ysize;

	RandomGenerator *deryaRNG = simStatus.getRandomGenerator();

	for (int qq = 0; qq < ages.size(); qq++)  {
		for (int kk = 0; kk < sex.size(); kk++)  {
			sprintf( ch_filename, ("../Gridded/" + std::string(SHORTCOUNTRY) + "_%c_%d_"+ YEAR +".tif").c_str(), sex[kk], ages[qq] );
			filename = std::string( ch_filename );
			std::cout << "Handling file [" << filename << "]\n" << std::flush;
			GeoTiffReader reader( filename );
			data = reader.readHeader();
			if (qq == 0 && kk == 0)  {
				std::cout << "Map dimensions are " << data.ncols << " x " << data.nrows << " grid elements.\n";
				tmpMap.resize( data.ncols );
				for (int qq = 0; qq < data.ncols; qq++)  {
					tmpMap[qq].resize( data.nrows, 0 );
				}

				delta = 10*GRIDRES;
				xsize = data.ncols / delta + (data.ncols % delta == 0 ? 0 : 1);
				ysize = data.nrows / delta + (data.nrows % delta == 0 ? 0 : 1);

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
				}

				baseMap.resize( xsize );
				for (int xx = 0; xx < xsize; xx++)  {
					baseMap[xx].resize( ysize, 0 );
				}
			}

			reader.readFile( 1.0, tmpMap );

			aa = groups[qq];
			for (int ii = 0; ii < data.ncols; ii++)  {
				for (int jj = 0; jj < data.nrows; jj++)  {
					int xx = ii/delta;
					int yy = jj/delta;
					tmpMap[ii][jj] = deryaRNG->binomial( tmpMap[ii][jj], REDUCE_FACTOR );
					baseMap[xx][yy] += tmpMap[ii][jj];
					agegrMap[aa][xx][yy] += tmpMap[ii][jj];
					ageSizes[aa] += tmpMap[ii][jj];
//std::cout << ii << " " << jj << " " << baseMap[aa][ii][jj] << " " << tmpMap[ii][jj] << " " << REDUCE_FACTOR << "\n";
				}
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

	filename = std::string( COUNTRY ) + "_" + std::to_string(GRIDRES) + "km.asc";
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
			handler << baseMap[xx][yy];
		}
		handler << "\n";
	}
	handler.close();

	for (int aa = 0; aa < nAgeGroups; aa++)  {
		filename = std::string( COUNTRY ) + "_" + std::to_string(GRIDRES) + std::string("km_") + std::to_string(aa) + ".dat";
		std::ofstream handler = std::ofstream(filename.c_str());
		for (int yy = ysize-1; yy >= 0; yy--)  {
			for (int xx = 0; xx < xsize; xx++)  {
				if (xx > 0)  handler << " ";
				if (baseMap[xx][yy] > 0)  {
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
		filename = std::string( COUNTRY ) + "_" + std::to_string(GRIDRES) + std::string("km_") + std::to_string(aa) + ".asc";
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
				if (baseMap[xx][yy] > 0)  {
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
		filename = std::string( COUNTRY ) + "_" + std::to_string(GRIDRES) + "km_stats.dat";
		std::ofstream handler = std::ofstream(filename.c_str());
		for (int aa = 0; aa < nAgeGroups; aa++)  {
			handler << 5*aa << " " << 5*aa+4 << " " << ageSizes[aa] << "\n";
		}
		handler.close();
	}

}


