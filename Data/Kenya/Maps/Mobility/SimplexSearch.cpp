/* 
 * Program: nmsimplex.c
 * Author : Michael F. Hutt
 * http://www.mikehutt.com
 * Nov. 3, 1997
 *
 * An implementation of the Nelder-Mead simplex method.
 *
 * Copyright (c) 1997 Michael F. Hutt
 * Released under the MIT License
 *
 */

#include "SimplexSearch.h"

/* create the initial simplex */
void SimplexSearch::initialize_simplex(std::vector<double> &start, double scale)  {
	double pn,qn;   /* values used to create initial simplex */

	//pn = scale*(std::sqrt(nn+1.0)-1.0+nn)/(nn*std::sqrt(2.0));
	//qn = scale*(std::sqrt(nn+1.0)-1.0)/(nn*std::sqrt(2.0));

	pn = 0.5*scale;
	qn = 1.0*scale;

//	pn = 0.5*(std::sqrt(5.0)-1.0)*scale;
//	qn = 1.0*scale;

	for (int ii=0; ii<nn; ii++) {
		vv[0][ii] = start[ii];
	}
	
	for (int ii=1; ii <= nn; ii++)  {
		for (int jj=0; jj < nn; jj++) {
			if (ii-1 == jj) {
				vv[ii][jj] = (1+pn) * start[jj];
			} else  {
				vv[ii][jj] = (1+qn) * start[jj];
			}
		}
	}
}



/* print out the initial values */
void SimplexSearch::print_initial_simplex()  {
	printf("Initial Values\n");
	for (int jj=0; jj <= nn; jj++) {
		for (int ii=0; ii<nn; ii++) {
			std::cout << vv[jj][ii] << ", ";
		}
		std::cout << "value " << ff[jj] << "\n";
	}
}


/* print out the value at each iteration */
void SimplexSearch::print_iteration(int itr)  {
	std::cout << "Iteration " << itr << "\n";
	for (int jj=0; jj <= nn; jj++) {
		for (int ii=0; ii < nn; ii++) {
			std::cout << vv[jj][ii] << " " << ff[jj] << "\n";
		}
	}
}


/* find the index of the largest value */
int SimplexSearch::vg_index(int vg)  {
	for (int jj=0; jj <= nn; jj++) {
		if (ff[jj] > ff[vg]) {
			vg = jj;
		}
	}
	return vg;
}


/* find the index of the smallest value */
int SimplexSearch::vs_index(int vs)  {
	for (int jj=0; jj<=nn; jj++) {
		if (ff[jj] < ff[vs]) {
			vs = jj;
		}
	}
	return vs;
}


/* find the index of the second largest value */
int SimplexSearch::vh_index(int vh, int vg)  {
	for (int jj=0; jj <= nn; jj++) {
		if (ff[jj] > ff[vh] && ff[jj] < ff[vg]) {
			vh = jj;
		}
	}
	return vh;
}


/* calculate the centroid */
void SimplexSearch::centroid(std::vector<double> &vm, int vg)  {
	double cent;

	for (int jj=0; jj<nn;jj++) {
		cent=0.0;
		for (int mm=0; mm<=nn; mm++) {
			if (mm != vg) {
				cent += vv[mm][jj];
			}
		}
		vm[jj] = cent/nn;
	}
}



SimplexSearch::SimplexSearch(int v_nn)  {nn = v_nn;}

SimplexSearch::~SimplexSearch()  {}


void   SimplexSearch::evalStats()  {
	double tmp;
	mean_vv.assign(nn, 0);
	error_vv.assign(nn, 0);
	mean_ff = 0;
	error_ff = 0;

	for (int jj=0; jj<nn;jj++) {
		for (int ii=0; ii<=nn; ii++) {
			mean_vv[jj] += vv[ii][jj];
		}
		mean_vv[jj] /= (nn+1);
	}
	for (int jj=0; jj<nn;jj++) {
		for (int ii=0; ii<=nn; ii++) {
			tmp = (vv[ii][jj] - mean_vv[jj]);
			error_vv[jj] += tmp*tmp;
		}
		error_vv[jj] = 1.96*std::sqrt( error_vv[jj]/(nn) )/sqrt(nn+1);
	}

	for (int ii=0; ii<=nn; ii++) {
		mean_ff += ff[ii];
	}
	mean_ff /= (nn+1);
	for (int ii=0; ii<=nn; ii++) {
		tmp += (ff[ii]-mean_ff);
		error_ff += tmp*tmp;
	}
	error_ff = 1.96*std::sqrt(error_ff/nn)/std::sqrt(nn+1);
	error_rel = error_ff / std::fabs(mean_ff);
}



void   SimplexSearch::printStats(int itr)  {
	std::cout << "========================================================================\n";
	std::cout << "Current stats for iteration [" << itr << "]:\n";
	for (int ii = 0; ii < nn; ii++)  {
		std::cout << "  Param " << ii << ": [" << mean_vv[ii] << " \u00b1 " << error_vv[ii] << "]\n";
	}
	std::cout << "\n";
	std::cout << "  Evaluation at point: [" << mean_ff << " \u00b1 " << error_ff << "]\n";
	std::cout << "  Relative error is:   [ " << error_rel << "]\n";
	std::cout << "------------------------------------------------------------------------\n";
	std::cout << "Value of points:\n";
	for (int ii=0; ii <= nn; ii++) {
		for (int jj=0; jj < nn; jj++) {
			std::cout << "Point " << ii << ", Param " << jj << ": [" << vv[ii][jj] << "]";
			if (jj < nn-1)  {
				std::cout << "\n";
			} else {
				std::cout << " -> [" << ff[ii] << "]\n";
			}
		}
	}
	std::cout << "========================================================================\n";
}


std::vector<double>  SimplexSearch::getParamsMean()  {
	return mean_vv;
}


std::vector<double>  SimplexSearch::getParamsError()  {
	return error_vv;
}



double SimplexSearch::search(double (*objfunc)(const std::vector<double>&), std::vector<double> &start, double EPSILON, double scale, void (*constrain)(std::vector<double> &), bool doprint)  {

	int vs;         /* vertex with smallest value */
	int vh;         /* vertex with next smallest value */
	int vg;         /* vertex with largest value */

	int ii,jj,row;
	int kk;   	  /* track the number of function evaluations */
	int itr;	  /* track the number of iterations */

	double fr;      /* value of function at reflection point */
	double fe;      /* value of function at expansion point */
	double fc;      /* value of function at contraction point */
	std::vector<double> vr;     /* reflection - coordinates */
	std::vector<double> ve;     /* expansion - coordinates */
	std::vector<double> vc;     /* contraction - coordinates */
	std::vector<double> vm;     /* centroid - coordinates */
	double min;

	double fsum,favg,ss;

	/* dynamically allocate arrays */

	/* allocate the rows of the arrays */
	vv.resize(nn+1);
	ff.resize(nn+1);
	vr.resize(nn);
	ve.resize(nn);
	vc.resize(nn);
	vm.resize(nn);

	/* allocate the columns of the arrays */
	for (int ii=0; ii<=nn; ii++) {
		vv[ii].resize(nn);
	}

	/* create the initial simplex */
	initialize_simplex(start,scale);
	print_initial_simplex();

	/* impose constraints */
	if (constrain != nullptr) {
		for (int jj=0; jj<=nn; jj++) {
			constrain(vv[jj]);
		}
	}

	/* find the initial function values */
	for (int jj=0; jj<=nn; jj++) {
		ff[jj] = objfunc(vv[jj]);
	}
	kk = nn+1;

	/* print out the initial values */
	print_initial_simplex();

	/* begin the main loop of the minimization */
	for (itr=1; itr<=MAX_IT; itr++) {
		/* find the index of the largest value */
		vg = vg_index(0);

		/* find the index of the smallest value */
		vs = vs_index(0);
		
		/* find the index of the second largest value */
		vh = vh_index(vs,vg);
		
		/* calculate the centroid */
		centroid(vm,vg);

		/* reflect vg to new vertex vr */
		for (jj=0; jj<=nn-1; jj++) {
			/*vr[j] = (1+ALPHA)*vm[j] - ALPHA*v[vg][j];*/
			vr[jj] = vm[jj]+ALPHA*(vm[jj]-vv[vg][jj]);
		}
		if (constrain != nullptr) {
			constrain(vr);
		}
		fr = objfunc(vr);
		kk++;
		
		if (fr < ff[vh] && fr >= ff[vs]) {
			for (jj=0; jj<=nn-1; jj++) {
				vv[vg][jj] = vr[jj];
			}
			ff[vg] = fr;
		}

		/* investigate a step further in this direction */
		if (fr < ff[vs]) {
			for (jj=0; jj<=nn-1; jj++) {
				/*ve[j] = GAMMA*vr[j] + (1-GAMMA)*vm[j];*/
				ve[jj] = vm[jj]+GAMMA*(vr[jj]-vm[jj]);
			}
			if (constrain != nullptr) {
				constrain(ve);
			}
			fe = objfunc(ve);
			kk++;
			
      /* 
	 by making fe < fr as opposed to fe < f[vs], 			   
	 Rosenbrocks function takes 63 iterations as opposed 
	 to 64 when using double variables. 
      */
			
			if (fe < fr) {
				for (jj=0; jj<=nn-1; jj++) {
					vv[vg][jj] = ve[jj];
				}
				ff[vg] = fe;
			} else {
				for (jj=0; jj<=nn-1; jj++) {
					vv[vg][jj] = vr[jj];
				}
				ff[vg] = fr;
			}
		}

		/* check to see if a contraction is necessary */
		if (fr >= ff[vh]) {
			if (fr < ff[vg] && fr >= ff[vh]) {
				/* perform outside contraction */
				for (jj=0; jj<=nn-1; jj++) {
					/*vc[j] = BETA*v[vg][j] + (1-BETA)*vm[j];*/
					vc[jj] = vm[jj]+BETA*(vr[jj]-vm[jj]);
				}
				if (constrain != nullptr) {
					constrain(vc);
				}
				fc = objfunc(vc);
				kk++;
			} else  {
				/* perform inside contraction */
				for (jj=0; jj<=nn-1; jj++) {
					/*vc[j] = BETA*v[vg][j] + (1-BETA)*vm[j];*/
					vc[jj] = vm[jj]-BETA*(vm[jj]-vv[vg][jj]);
				}
				if (constrain != nullptr) {
					constrain(vc);
				}
				fc = objfunc(vc);
				kk++;
			}


			if (fc < ff[vg])  {
				for (jj=0; jj<= nn-1; jj++) {
					vv[vg][jj] = vc[jj];
				}
				ff[vg] = fc;
			} else {
	/* 
	   at this point the contraction is not successful,
	   we must halve the distance from vs to all the 
	   vertices of the simplex and then continue.
	   1997-10-31 - modified to account for ALL vertices. 
	*/
	
				for (row=0; row<=nn; row++) {
					if (row != vs) {
						for (jj=0; jj<=nn-1; jj++) {
							vv[row][jj] = vv[vs][jj]+(vv[row][jj]-vv[vs][jj])/2.0;
						}
					}
				}

				/* re-evaluate all the vertices */
				for (jj=0; jj<=nn; jj++) {
					ff[jj] = objfunc(vv[jj]);
				}
	
				/* find the index of the largest value */
				vg = vg_index(0);
	
				/* find the index of the smallest value */
				vs = vs_index(0);
	
				/* find the index of the second largest value */
				vh = vh_index(vs,vg);

				if (constrain != nullptr) {
					constrain(vv[vg]);
				}
				ff[vg] = objfunc(vv[vg]);
				kk++;
				if (constrain != nullptr) {
					constrain(vv[vh]);
				}
				ff[vh] = objfunc(vv[vh]);
				kk++;
			}
		}

		/* print out the value at each iteration */
		//print_iteration(itr);
		evalStats();
		if (doprint)  {
			printStats(itr);
		}

		/* test for convergence */
		if (error_rel < EPSILON) break;
	}
	/* end main loop of the minimization */

	/* find the index of the smallest value */
	vs = vs_index(0);
	
	/*printf("The minimum was found at\n"); */
	for (jj=0; jj<nn; jj++) {
		/*printf("%e\n",v[vs][j]);*/
		start[jj] = vv[vs][jj];
	}
	min=objfunc(vv[vs]);
	kk++;
  	std::cout << kk << " Function Evaluations\n";
	std::cout << itr <<" Iterations through program\n";

	return min;
}

