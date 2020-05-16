/*
 * header for nmsimplex.c
 * Author : Michael F. Hutt
 * http://www.mikehutt.com
 * 
 * An implementation of the Nelder-Mead simplex method.
 *
 * Copyright (c) 1997 Michael F. Hutt
 * Released under the MIT License
 *
 */

#ifndef NM_SIMPLEX_H
#define NM_SIMPLEX_H

#include <cstdio>
#include <iostream>
#include <cmath>
#include <vector>





class  SimplexSearch  {
public:
	SimplexSearch(int v_nn);
	~SimplexSearch();
	double search(double (*objfunc)(const std::vector<double> &), std::vector<double> &start, double, double, void (*constrain)(std::vector<double> &), bool = true);
	std::vector<double>  getParamsMean();
	std::vector<double>  getParamsError();
private:
//	void my_constraints(std::vector<double> &);
	void initialize_simplex(std::vector<double> &, double scale);
	void print_initial_simplex();
	void print_iteration(int);
	int vg_index(int);
	int vs_index(int);
	int vh_index(int, int);
	void centroid(std::vector<double> &, int);
	void evalStats();
	void printStats(int);
	static constexpr int    MAX_IT = 1000;      /* maximum number of iterations */
	static constexpr double ALPHA  = 1.0;       /* reflection coefficient */
//	static constexpr double BETA   = 0.5;       /* contraction coefficient */
//	static constexpr double GAMMA  = 2.0;       /* expansion coefficient */
	static constexpr double BETA   = 0.5*(std::sqrt(5.0)-1.0);       /* contraction coefficient */
	static constexpr double GAMMA  = 0.5*(std::sqrt(5.0)+1.0);       /* expansion coefficient */
	std::vector< std::vector<double> > vv;     /* holds vertices of simplex */
	std::vector<double> ff;      /* value of function at each vertex */
	std::vector<double> mean_vv;
	std::vector<double> error_vv;
	double mean_ff;
	double error_ff;
	double error_rel;
	int nn;
};


#endif

