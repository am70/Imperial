#pragma once
#ifndef lModH_H  
#define lModH_H

#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <string>
#include<fstream>
#include <numeric>
#include <boost/regex.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/binomial_distribution.hpp>
#include<boost/random/mersenne_twister.hpp>
#include<boost/random/poisson_distribution.hpp>
#include<boost/math/distributions/gamma.hpp>
#include<boost/math/special_functions/binomial.hpp>
#include <boost/range/numeric.hpp>
#include <functional>
#include <boost/math/distributions/inverse_gaussian.hpp>
#include <boost/math/distributions/uniform.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random.hpp>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include<boost/math/distributions/normal.hpp>
#include<thread>

#endif
using namespace std;
typedef long unsigned int luint;

vector<int> txtReader(string);


//model parameters
struct modParms {
	double dE = 0.150; double dL = 0.269; double dP = 1.563; double uoE = 0.034; double uoL = 0.035; double uP = 0.25; double uM = 0.096; double Y = 13.25; int S = 3;
	int tr = 7; double sf1 = 20; double sf2 = 20; double sf3 = 20; double sf4 = 20; double sf5 = 20; double sf6 = 20; double dt = 0.25; double n = 50; double Emax = 93.6;
	int E0 = 177; int L0 = 8; int P0 = 1; int M0 = 7; double z1 = 5000; double z2 = 5000; double z3 = 5000; double z4 = 5000; double z5 = 5000; double z6 = 5000; double B = 21.19;
	int startTime = 0; int endTime = 2000; vector<int> rF; double w = 0.01; double sf = 20; double z = 20; int fxdPrm;
};

//obs dat data struct
struct obsDatX {
	vector <tuple<int, int>> garki408;
	vector <tuple<int, int>> garki202;
	vector <tuple<int, int>> garki154;
	vector <tuple<int, int>> garki155;
	vector <tuple<int, int>> garki218;
	vector <tuple<int, int>> garki304;
};


//vector<int> txtReader(string file);
int binom(int n, double p, boost::mt19937 rd);
luint rpois(luint lambda, boost::mt19937 rd);
double lbeta(double a, double b);
double betaBinom(double k, double n, double p, double w);
vector<tuple<int, int, int, int>> mPmod(modParms, boost::mt19937 rd);
tuple<int, int, int, int, double> modStepFnc(modParms wp, int obsData, boost::mt19937 rd);
double llFunc(int particles, modParms prms, obsDatX obsDat, int fixedParam);

//pfilter
double pFilt(int n,
	vector< tuple<int, int> > obsData,
	modParms prms,
	bool resM,
	int fxdParams);

//istate
vector<tuple<int, int, int, int, double>> iState(int N, int time, modParms iParms, int fxdParm);

//rand 0-1 generator
double rn01(void);

//data structure for MMH output
struct pMMHres {
	vector<double> uoE;
	vector<double> uoL;
	vector<double> uP;
	vector<double> Y;
	vector<double> w;
	vector<double> n;
	vector<double> z1;
	vector<double> z2;
	vector<double> z3;
	vector<double> z4;
	vector<double> z5;
	vector<double> z6;

	vector<double> sf1;
	vector<double> sf2;
	vector<double> sf3;
	vector<double> sf4;
	vector<double> sf5;
	vector<double> sf6;

	vector<double> dE;
	vector<double> dL;
	vector<double> dP;

	vector<double> ll;
};

//pMMH sampler 
pMMHres pMMHSampler(
	modParms initParams,
	int fixedParam,
	vector<double> sdProps,
	vector<double> acptRs,
	vector<tuple<string, double>> fitPrms,
	vector<double> maxSdProps,
	int niter,
	int particles,
	int nburn,
	bool monitoring,
	int startAdapt,
	int	tell,
	obsDatX	oDat);