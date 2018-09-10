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
#include <windows.h>
#include <stdio.h>
#include <stdio.h>

#endif
using namespace std;
using std::cerr;

typedef long unsigned int luint;

vector<double> txtReader(string, double);
double medianFnc(vector<double> vec );

//model parameters
struct modParms {
	double dE = 0.150; double dL = 0.269; double dP = 1.563; double uoE = 0.034; double uoL = 0.035; double uP = 0.25; double uM = 0.096; double Y = 13.25; int S = 3;
	int tr = 7; double sf1 = 20; double sf2 = 20; double sf3 = 20; double sf4 = 20; double sf5 = 20; double sf6 = 20; double sf7 = 20; double sf8 = 20; double dt = 0.25; double n = 50; double Emax = 93.6;
	double E0 = 177; double L0 = 8; double P0 = 1; double M0 = 7; double z1 = 5000; double z2 = 5000; double z3 = 5000; double z4 = 5000; double z5 = 5000; double z6 = 5000; double z7 = 5000; double z8 = 5000; double B = 21.19;
	int startTime = 0; int endTime = 2000; vector<double> rF; double w = 0.01; double sf = 20; double z = 20; int fxdPrm; double o; double Mg; double p; double tau; double lK; double lKs; double lKm;
};

//obs dat data struct
struct obsDatX {
	vector <tuple<int, int>> garki154;
	vector <tuple<int, int>> garki218;
	vector <tuple<int, int>> garki304;
	vector <tuple<int, int>> garki202;
	vector <tuple<int, int>> garki553;
	vector <tuple<int, int>> garki802;
	vector <tuple<int, int>> garki408;
	vector <tuple<int, int>> garki801;
	vector <tuple<int, int>> garki801_2;
	vector <tuple<int, int>> garki802_2;
	vector <tuple<int, int>> garki553_2;
};


double nBgP(double k, double n, double p);

//vector<int> txtReader(string file);
vector<tuple<double, double, double, double, double>> normalise(vector<tuple<double, double, double, double, double>> particles, double llSum);
vector<tuple<double, double, double, double, double>> rSamp(vector<std::tuple<double, double, double, double, double>>& samp);
int binom(int n, double p, boost::mt19937 rd);
luint rpois(luint lambda, boost::mt19937 rd);
double lbeta(double a, double b);
double betaBinom(double k, double n, double p, double w);
vector<tuple<double, double, double, double,double>> mPmod(modParms, boost::mt19937_64 rd, string dFunc);
tuple<double, double, double, double, double> modStepFnc(modParms wp, int obsData, boost::mt19937_64 rd, string dFunc);
double llFunc(int particles, modParms prms, obsDatX obsDat, string dFunc);
double dbinom(double k, double n, double p);
double nB(double k, double n, double r, double p);

//pfilter
double pFilt(int n,
	vector<tuple<int, int>> obsData,
	modParms prms,
	bool resM,
	string outputFile,
	bool reff,
	string dFunc
);

//istate
vector<tuple<double, double, double, double, double>> iState(int N, int time, modParms iParms, string dFunc);

//rand 0-1 generator
double rn01(void);

//data structure for MMH output
struct pMMHres {
	vector<double> uoE;
	vector<double> uoL;
	vector<double> uP;
	vector<double> uM;
	vector<double> Y;
	vector<double> w;
	vector<double> n;
	vector<double> z1;
	vector<double> z2;
	vector<double> z3;
	vector<double> z4;
	vector<double> z5;
	vector<double> z6;
	vector<double> z7;
	vector<double> z8;

	vector<double> sf1;
	vector<double> sf2;
	vector<double> sf3;
	vector<double> sf4;
	vector<double> sf5;
	vector<double> sf6;
	vector<double> sf7;
	vector<double> sf8;

	vector<double> dE;
	vector<double> dL;
	vector<double> dP;
	vector<double> o;
	vector<double> Mg;
	vector<double> p;
	vector<double> tau;
	vector<double> lK;
	vector<double> lKs;
	vector<double> lKm;

	vector<double> ll;
};

//pMMH sampler 
pMMHres pMMHSampler(
	modParms initParams,
	string dFunc,
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

void pFitFunc(int particles, pMMHres results, obsDatX obsDat, modParms prms, string outputFile, string dFunc);

double resultsWriter(string fileName, string outputFolder, pMMHres results);

struct pmcmcOptions {
	string outputFolder;
	string dFunc;
	int iterations;
	int particles;
	int nburn;
	int monitoring;
	int startAdapt;
	int tell;
	string initParamsLoc;
};


pmcmcOptions optionsReader(string optionsTextFile, pmcmcOptions pmcmcOpt);
modParms initParamsReader(string initParamsTxtFile, modParms initParams);