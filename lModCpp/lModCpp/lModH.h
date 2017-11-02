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
#endif
using namespace std;
typedef long unsigned int luint;


vector<int> txtReader(string);


//model parameters
struct modParms {
	double dE = 0.150; double dL=0.269; double dP=1.563; double uoE=0.034; double uoL=0.035; double uP=0.25; double uM=0.096; double Y=13.25; int S=3;
	int tr = 7; double sf1 = 20; double sf2 = 20; double sf3 = 20; double sf4 = 20; double dt = 0.25; double n = 50; double Emax = 93.6; 
	int E0 = 177; int L0 = 8; int P0 = 1; int M0 = 7; double z1 = 5000; double z2 = 5000; double z3 = 5000; double z4 = 5000; double B = 21.19;
	int startTime = 0; int endTime = 2000; vector<int> rF; double w = 0.01; double sf =20; double z =20;
};

//obs dat data struct
struct obsDatX {
	vector <tuple<int, int>> garki101;
	vector <tuple<int, int>> garki104;
	vector <tuple<int, int>> garki219;
	vector <tuple<int, int>> garki220;
};

//data structure for particles to be weighted
struct wpStruct {
	int E0;
	int L0;
	int P0;
	int M0;
	int startTime;
	int endTime;
	double uoE;
	double uoL;
	double uP;
	double Y;
	double sf;
	vector<int> rFclust;
	double n;
	double fxdPrm;
	double w;
};

//vector<int> txtReader(string file);
int binom(int n, double p);
luint rpois(luint lambda);
double lbeta(double a, double b);
double betaBinom(double k, double n, double p, double w);
vector<tuple<int, int, int, int>> mPmod(modParms);
tuple<int, int, int, int, double> modStepFnc(wpStruct wp, int obsData);

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
	vector<double> sf1;
	vector<double> sf2;
	vector<double> sf3;
	vector<double> sf4;
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
	int monitoring,
	int startAdapt,
	int	tell,
	obsDatX	oDat);