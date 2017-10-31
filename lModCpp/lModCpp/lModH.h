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

//vector<tuple<int, int>>garki101 = { { 353,1 },{ 368,18 },{ 382,3 },{ 396,3 },{ 410,7 },{ 424,24 },{ 440,40 },{ 454,1 },{ 468,2 },{ 482,0 },{ 496,0 },{ 510,0 },{ 524,0 },{ 538,0 },{ 552,0 } };

vector<int> txtReader(string);


//model parameters
struct modParms {
	double dE = 0.150; double dL=0.269; double dP=1.563; double uoE=0.034; double uoL=0.035; double uP=0.25; double uM=0.096; double Y=13.25; int S=3;
	int tr = 7; double sf1 = 1000; double sf2 = 1000; double sf3 = 1000; double sf4 = 1000; double dt = 0.25; double n = 50; double Emax = 93.6; 
	int E0 = 177; int L0 = 8; int P0 = 1; int M0 = 7; double z1 = 5000; double z2 = 5000; double z3 = 5000; double z4 = 5000; double B = 21.19;
	int startTime = 500; int endTime = 1000; vector<int> rF; double w = 0.01;
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
	double uoE;
	double uoL;
	double uP;
	double Y;
	double w;
	double n;
	double z1;
	double z2;
	double z3;
	double z4;
	double sf1;
	double sf2;
	double sf3;
	double sf4;
	double ll;
};