#include<iostream>
#include"lModH.h"

#include <iostream>
#include <map>
#include <random>
#include <chrono>

int main()
{


	double inf = std::numeric_limits<double>::infinity();



	// read in garki village data
	obsDatX garkiDat;
	garkiDat.garki101 = { { 353,1 },{ 368,18 },{ 382,3 },{ 396,3 },{ 410,7 },{ 424,24 },{ 440,40 },{ 454,1 },{ 468,2 },{ 482,0 },{ 496,0 },{ 510,0 },{ 524,0 },{ 538,0 },{ 552,0 } };
	garkiDat.garki104 = { { 353,0 },{ 368,20 },{ 382,0 },{ 396,11 },{ 410,4 },{ 424,3 },{ 440,6 },{ 454,3 },{ 468,1 },{ 482,1 },{ 496,0 },{ 510,0 },{ 524,0 },{ 538,0 },{ 552,0 } };
	garkiDat.garki219 = { { 346,5 },{ 361,8 },{ 375,18 },{ 389,8 },{ 403,7 },{ 417,9 },{ 433,35 },{ 447,18 },{ 461,2 },{ 475,0 },{ 489,0 },{ 503,0 },{ 517,0 },{ 531,0 },{ 545,0 } };
	garkiDat.garki220 = { { 346,1 },{ 361,3 },{ 375,26 },{ 389,27 },{ 403,6 },{ 417,14 },{ 433,20 },{ 447,17 },{ 461,3 },{ 475,0 },{ 489,0 },{ 503,0 },{ 517,0 },{ 531,0 },{ 545,0 } };


	modParms initParams;

	initParams.uoE = 0.03; 
	initParams.uoL = 0.03156214;
	initParams.uP = 0.2499843;
	initParams.Y = 11.5012;
	initParams.z1 = 1001;
	initParams.z2 = 1002;
	initParams.z3 = 1003;
	initParams.z4 = 1004;
	initParams.w = 0.01;
	initParams.sf1 = 57110;
	initParams.sf2 = 12660;
	initParams.sf3 = 8097;
	initParams.sf4 = 64342;
	initParams.n = 50;

	initParams.startTime = 1384;
	initParams.endTime = 2180;
	initParams.sf = 57110;

	initParams.rF = txtReader("C:\\Imperial\\rf2.txt");

	initParams.z = 1003;
	initParams.sf = 8097;
	vector<tuple<int, int, int, int>> modRu = mPmod(initParams);
//	cin.get();

	vector<double> sdProps = { 0.01, 0.01, 0.1, 1, 1, 1, 1, 1, 0.1, 3, 3, 3, 3, 3 };
	vector<double> maxSdProps = { 0.1, 0.1, 0.1, inf, inf, inf, inf, inf, 0.1, 1500000, 1500000, 1500000, 1500000, inf };
	vector<double> acptRs = { 0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.5,0.5,0.6,0.45,0.6,0.6,0.4 };
	vector<tuple<string, double>> fitPrms = { {"uoE", 0.03},{ "uoL", 0.031 },{ "uP", 0.249 },{ "Y", 11.5 },{ "w", 0.01 },{ "z1", 1001 },{ "z2", 1002 },{ "z3", 1003 },{ "z4", 1004 },
	{"sf1", 57110} ,{ "sf2", 12660 } ,{ "sf3", 8097 } ,{ "sf4", 64342 } ,{ "n", 50 } };
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	pMMHres results = pMMHSampler(
		initParams,
		7,
		sdProps,
		acptRs,
		fitPrms,
		maxSdProps,
		200,//iterations
		50,//particles
		20,//nburn
		1,//monitoring
		150,//start adapt
		20,//tell
		garkiDat);
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << std::endl;

	cout << "fin";
	cin.get();
	//pMMH sampler 
	//pMMHres pMMHSampler(
	//	modParms initParams,
	//	int fixedParam,
	//	vector<double> sdProps,
	//	vector<double> acptRs,
	//	vector<tuple<string, double>> fitPrms,
	//	vector<double> maxSdProps,
	//	int niter,
	//	int particles,
	//	int nburn,
	//	int monitoring,
	//	int startAdapt,
	//	int	tell,
	//	obsDatX	oDat);
	//

	}





		
	
