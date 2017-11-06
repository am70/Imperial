#include<iostream>
#include"lModH.h"

#include <iostream>
#include <map>
#include <random>
#include <chrono>
#include <fstream>

int main()
{
	//std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	//std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	//std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << std::endl;

	double inf = std::numeric_limits<double>::infinity();
	modParms initParams;
	// read in garki village data
	obsDatX garkiDat;
	garkiDat.garki101 = { { 353,1 },{ 368,18 },{ 382,3 },{ 396,3 },{ 410,7 },{ 424,24 },{ 440,40 },{ 454,1 },{ 468,2 },{ 482,0 },{ 496,0 },{ 510,0 },{ 524,0 },{ 538,0 },{ 552,0 } };
	garkiDat.garki104 = { { 353,0 },{ 368,20 },{ 382,0 },{ 396,11 },{ 410,4 },{ 424,3 },{ 440,6 },{ 454,3 },{ 468,1 },{ 482,1 },{ 496,0 },{ 510,0 },{ 524,0 },{ 538,0 },{ 552,0 } };
	garkiDat.garki219 = { { 346,5 },{ 361,8 },{ 375,18 },{ 389,8 },{ 403,7 },{ 417,9 },{ 433,35 },{ 447,18 },{ 461,2 },{ 475,0 },{ 489,0 },{ 503,0 },{ 517,0 },{ 531,0 },{ 545,0 } };
	garkiDat.garki220 = { { 346,1 },{ 361,3 },{ 375,26 },{ 389,27 },{ 403,6 },{ 417,14 },{ 433,20 },{ 447,17 },{ 461,3 },{ 475,0 },{ 489,0 },{ 503,0 },{ 517,0 },{ 531,0 },{ 545,0 } };

	//modParms initParams;
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
	initParams.rF = txtReader("C:\\Imperial\\rf2.txt");

	vector<double> sdProps = { 0.01, 0.01, 0.1, 1, 0.01, 1, 1, 1, 0.1, 3, 3, 3, 3, 3 };
	vector<double> maxSdProps = { 0.1, 0.1, 0.1, inf, inf, inf, 1500000, 1500000, 1500000, 1500000, 1500000, 1500000, 1500000, 1500000 };
	vector<double> acptRs = { 0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.5,0.5,0.6,0.45,0.6,0.6,0.4 };
	vector<tuple<string, double>> fitPrms = { {"uoE", 0.03},{ "uoL", 0.031 },{ "uP", 0.249 },{ "Y", 11.5 },{ "w", 0.01 },{ "n", 50 },{ "z1", 1001 },{ "z2", 1002 },{ "z3", 1003 },{ "z4", 1004 },
	{"sf1", 57110} ,{ "sf2", 12660 } ,{ "sf3", 8097 } ,{ "sf4", 64342 }  };


	pMMHres results = pMMHSampler(
		initParams,//initial parameters
		4,//fixed parameter for days of rainfall
		sdProps,//initial sf for param proposals
		acptRs,//acceptance ratios
		fitPrms,//tuple of initial parm values plus names - needed as no reflection - maybe can be coded better
		maxSdProps,//max sd for each parameter proposal in tuner
		500000,//iterations
		75,//particles
		1,//nburn
		1,//monitoring
		150,//start adapt
		100,//tell
		garkiDat//observed data
	);
//write results to csv
	ofstream myfile;
	myfile.open("C:\\Imperial\\examp.txt");
		for (auto iter = 0; iter != size(results.ll); ++iter) {
			myfile << results.uoE.at(iter)<<" " << results.uoL.at(iter) << " " << results.uP.at(iter) << " " << results.Y.at(iter) << " " << results.w.at(iter) 
			<< " " << results.n.at(iter) << " " << results.z1.at(iter) << " " << results.z2.at(iter) << " " << results.z3.at(iter) << " " << results.z4.at(iter) 
			<< " " << results.sf1.at(iter) << " " << results.sf2.at(iter) << " " << results.sf3.at(iter) << " " << results.sf4.at(iter) << " " << results.ll.at(iter) << std::endl;
		}

		cout << "fin";

	cin.get();


	}





		
	
