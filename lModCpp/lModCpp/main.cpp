#include<iostream>
#include"lModH.h"

#include <iostream>
#include <map>
#include <random>
#include <chrono>
#include <fstream>

int main()
{



	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();


	double inf = std::numeric_limits<double>::infinity();
	modParms initParams;
	// read in garki village data
	obsDatX garkiDat;
	garkiDat.garki101 = { { 353,1 },{ 368,18 },{ 382,3 },{ 396,3 },{ 410,7 },{ 424,24 },{ 440,40 },{ 454,1 },{ 468,2 },{ 482,0 },{ 496,0 },{ 510,0 },{ 524,0 },{ 538,0 },{ 552,0 } };
	garkiDat.garki104 = { { 353,0 },{ 368,20 },{ 382,0 },{ 396,11 },{ 410,4 },{ 424,3 },{ 440,6 },{ 454,3 },{ 468,1 },{ 482,1 },{ 496,0 },{ 510,0 },{ 524,0 },{ 538,0 },{ 552,0 } };
	garkiDat.garki219 = { { 346,5 },{ 361,8 },{ 375,18 },{ 389,8 },{ 403,7 },{ 417,9 },{ 433,35 },{ 447,18 },{ 461,2 },{ 475,0 },{ 489,0 },{ 503,0 },{ 517,0 },{ 531,0 },{ 545,0 } };
	garkiDat.garki220 = { { 346,1 },{ 361,3 },{ 375,26 },{ 389,27 },{ 403,6 },{ 417,14 },{ 433,20 },{ 447,17 },{ 461,3 },{ 475,0 },{ 489,0 },{ 503,0 },{ 517,0 },{ 531,0 },{ 545,0 } };
	garkiDat.garki782 = { { 1109,0 },{ 1123,0 },{ 1137,2 },{ 1152,6 },{ 1165,39 },{ 1179,14 },{ 1195,4 },{ 1206,14 },{ 1221,34 },{ 1235,16 },{ 1249,6 },{ 1263,6 },{ 1279,0 }};
	garkiDat.garki783 = { { 1109,0 },{ 1123,0 },{ 1137,3 },{ 1152,39 },{ 1165,34 },{ 1179,8 },{ 1195,2 },{ 1206,6 },{ 1221,3 },{ 1235,3 },{ 1249,0 },{ 1263,2 },{ 1279,0 } };




	//modParms initParams;
	initParams.uoE = 0.035; 
	initParams.uoL = 0.03156214;
	initParams.uP = 0.2499843;
	initParams.Y = 11.5012;
	initParams.z1 = 3;
	initParams.z2 = 3;
	initParams.z3 = 3;
	initParams.z4 = 3;
	initParams.z5 = 3;
	initParams.z6 = 3;
	initParams.w = 0.01;
	initParams.sf1 = 50000;
	initParams.sf2 = 40000;
	initParams.sf3 = 40000;
	initParams.sf4 = 50000;
	initParams.sf5 = 50000;
	initParams.sf6 = 50000;
	initParams.n = 50;
	initParams.rF = txtReader("Q:\\Imperial\\rf2.txt");

	vector<double> sdProps = { 0.01, 0.01, 0.1, 1, 0.001, 1, 1, 1,1,1,1,1, 10, 10, 10, 10, 10, 10 };
	vector<double> maxSdProps = { 0.1, 0.1, 0.1, 20, 0.1, 30, 1, 1, 1, 1,1,1, 1e+5, 1e+5, 1e+5, 1e+5,1e+5, 1e+5 };
	vector<double> acptRs = { 0.3,0.3,0.3,0.3,0.3,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6 };
	vector<tuple<string, double>> fitPrms = { {"uoE", initParams.uoE },{ "uoL", initParams.uoL },{ "uP", initParams.uP },{ "Y", initParams.Y },
	{ "w", initParams.w },{ "n", initParams.n },{ "z1", initParams.z1 },{ "z2", initParams.z2 },{ "z3", initParams.z3 },{ "z4", initParams.z4 },{ "z5", initParams.z4 },{ "z6", initParams.z4 },
	{"sf1", initParams.sf1 } ,{ "sf2", initParams.sf2 } ,{ "sf3", initParams.sf3 } ,{ "sf4", initParams.sf4 } ,{ "sf5", initParams.sf5 }  ,{ "sf6", initParams.sf6 } };
	






	pMMHres results = pMMHSampler(
		initParams,//initial parameters
		7,//fixed parameter for days of rainfall
		sdProps,//initial sf for param proposals
		acptRs,//acceptance ratios
		fitPrms,//tuple of initial parm values plus names - needed as no reflection - maybe can be coded better
		maxSdProps,//max sd for each parameter proposal in tuner
		50000,//iterations
		40,//particles
		50000,//nburn 
		1,//monitoring
		500,//start adapt
		100,//tell
		garkiDat//observed data
	);
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << std::endl;

//write results to csv
	ofstream myfile;
	myfile.open("Q:\\Imperial\\500k6vill.txt");
		for (auto iter = 0; iter != size(results.ll); ++iter) {
			myfile << results.uoE.at(iter)<<" " << results.uoL.at(iter) << " " << results.uP.at(iter) << " " << results.Y.at(iter) << " " << results.w.at(iter) 
			<< " " << results.n.at(iter) << " " << results.z1.at(iter) << " " << results.z2.at(iter) << " " << results.z3.at(iter) << " " << results.z4.at(iter) << " " << results.z5.at(iter) << " " << results.z6.at(iter)
			<< " " << results.sf1.at(iter) << " " << results.sf2.at(iter) << " " << results.sf3.at(iter) << " " << results.sf4.at(iter) << " " << results.sf5.at(iter)<< " " << results.sf6.at(iter)<< " " << results.ll.at(iter) << std::endl;
		}
		cout << "end";
	cin.get();

	}





		
	
