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



	modParms prms;


	prms.uoE = 0.03397;
	prms.uoL = 0.03808;
	prms.uP = 0.1908;
	prms.Y = 17.7535;
	prms.w = 0.0100;
	prms.n = 2.275;
	prms.z1 = pow(10, 1.986);
	prms.z2 = pow(10, 2.61);
	prms.z3 = pow(10, 2.85);
	prms.z4 = pow(10, 2.28);
	prms.z = pow(10, 1.94375);
	prms.sf1 = pow(10, 4.714);
	prms.sf2 = pow(10, 3.94);
	prms.sf3 = pow(10, 4.06);
	prms.sf4 = pow(10, 4.111);
	prms.sf = 10837.5;
	prms.rF = txtReader("Q:\\Imperial\\rf1.txt");
	prms.E0 = 87;
	prms.L0 = 28;
	prms.P0 = 13;
	prms.M0 = 111;
	prms.startTime = 1408;
	prms.endTime = 2204;
	prms.tr = 4;
	prms.dt = 0.25;

	boost::mt19937 mrandThread(std::random_device{}());


	cout << "start" << endl;
	vector<tuple<int, int, int, int>> test2 = mPmod(prms, mrandThread);


	ofstream myfile1;
	myfile1.open("Q:\\Imperial\\test.txt");
	for (auto iter = 0; iter != prms.endTime - prms.startTime; ++iter) {
		myfile1 << get<3>(test2[iter]) << endl;
	}


	double test = pFilt(1,
		garkiDat.garki101,
		prms,
		false,
		7);
	cout << test;
	//cin.get();



	//modParms initParams;
	initParams.uoE = 0.035;
	initParams.uoL = 0.03156214;
	initParams.uP = 0.2499843;
	initParams.Y = 11.5012;
	initParams.z1 = 3;
	initParams.z2 = 3;
	initParams.z3 = 3;
	initParams.z4 = 3;

	initParams.w = 0.08;
	initParams.sf1 = 20000;
	initParams.sf2 = 20000;
	initParams.sf3 = 20000;
	initParams.sf4 = 20000;

	initParams.dE = 0.15;
	initParams.dL = 0.2;
	initParams.dP = 1;


	initParams.n = 20;
	initParams.rF = txtReader("Q:\\Imperial\\rf2.txt");

	vector<double> sdProps = { 0.01, 0.01, 0.1, 1, 0.001, 1, 1, 1,1,1, 10, 10, 10, 10,0.05,0.05,0.1 };
	vector<double> maxSdProps = { 0.1, 0.1, 0.1, 20, 0.1, 5, 1, 1, 1, 1, 1e+5, 1e+5, 1e+5, 1e+5,0.1,0.1,2 };
	vector<double> acptRs = { 0.3,0.3,0.3,0.3,0.3,0.5,0.5,0.5,0.5,0.5,0.6,0.6,0.6,0.6,0.5,0.5,0.5, };
	vector<tuple<string, double>> fitPrms = { { "uoE", initParams.uoE },{ "uoL", initParams.uoL },{ "uP", initParams.uP },{ "Y", initParams.Y },
	{ "w", initParams.w },{ "n", initParams.n },{ "z1", initParams.z1 },{ "z2", initParams.z2 },{ "z3", initParams.z3 },{ "z4", initParams.z4 },
	{ "sf1", initParams.sf1 } ,{ "sf2", initParams.sf2 } ,{ "sf3", initParams.sf3 } ,{ "sf4", initParams.sf4 } ,{ "dE", initParams.dE },{ "dL", initParams.dL } ,{ "dP", initParams.dP } };



	pMMHres results = pMMHSampler(
		initParams,//initial parameters
		7,//fixed parameter for days of rainfall   
		sdProps,//initial sf for param proposals
		acptRs,//acceptance ratios
		fitPrms,//tuple of initial parm values plus names - needed as no reflection - maybe can be coded better
		maxSdProps,//max sd for each parameter proposal in tuner
		1000000,//iterations
		75,//particles
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
	myfile.open("Q:\\Imperial\\1milTestTr7.txt");
	for (auto iter = 0; iter != size(results.ll); ++iter) {
		myfile << results.uoE.at(iter) << " " << results.uoL.at(iter) << " " << results.uP.at(iter) << " " << results.Y.at(iter) << " " << results.w.at(iter)
			<< " " << results.n.at(iter) << " " << results.z1.at(iter) << " " << results.z2.at(iter) << " " << results.z3.at(iter) << " " << results.z4.at(iter) 
			<< " " << results.sf1.at(iter) << " " << results.sf2.at(iter) << " " << results.sf3.at(iter) << " " << results.sf4.at(iter) << " " << results.dE.at(iter) 
			<< " " << results.dL.at(iter) << " " << results.dP.at(iter) << " " << results.ll.at(iter) << std::endl;
	}
	cout << "end";
	cin.get();

}







