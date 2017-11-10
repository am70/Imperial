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
	initParams.uoE = 0.035; 
	initParams.uoL = 0.03156214;
	initParams.uP = 0.2499843;
	initParams.Y = 11.5012;
	initParams.z1 = 5000;
	initParams.z2 = 5000;
	initParams.z3 = 5000;
	initParams.z4 = 5000;
	initParams.w = 0.01;
	initParams.sf1 = 57110;
	initParams.sf2 = 12660;
	initParams.sf3 = 8097;
	initParams.sf4 = 64342;
	initParams.n = 90;
	initParams.rF = txtReader("Q:\\Imperial\\rf2.txt");

	vector<double> sdProps = { 0.01, 0.01, 0.1, 1, 0.01, 0, 1, 1, 0.1, 3, 3, 3, 3, 3 };
	vector<double> maxSdProps = { 0.1, 0.1, 0.1, 20, 0.1, 0, 5e+5, 5e+5, 5e+5, 5e+5, 5e+5, 5e+5, 5e+5, 5e+5, };
	vector<double> acptRs = { 0.3,0.3,0.3,0.3,0.3,0.5,0.5,0.5,0.5,0.6,0.6,0.6,0.6,0.6 };
	vector<tuple<string, double>> fitPrms = { {"uoE", initParams.uoE },{ "uoL", initParams.uoL },{ "uP", initParams.uP },{ "Y", initParams.Y },
	{ "w", initParams.w },{ "n", initParams.n },{ "z1", initParams.z1 },{ "z2", initParams.z2 },{ "z3", initParams.z3 },{ "z4", initParams.z4 },
	{"sf1", initParams.sf1 } ,{ "sf2", initParams.sf2 } ,{ "sf3", initParams.sf3 } ,{ "sf4", initParams.sf4 }};


	pMMHres results = pMMHSampler(
		initParams,//initial parameters
		7,//fixed parameter for days of rainfall
		sdProps,//initial sf for param proposals
		acptRs,//acceptance ratios
		fitPrms,//tuple of initial parm values plus names - needed as no reflection - maybe can be coded better
		maxSdProps,//max sd for each parameter proposal in tuner
		100000,//iterations
		75,//particles
		10000,//nburn
		1,//monitoring
		500,//start adapt
		100,//tell
		garkiDat//observed data
	);
//write results to csv
	ofstream myfile;
	myfile.open("Q:\\Imperial\\100kn90.txt");
		for (auto iter = 0; iter != size(results.ll); ++iter) {
			myfile << results.uoE.at(iter)<<" " << results.uoL.at(iter) << " " << results.uP.at(iter) << " " << results.Y.at(iter) << " " << results.w.at(iter) 
			<< " " << results.n.at(iter) << " " << results.z1.at(iter) << " " << results.z2.at(iter) << " " << results.z3.at(iter) << " " << results.z4.at(iter) 
			<< " " << results.sf1.at(iter) << " " << results.sf2.at(iter) << " " << results.sf3.at(iter) << " " << results.sf4.at(iter) << " " << results.ll.at(iter) << std::endl;
		}



		initParams.n = 10;
		vector<tuple<string, double>> fitPrms2 = { { "uoE", initParams.uoE },{ "uoL", initParams.uoL },{ "uP", initParams.uP },{ "Y", initParams.Y },
		{ "w", initParams.w },{ "n", initParams.n },{ "z1", initParams.z1 },{ "z2", initParams.z2 },{ "z3", initParams.z3 },{ "z4", initParams.z4 },
		{ "sf1", initParams.sf1 } ,{ "sf2", initParams.sf2 } ,{ "sf3", initParams.sf3 } ,{ "sf4", initParams.sf4 } };

		vector<double> sdProps2 = { 0.01, 0.01, 0.1, 1, 0.01, 0, 1, 1, 0.1, 3, 3, 3, 3, 3 };
		vector<double> maxSdProps2 = { 0.1, 0.1, 0.1, 20, 0.1, 0, 5e+5, 5e+5, 5e+4, 5e+5, 5e+5, 5e+5, 5e+5, 5e+5, };
		pMMHres results2 = pMMHSampler(
			initParams,//initial parameters
			7,//fixed parameter for days of rainfall
			sdProps2,//initial sf for param proposals
			acptRs,//acceptance ratios
			fitPrms2,//tuple of initial parm values plus names - needed as no reflection - maybe can be coded better
			maxSdProps2,//max sd for each parameter proposal in tuner
			100000,//iterations
			75,//particles
			10000,//nburn
			1,//monitoring
			500,//start adapt
			100,//tell
			garkiDat//observed data
		);
		//write results to csv
		ofstream myfile2;
		myfile2.open("Q:\\Imperial\\100kn10.txt");
		for (auto iter = 0; iter != size(results2.ll); ++iter) {
			myfile2 << results2.uoE.at(iter) << " " << results2.uoL.at(iter) << " " << results2.uP.at(iter) << " " << results2.Y.at(iter) << " " << results2.w.at(iter)
				<< " " << results2.n.at(iter) << " " << results2.z1.at(iter) << " " << results2.z2.at(iter) << " " << results2.z3.at(iter) << " " << results2.z4.at(iter)
				<< " " << results2.sf1.at(iter) << " " << results2.sf2.at(iter) << " " << results2.sf3.at(iter) << " " << results2.sf4.at(iter) << " " << results2.ll.at(iter) << std::endl;
		}


		initParams.n = 30;

		vector<tuple<string, double>> fitPrms3 = { { "uoE", initParams.uoE },{ "uoL", initParams.uoL },{ "uP", initParams.uP },{ "Y", initParams.Y },
		{ "w", initParams.w },{ "n", initParams.n },{ "z1", initParams.z1 },{ "z2", initParams.z2 },{ "z3", initParams.z3 },{ "z4", initParams.z4 },
		{ "sf1", initParams.sf1 } ,{ "sf2", initParams.sf2 } ,{ "sf3", initParams.sf3 } ,{ "sf4", initParams.sf4 } };

		vector<double> sdProps3 = { 0.01, 0.01, 0.1, 1, 0.01, 0, 1, 1, 0.1, 3, 3, 3, 3, 3 };
		vector<double> maxSdProps3 = { 0.1, 0.1, 0.1, 20, 0.1, 0, 5e+5, 5e+5, 5e+5, 5e+5, 5e+5, 5e+5, 5e+5, 5e+5, };
		pMMHres results3 = pMMHSampler(
			initParams,//initial parameters
			7,//fixed parameter for days of rainfall
			sdProps3,//initial sf for param proposals
			acptRs,//acceptance ratios
			fitPrms3,//tuple of initial parm values plus names - needed as no reflection - maybe can be coded better
			maxSdProps3,//max sd for each parameter proposal in tuner
			100000,//iterations
			75,//particles
			10000,//nburn
			1,//monitoring
			500,//start adapt
			100,//tell
			garkiDat//observed data
		);
		//write results to csv
		ofstream myfile3;
		myfile3.open("Q:\\Imperial\\100kn30.txt");
		for (auto iter = 0; iter != size(results3.ll); ++iter) {
			myfile3 << results3.uoE.at(iter) << " " << results3.uoL.at(iter) << " " << results3.uP.at(iter) << " " << results3.Y.at(iter) << " " << results3.w.at(iter)
				<< " " << results3.n.at(iter) << " " << results3.z1.at(iter) << " " << results3.z2.at(iter) << " " << results3.z3.at(iter) << " " << results3.z4.at(iter)
				<< " " << results3.sf1.at(iter) << " " << results3.sf2.at(iter) << " " << results3.sf3.at(iter) << " " << results3.sf4.at(iter) << " " << results3.ll.at(iter) << std::endl;
		}



		vector<tuple<string, double>> fitPrms4 = { { "uoE", initParams.uoE },{ "uoL", initParams.uoL },{ "uP", initParams.uP },{ "Y", initParams.Y },
		{ "w", initParams.w },{ "n", initParams.n },{ "z1", initParams.z1 },{ "z2", initParams.z2 },{ "z3", initParams.z3 },{ "z4", initParams.z4 },
		{ "sf1", initParams.sf1 } ,{ "sf2", initParams.sf2 } ,{ "sf3", initParams.sf3 } ,{ "sf4", initParams.sf4 } };
		initParams.n = 50;
		vector<double> sdProps4 = { 0.01, 0.01, 0.1, 1, 0.01, 0, 1, 1, 0.1, 3, 3, 3, 3, 3 };
		vector<double> maxSdProps4 = { 0.1, 0.1, 0.1, 20, 0.1, 0, 5e+5, 5e+5, 5e+5, 5e+5, 5e+5, 5e+5, 5e+5, 5e+5, };
		pMMHres results4 = pMMHSampler(
			initParams,//initial parameters
			7,//fixed parameter for days of rainfall
			sdProps4,//initial sf for param proposals
			acptRs,//acceptance ratios
			fitPrms4,//tuple of initial parm values plus names - needed as no reflection - maybe can be coded better
			maxSdProps4,//max sd for each parameter proposal in tuner
			100000,//iterations
			75,//particles
			10000,//nburn
			1,//monitoring
			500,//start adapt
			100,//tell
			garkiDat//observed data
		);
		//write results to csv
		ofstream myfile4;
		myfile4.open("Q:\\Imperial\\100kn50.txt");
		for (auto iter = 0; iter != size(results4.ll); ++iter) {
			myfile4 << results4.uoE.at(iter) << " " << results4.uoL.at(iter) << " " << results4.uP.at(iter) << " " << results4.Y.at(iter) << " " << results4.w.at(iter)
				<< " " << results4.n.at(iter) << " " << results4.z1.at(iter) << " " << results4.z2.at(iter) << " " << results4.z3.at(iter) << " " << results4.z4.at(iter)
				<< " " << results4.sf1.at(iter) << " " << results4.sf2.at(iter) << " " << results4.sf3.at(iter) << " " << results4.sf4.at(iter) << " " << results4.ll.at(iter) << std::endl;
		}


		initParams.n = 70;
		vector<tuple<string, double>> fitPrms5 = { { "uoE", initParams.uoE },{ "uoL", initParams.uoL },{ "uP", initParams.uP },{ "Y", initParams.Y },
		{ "w", initParams.w },{ "n", initParams.n },{ "z1", initParams.z1 },{ "z2", initParams.z2 },{ "z3", initParams.z3 },{ "z4", initParams.z4 },
		{ "sf1", initParams.sf1 } ,{ "sf2", initParams.sf2 } ,{ "sf3", initParams.sf3 } ,{ "sf4", initParams.sf4 } };
		vector<double> sdProps5 = { 0.01, 0.01, 0.1, 1, 0.01, 0, 1, 1, 0.1, 3, 3, 3, 3, 3 };
		vector<double> maxSdProps5 = { 0.1, 0.1, 0.1, 20, 0.1, 0, 5e+5, 5e+5, 5e+5, 5e+5, 5e+5, 5e+5, 5e+5, 5e+5, };
		pMMHres results5 = pMMHSampler(
			initParams,//initial parameters
			7,//fixed parameter for days of rainfall
			sdProps5,//initial sf for param proposals
			acptRs,//acceptance ratios
			fitPrms5,//tuple of initial parm values plus names - needed as no reflection - maybe can be coded better
			maxSdProps5,//max sd for each parameter proposal in tuner
			100000,//iterations
			75,//particles
			10000,//nburn
			1,//monitoring
			500,//start adapt
			100,//tell
			garkiDat//observed data
		);
		//write results to csv
		ofstream myfile5;
		myfile5.open("Q:\\Imperial\\100kn70.txt");
		for (auto iter = 0; iter != size(results5.ll); ++iter) {
			myfile5 << results5.uoE.at(iter) << " " << results5.uoL.at(iter) << " " << results5.uP.at(iter) << " " << results5.Y.at(iter) << " " << results5.w.at(iter)
				<< " " << results5.n.at(iter) << " " << results5.z1.at(iter) << " " << results5.z2.at(iter) << " " << results5.z3.at(iter) << " " << results5.z4.at(iter)
				<< " " << results5.sf1.at(iter) << " " << results5.sf2.at(iter) << " " << results5.sf3.at(iter) << " " << results5.sf4.at(iter) << " " << results5.ll.at(iter) << std::endl;
		}





		cout << "fin";

	cin.get();


	}





		
	
