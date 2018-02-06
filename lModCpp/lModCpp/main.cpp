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

	modParms initParams;
	obsDatX garkiDatGam;

// read in garki village data anoph gamb 
	garkiDatGam.garki154 = { { 15,3 },{ 29,4 },{ 43,4 },{ 57,6 },{ 71,17 },{ 85,7 },{ 99,20 },{ 113,34 },{ 141,18 },{ 155,5 },{ 169,4 },{ 183,0 },{ 196,0 } , { 210,0 }, { 224,0 }, { 238,1 }, { 252,0 }, { 266,0 }};
	garkiDatGam.garki202 = { { 8,18 },{ 22,15 },{ 36,20 },{ 50,167 },{ 64,155 },{ 78,370 },{ 92,178 },{ 106,127 },{ 120,75 },{ 134,25 },{ 147,5 },{ 162,3 },{ 175,0 } , { 189,0 }, { 203,1 }, { 216,0 }, { 231,0 }, { 245,0 } };
	garkiDatGam.garki218 = { { 19,99 },{ 33,266 },{ 47,211 },{ 61,301 },{ 75,475 },{ 89,309 },{ 103,279 },{ 117,52 },{ 131,14 },{ 147,6 },{ 156,10 }, { 172,0 }, { 186,0 }, { 200,0 }, { 215,0 }, { 228,0 }, { 242,0 }};
	garkiDatGam.garki304 = { { 18,3 },{ 32,20 },{ 46,19 },{ 60,80 },{ 74,69 },{ 88,357 },{ 102,310 },{ 116,151 },{ 130,141 },{ 144,129 },{ 158,15 },{ 172,2 },{ 186,1 },{ 197,0 },{ 211,0 } , { 225,0 }, { 239,0 }, { 253,0 }, { 267,0 } };
	garkiDatGam.garki553 = { { 12,112 },{ 26,143 },{ 40,656 },{ 54,250 },{ 68,673 },{ 82,683 },{ 96,336 },{ 110,124 },{ 124,19 },{ 138,2 },{ 153,7 },{ 165,13 },{ 179,1 },{ 193,0 } , { 207,0 }, { 221,0 } };
	garkiDatGam.garki802 = { { 18,69 },{ 32,285 },{ 46,66 },{ 60,20 },{ 74,103 },{ 88,67 },{ 102,6 },{ 116,8 },{ 130,3 },{ 144,2 },{ 158,0 },{ 171,2 } , { 185,0 }};


	//modParms initParams;
	initParams.uoE = 0.0598336;
	initParams.uoL = 0.0398348;
	initParams.uP = 0.18053;
	initParams.uM = 0.0613272;
	initParams.Y = 13.44607;
	initParams.w = 0.0001;
	initParams.n = 11.6496;

	initParams.z1 = 7.58718;
	initParams.z2 = 4.91746;
	initParams.z3 = 5.87573;
	initParams.z4 = 4.90541;
	initParams.z5 = 7.81856;
	initParams.z6 = 6.13859;

	initParams.sf1 = 10.6058;
	initParams.sf2 = 7.30038;
	initParams.sf3 = 7.53952;
	initParams.sf4 = 13.4258;
	initParams.sf5 = 7.34137;
	initParams.sf6 = 6.61745;

	initParams.dE = 0.00751125;
	initParams.dL = 0.420767;
	initParams.dP = 1.43223;
	initParams.o = 6.11526;
	initParams.dt = 0.25;
	initParams.Mg = 184.101;

	initParams.p = 0.0001;

	initParams.rF = txtReader("Q:\\Imperial\\lModCpp\\Data\\rf01.txt", 0.25);





	vector<double> sdProps = { 
		0.001, 0.001, 0.01,0.01,5,0.001, 
		1, 1, 1,1,1,1,
		1, 1, 1, 1, 1,1,
		1,0.05,0.05,0.05,1,40,0.001};
	vector<double> maxSdProps = {
		0.05, 0.05, 0.8, 0.5,5,0.01,
		3, 3, 3, 3,3,3,
		4, 4, 4, 4,4,4,
		0.1,0.1,0.1,0.1,1,250,0.001};
	vector<double> acptRs = {
		0.25,0.25,0.25,0.25,0.25,0.25,
		0.25,0.25,0.25,0.25,0.25,0.25,
		0.25,0.25,0.25,0.25,0.25,0.25,
		0.25,0.25,0.25,0.25,0.25,0.25,0.25 };

	vector<tuple<string, double>> fitPrms = { { "uoE", initParams.uoE },{ "uoL", initParams.uoL },{ "uP", initParams.uP },{ "uM", initParams.uM },{ "Y", initParams.Y },
	{ "w", initParams.w },{ "n", initParams.n },{ "z1", initParams.z1 },{ "z2", initParams.z2 },{ "z3", initParams.z3 },{ "z4", initParams.z4 },{ "z5", initParams.z5 },{ "z6", initParams.z6 },
	{ "sf1", initParams.sf1 } ,{ "sf2", initParams.sf2 } ,{ "sf3", initParams.sf3 } ,{ "sf4", initParams.sf4 } ,{ "sf5", initParams.sf5 } ,{ "sf6", initParams.sf6 }
	,{ "dE", initParams.dE },{ "dL", initParams.dL } ,{ "dP", initParams.dP } ,{ "o", initParams.o },{ "Mg", initParams.Mg },{"p",initParams.p } };



	pMMHres results = pMMHSampler(
		initParams,//initial parameters
		6,//fixed parameter(s) if any - obsolete? 
		sdProps,//initial sf for param proposals
		acptRs,//acceptance ratios
		fitPrms,//tuple of initial parm values plus names - needed as no reflection - maybe can be coded better
		maxSdProps,//max sd for each parameter proposal in tuner
		50000,//iterations
		45,//particles
		10000,//nburn 
		1,//monitoring
		1000,//start adapt
		25,//tell
		garkiDatGam//observed data
	);

	//std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	//std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << std::endl;

	//write results to csv
	ofstream myfile;
	myfile.open("Q:\\Imperial\\1milTestTr7_TEST3.txt");
	for (auto iter = 0; iter != size(results.ll); ++iter) {
		myfile << results.uoE.at(iter) << " " << results.uoL.at(iter) << " " << results.uP.at(iter) << " " << results.Y.at(iter) << " " << results.w.at(iter)
			<< " " << results.n.at(iter) << " " << results.z1.at(iter) << " " << results.z2.at(iter) << " " << results.z3.at(iter) << " " << results.z4.at(iter) << " " << results.z5.at(iter) << " " << results.z6.at(iter)
			<< " " << results.sf1.at(iter) << " " << results.sf2.at(iter) << " " << results.sf3.at(iter) << " " << results.sf4.at(iter) << " " << results.sf5.at(iter) << " " << results.sf6.at(iter) << 
			" " << results.dE.at(iter) << " " << results.dL.at(iter) << " " << results.dP.at(iter) << " "  << results.o.at(iter)<< " " << results.uM.at(iter) << " " << results.Mg.at(iter) << " " << results.p.at(iter)<< " " << results.ll.at(iter) <<endl;
	}
	

	pFitFunc(75, results, garkiDatGam, 20, initParams);
	cout << "End";
	cin.get();

	//std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

}







