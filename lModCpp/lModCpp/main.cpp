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

	garkiDatGam.garki408 = { { 5,7 },{ 19,99 },{ 33,670 },{ 47,1266 },{ 61,876 },{ 75,1340 },{ 89,227 },{ 103,513 },{ 117,49 },{ 131,19 },{ 145,15 },{ 158,7 } ,{ 172,15 },{186,2} };

	garkiDatGam.garki202 = { { 8,18 },{ 22,15 },{ 36,20 },{ 50,167 },{ 64,155 },{ 78,370 },{ 92,178 },{ 106,127 },{ 120,75 },{ 134,25 },{ 147,5 },{ 162,3 },{ 175,0 } , { 189,0 }, { 203,1 }, { 216,0 }, { 231,0 }, { 245,0 } };
	garkiDatGam.garki218 = { { 19,99 },{ 33,266 },{ 47,211 },{ 61,301 },{ 75,475 },{ 89,309 },{ 103,279 },{ 117,52 },{ 131,14 },{ 147,6 },{ 156,10 }, { 172,0 }, { 186,0 }, { 200,0 }, { 215,0 }, { 228,0 }, { 242,0 }};
	garkiDatGam.garki304 = { { 18,3 },{ 32,20 },{ 46,19 },{ 60,80 },{ 74,69 },{ 88,357 },{ 102,310 },{ 116,151 },{ 130,141 },{ 144,129 },{ 158,15 },{ 172,2 },{ 186,1 },{ 197,0 },{ 211,0 } , { 225,0 }, { 239,0 }, { 253,0 }, { 267,0 } };
	garkiDatGam.garki553 = { { 12,112 },{ 26,143 },{ 40,656 },{ 54,250 },{ 68,673 },{ 82,683 },{ 96,336 },{ 110,124 },{ 124,19 },{ 138,2 },{ 153,7 },{ 165,13 },{ 179,1 },{ 193,0 } , { 207,0 }, { 221,0 } };
	garkiDatGam.garki802 = { { 18,69 },{ 32,285 },{ 46,66 },{ 60,20 },{ 74,103 },{ 88,67 },{ 102,6 },{ 116,8 },{ 130,3 },{ 144,2 },{ 158,0 },{ 171,2 } , { 185,0 }};


	//modParms initParams;
	initParams.uoE = 0.0342128;
	initParams.uoL = 0.0380283;
	initParams.uP = 0.243781;
	initParams.uM = 0.0887209;
	initParams.Y = 9.64374;
	initParams.w = 175;
	initParams.n = 1.80589;

	initParams.z1 = 3.81073;
	initParams.z2 = 2.84467;
	initParams.z3 = 3.79285;
	initParams.z4 = 2.3262;
	initParams.z5 = 4.01179;
	initParams.z6 = 4.1425;

	initParams.sf1 = 5.14794;
	initParams.sf2 = 4.87898;
	initParams.sf3 = 5.17532;
	initParams.sf4 = 5.59073;
	initParams.sf5 = 4.93854;
	initParams.sf6 = 3.82135;

	initParams.dE = 0.166911;
	initParams.dL = 0.634255;
	initParams.dP = 0.778013;
	initParams.o = 5.87151;
	initParams.dt = 0.25;
	initParams.Mg = 5.56342;

	initParams.p = 0.01;

	initParams.rF = txtReader("Q:\\Imperial\\lModCpp\\Data\\rf01.txt", 0.25);


	vector<double> sdProps = { 
		0.001, 0.001, 0.01,0.01,1,0, 
		0.1, 0.1, 0.1,0.1,0.1,0.1,
		0.1, 0.1, 0.1, 0.1, 0.1,0.1,
		1,0.05,0.05,0.05,1,5,0};

	vector<double> maxSdProps = {
		0.05, 0.05, 0.8, 0.5,5,0,
		0.5, 0.5, 0.5, 0.5,0.5,0.5,
		0.5, 0.5, 0.5, 0.5,0.5,0.5,
		0.1,0.1,0.1,0.1,1,10,0};
	vector<double> acptRs = {
		0.2,0.2,0.2,0.2,0.2,0.2,
		0.2,0.2,0.2,0.2,0.2,0.2,
		0.2,0.2,0.2,0.2,0.2,0.2,
		0.2,0.2,0.2,0.2,0.2,0.2,0.2 };


	//vector<double> sdProps2 = {
	//	0, 0, 0,0,0,0,
	//	0, 0, 0,0,0,0,
	//	0, 0, 0, 0, 0,0,
	//	0,0,0,0,0,0,0 };

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
		300000,//iterations
		150,//particles
		75000,//nburn 
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







