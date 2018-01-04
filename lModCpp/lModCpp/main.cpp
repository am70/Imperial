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

	/*for (int j = 0; j < boost::size(rainfall_04); j++) {
		cout << rainfall_04.at(j) << ",";
	}
	cin.get();*/

	double inf = std::numeric_limits<double>::infinity();
	modParms initParams;
	// read in garki village data gamb and fun
	//obsDatX garkiDatAll;
	obsDatX garkiDatGam;
	//obsDatX garkiDatFun;

	//garkiDatAll.garki154 = { { 15,3 },{ 29,4 },{ 43,4 },{ 57,6 },{ 71,18 },{ 85,7 },{ 99,21 },{ 113,36 },{ 141,24 },{ 155,5 },{ 169,4 },{ 183,0 },{ 196,0 },{ 210,0 },{ 224,0 },{ 238,1 },{ 252,0 },{ 266,0 } };
	//garkiDatAll.garki202 = { { 8,18 },{ 22,15 },{ 36,20 },{ 50,168 },{ 64,158 },{ 78,376 },{ 92,187 },{ 106,178 },{ 120,94 },{ 134,26 },{ 147,6 },{ 162,4 },{ 175,0 },{ 189,0 },{ 203,1 },{ 216,0 },{ 231,0 },{ 245,0 } };
	//garkiDatAll.garki218 = { {5,31} ,{ 19,99 },{ 33,266 },{ 47,212 },{ 61,311 },{ 75,489 },{ 89,323 },{ 103,293 },{ 117,72 },{ 131,15 },{ 147,6 },{ 156,11 },{ 172,0 },{ 186,0 },{ 200,0 } ,{ 215,0 } ,{ 228,0 } ,{ 242,0 } };
	//garkiDatAll.garki304 = { { 18,6 },{ 32,20 },{ 46,20 },{ 60,80 },{ 74,72 },{ 88,357 },{ 102,317 },{ 116,156 },{ 130,158 },{ 144,145 },{ 158,15 },{ 172,2 },{ 186,1 },{ 197,0 },{ 211,0 } ,{ 225,0 } ,{ 239,0 } ,{ 253,0 },{ 267,0 } };	
	//garkiDatAll.garki553 = { { 12,113 },{ 26,144 },{ 40,659 },{ 54,252 },{ 68,685 },{ 82,702 },{ 96,378 },{ 110,145 },{ 124,22 },{ 138,5 },{ 153,8 },{ 165,17 },{ 179,1 },{ 193,0 },{ 207,0 },{ 221,1 } };
	//garkiDatAll.garki802 = { {4,9} ,{ 18,41 },{ 32,260 },{ 46,140 },{ 60,89 },{ 74,295 },{ 88,232 },{ 102,204 },{ 116,108 },{ 130,27 },{ 144,10 },{ 158,2 },{ 171,5 },{ 185,0 },{ 199,0 } ,{ 212,2 } };

// read in garki village data gamb 
	garkiDatGam.garki154 = { { 15,3 },{ 29,4 },{ 43,4 },{ 57,6 },{ 71,17 },{ 85,7 },{ 99,20 },{ 113,34 },{ 141,18 },{ 155,5 },{ 169,4 },{ 183,0 },{ 196,0 } , { 210,0 }, { 224,0 }, { 238,1 }, { 252,0 }, { 266,0 }};
	garkiDatGam.garki202 = { { 8,18 },{ 22,15 },{ 36,20 },{ 50,167 },{ 64,155 },{ 78,370 },{ 92,178 },{ 106,127 },{ 120,75 },{ 134,25 },{ 147,5 },{ 162,3 },{ 175,0 } , { 189,0 }, { 203,1 }, { 216,0 }, { 231,0 }, { 245,0 } };
	garkiDatGam.garki218 = { { 19,99 },{ 33,266 },{ 47,211 },{ 61,301 },{ 75,475 },{ 89,309 },{ 103,279 },{ 117,52 },{ 131,14 },{ 147,6 },{ 156,10 },{ 172,0 },{ 186,0 } , { 200,0 }, { 215,0 }, { 228,0 }, { 242,0 } };
	garkiDatGam.garki304 = { { 18,3 },{ 32,20 },{ 46,19 },{ 60,80 },{ 74,69 },{ 88,357 },{ 102,310 },{ 116,151 },{ 130,141 },{ 144,129 },{ 158,15 },{ 172,2 },{ 186,1 },{ 197,0 },{ 211,0 } , { 225,0 }, { 239,0 }, { 253,0 }, { 267,0 } };
	garkiDatGam.garki553 = { { 12,112 },{ 26,143 },{ 40,656 },{ 54,250 },{ 68,673 },{ 82,683 },{ 96,336 },{ 110,124 },{ 124,19 },{ 138,2 },{ 153,7 },{ 165,13 },{ 179,1 },{ 193,0 } , { 207,0 }, { 221,0 } };
	garkiDatGam.garki802 = { { 18,36 },{ 32,249 },{ 46,118 },{ 60,69 },{ 74,250 },{ 88,93 },{ 102,68 },{ 116,24 },{ 130,11 },{ 144,5 },{ 158,0 },{ 171,0 } , { 185,0 }, { 199,0 }, { 212,0 } };


	//modParms initParams;
	initParams.uoE = 0.0781;
	initParams.uoL = 0.0084;
	initParams.uP = 0.31;
	initParams.Y = 8,37;
	initParams.z1 = 4.309;
	initParams.z2 = 4.367;
	initParams.z3 = 4.08;
	initParams.z4 = 3.64;
	initParams.z5 = 4.68;
	initParams.z6 = 4.14;

	initParams.w = 0.000495889;
	initParams.sf1 = 4.909;
	initParams.sf2 = 7.87;
	initParams.sf3 = 6.17;
	initParams.sf4 = 5.72;
	initParams.sf5 = 5.74;
	initParams.sf6 = 5.079;


	initParams.dE = 0.144;
	initParams.dL = 0.317;
	initParams.dP = 1.4492;
	initParams.o = 11;


	initParams.n = 1.27;
	initParams.rF = txtReader("Q:\\Imperial\\lModCpp\\Data\\rf05.txt",0.25);

	vector<double> sdProps = { 0.1, 0.1, 0.1, 1, 0.001, 10, 2, 2,2,2,2,2, 2, 2, 2, 2,2,2,0.1,0.1,0.1};
	vector<double> maxSdProps = { 0.1, 0.1, 0.1, 2, 0.1, 2, 1, 1, 1, 1,1,1, 1, 1, 1, 1,1,1,0.1,0.1,0.1,3};
	vector<double> acptRs = { 0.3,0.3,0.3,0.3,0.3,0.3,
		0.3,0.3,0.3,0.3,0.3,0.3,
		0.3,0.3,0.3,0.3,0.3,0.3,
		0.3,0.3,0.3,0.3 };
	vector<tuple<string, double>> fitPrms = { { "uoE", initParams.uoE },{ "uoL", initParams.uoL },{ "uP", initParams.uP },{ "Y", initParams.Y },
	{ "w", initParams.w },{ "n", initParams.n },{ "z1", initParams.z1 },{ "z2", initParams.z2 },{ "z3", initParams.z3 },{ "z4", initParams.z4 },{ "z5", initParams.z5 },{ "z6", initParams.z6 },
	{ "sf1", initParams.sf1 } ,{ "sf2", initParams.sf2 } ,{ "sf3", initParams.sf3 } ,{ "sf4", initParams.sf4 } ,{ "sf5", initParams.sf5 } ,{ "sf6", initParams.sf6 } 
	,{ "dE", initParams.dE },{ "dL", initParams.dL } ,{ "dP", initParams.dP } ,{ "o", initParams.o }};



	pMMHres results = pMMHSampler(
		initParams,//initial parameters
		10.7714,//fixed parameter for days of rainfall   
		sdProps,//initial sf for param proposals
		acptRs,//acceptance ratios
		fitPrms,//tuple of initial parm values plus names - needed as no reflection - maybe can be coded better
		maxSdProps,//max sd for each parameter proposal in tuner
		50000,//iterations
		75,//particles
		20000,//nburn 
		1,//monitoring
		3000,//start adapt
		25,//tell
		garkiDatGam//observed data
	);
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << std::endl;

	//write results to csv
	ofstream myfile;
	myfile.open("Q:\\Imperial\\1milTestTr7_M_1.txt");
	for (auto iter = 0; iter != size(results.ll); ++iter) {
		myfile << results.uoE.at(iter) << " " << results.uoL.at(iter) << " " << results.uP.at(iter) << " " << results.Y.at(iter) << " " << results.w.at(iter)
			<< " " << results.n.at(iter) << " " << results.z1.at(iter) << " " << results.z2.at(iter) << " " << results.z3.at(iter) << " " << results.z4.at(iter) << " " << results.z5.at(iter) << " " << results.z6.at(iter)
			<< " " << results.sf1.at(iter) << " " << results.sf2.at(iter) << " " << results.sf3.at(iter) << " " << results.sf4.at(iter) << " " << results.sf5.at(iter) << " " << results.sf6.at(iter) << 
			" " << results.dE.at(iter) << " " << results.dL.at(iter) << " " << results.dP.at(iter) << " "  << results.o.at(iter) << " " << results.ll.at(iter) <<endl;
	}

	pFitFunc(75, results, garkiDatGam, 20, initParams);

	cout << "end";
	cin.get();

}







