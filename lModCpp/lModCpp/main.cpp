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
	garkiDat.garki154 = { { 15,3 },{ 29,5 },{ 43,4 },{ 57,6 },{ 71,18 },{ 85,7 },{ 99,21 },{ 113,36 },{ 141,24 },{ 155,5 },{ 169,4 },{ 183,0 },{ 196,0 },{ 210,0 },{ 224,0 },{ 238,1 },{ 252,0 },{ 266,0 } };
	garkiDat.garki202 = { { 8,18 },{ 22,16 },{ 36,20 },{ 50,168 },{ 64,158 },{ 78,377 },{ 92,191 },{ 106,188 },{ 120,96 },{ 134,28 },{ 147,6 },{ 162,4 },{ 175,0 },{ 189,0 },{ 203,1 },{ 216,0 },{ 231,0 },{ 245,1 } };

	garkiDat.garki218 = {{ 19,107 },{ 33,271 },{ 47,212 },{ 61,312 },{ 75,492 },{ 89,326 },{ 103,304 },{ 117,74 },{ 131,15 },{ 147,6 },{ 156,11 },{ 172,0 },{ 186,0 },{ 200,0 } ,{ 215,0 } ,{ 228,0 } ,{ 242,0 } };	

	
	garkiDat.garki304 = { { 18,9 },{ 32,21 },{ 46,23 },{ 60,88 },{ 74,72 },{ 88,358 },{ 102,325 },{ 116,165 },{ 130,224 },{ 144,160 },{ 158,15 },{ 172,2 },{ 186,1 },{ 197,0 },{ 211,0 } ,{ 225,0 } ,{ 239,0 } ,{ 253,0 },{ 267,0 } };	
	
	garkiDat.garki553 = { { 12,120 },{ 26,145 },{ 40,666 },{ 54,253 },{ 68,685 },{ 82,703 },{ 96,379 },{ 110,148 },{ 124,23 },{ 138,5 },{ 153,8 },{ 165,17 },{ 179,1 },{ 193,0 },{ 207,0 },{ 221,1 },{ 235,0 } };//,{ 249,5 },{ 263,4 } , { 276,1 }, { 291,0 }  , { 304,6 }, { 319,1 }, { 333,5 }, { 347,23 }, { 362,37 }, { 376,116 }, { 390,128 }, { 404,45 }, { 418,59 }, { 434,150 }, { 448,160 }, { 462,44 }, { 476,46 }, { 490,3 }, { 504,1 }, { 518,1 }, { 532,0 }, { 546,0 }, { 560,1 }, { 574,0 }, { 588,0 }, { 602,1 }};
	garkiDat.garki802 = { { 18,41 },{ 32,260 },{ 46,140 },{ 60,89 },{ 74,295 },{ 88,233 },{ 102,114 },{ 116,28 },{ 130,28 },{ 144,10 },{ 158,3 },{ 171,5 },{ 185,0 },{ 199,0 } ,{ 212,3 } ,{ 227,0 } ,{ 241,3 } };// ,{ 255,0 } ,{ 269,6 }  ,{ 283,14 } ,{ 297,7 } ,{ 311,0 } } , { 325,5 }, { 339,10 }, { 353,8 }, { 368,85 }, { 382,26 }, { 396,26 }, { 410,51 }, { 424,74 }, { 440,98 }, { 454,45 }, { 468,30 }, { 482,23 }, { 496,0 }, { 510,0 }, { 524,0 }, { 538,1 }, { 552,1 }, { 566,2 }, { 580,0 }, { 594,0 }, { 608,0 }};


	//modParms initParams;
	initParams.uoE = 0.0263704;
	initParams.uoL = 0.00174293;
	initParams.uP = 0.211117;
	initParams.Y = 11.5012;
	initParams.z1 = 3.47013;
	initParams.z2 = 5.0887;
	initParams.z3 = 5.97176;
	initParams.z4 = 1.71979;
	initParams.z5 = 5.72207;
	initParams.z6 = 5.71288;

	initParams.w = 0.0553242;
	initParams.sf1 = 542364;
	initParams.sf2 = 88225.7;
	initParams.sf3 = 529786;
	initParams.sf4 = 1.0527e+06;
	initParams.sf5 = 1.0577e+06;
	initParams.sf6 = 642920;


	initParams.dE = 0.188835;
	initParams.dL = 0.366656;
	initParams.dP = 2.25097;


	initParams.n = 21.1964;
	initParams.rF = txtReader("C:\\Imperial\\lModCpp\\Data\\rf05.txt");

	vector<double> sdProps = { 0.01, 0.01, 0.1, 1, 0.001, 1, 1, 1,1,1,1,1, 10000, 10000, 10000, 10000,10000,10000,0.05,0.05,0.1 };
	vector<double> maxSdProps = { 0.1, 0.1, 0.1, 20, 0.1, 5, 1, 1, 1, 1,1,1, 1e+5, 1e+5, 1e+5, 1e+5,1e+5,1e+5,0.1,0.1,2 };
	vector<double> acptRs = { 0.3,0.3,0.3,0.3,0.3,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.6,0.6,0.6,0.6,0.6,0.6,0.5,0.5,0.5, };
	vector<tuple<string, double>> fitPrms = { { "uoE", initParams.uoE },{ "uoL", initParams.uoL },{ "uP", initParams.uP },{ "Y", initParams.Y },
	{ "w", initParams.w },{ "n", initParams.n },{ "z1", initParams.z1 },{ "z2", initParams.z2 },{ "z3", initParams.z3 },{ "z4", initParams.z4 },{ "z5", initParams.z5 },{ "z6", initParams.z6 },
	{ "sf1", initParams.sf1 } ,{ "sf2", initParams.sf2 } ,{ "sf3", initParams.sf3 } ,{ "sf4", initParams.sf4 } ,{ "sf5", initParams.sf5 } ,{ "sf6", initParams.sf6 } 
	,{ "dE", initParams.dE },{ "dL", initParams.dL } ,{ "dP", initParams.dP } };



	pMMHres results = pMMHSampler(
		initParams,//initial parameters
		7,//fixed parameter for days of rainfall   
		sdProps,//initial sf for param proposals
		acptRs,//acceptance ratios
		fitPrms,//tuple of initial parm values plus names - needed as no reflection - maybe can be coded better
		maxSdProps,//max sd for each parameter proposal in tuner
		500000,//iterations
		45,//particles
		50000,//nburn 
		1,//monitoring
		1000,//start adapt
		200,//tell
		garkiDat//observed data
	);
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << std::endl;

	//write results to csv
	ofstream myfile;
	myfile.open("C:\\Imperial\\1milTestTr7.txt");
	for (auto iter = 0; iter != size(results.ll); ++iter) {
		myfile << results.uoE.at(iter) << " " << results.uoL.at(iter) << " " << results.uP.at(iter) << " " << results.Y.at(iter) << " " << results.w.at(iter)
			<< " " << results.n.at(iter) << " " << results.z1.at(iter) << " " << results.z2.at(iter) << " " << results.z3.at(iter) << " " << results.z4.at(iter) << " " << results.z5.at(iter) << " " << results.z6.at(iter)
			<< " " << results.sf1.at(iter) << " " << results.sf2.at(iter) << " " << results.sf3.at(iter) << " " << results.sf4.at(iter) << " " << results.sf5.at(iter) << " " << results.sf6.at(iter) << 
			" " << results.dE.at(iter) << " " << results.dL.at(iter) << " " << results.dP.at(iter) << " " << results.ll.at(iter) << std::endl;
	}
	cout << "end";
	cin.get();

}







