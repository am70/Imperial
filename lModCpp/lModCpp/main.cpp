#include<iostream>
#include"lModH.h"
#include <iostream>
#include <map>
#include <random>
#include <chrono>
#include <fstream>



int main()
{

	modParms initParms;
	obsDatX garkiDatGam;

	initParms = initParamsReader("Q:\\Imperial\\particleTestLinearN\\1.5milResults\\exp\\medianParms.txt", initParms);
	cin.get();


// read in garki village data anoph gamb 
	garkiDatGam.garki154 = { { 15,3 },{ 29,4 },{ 43,4 },{ 57,6 },{ 71,17 },{ 85,7 },{ 99,20 },{ 113,34 },{ 141,18 },{ 155,5 },{ 169,4 },{ 183,0 },{ 196,0 } , { 210,0 }, { 224,0 }, { 238,1 }, { 252,0 }, { 266,0 }};

	garkiDatGam.garki408 = { { 5,7 },{ 19,99 },{ 33,670 },{ 47,1266 },{ 61,876 },{ 75,1340 },{ 89,227 },{ 103,513 },{ 117,49 },{ 131,19 },{ 145,15 },{ 158,7 } ,{ 172,15 },{186,2} };

	garkiDatGam.garki202 = { { 8,18 },{ 22,15 },{ 36,20 },{ 50,167 },{ 64,155 },{ 78,370 },{ 92,178 },{ 106,127 },{ 120,75 },{ 134,25 },{ 147,5 },{ 162,3 },{ 175,0 } , { 189,0 }, { 203,1 }, { 216,0 }, { 231,0 }, { 245,0 } };
	garkiDatGam.garki218 = { { 19,99 },{ 33,266 },{ 47,211 },{ 61,301 },{ 75,475 },{ 89,309 },{ 103,279 },{ 117,52 },{ 131,14 },{ 147,6 },{ 156,10 }, { 172,0 }, { 186,0 }, { 200,0 }, { 215,0 }, { 228,0 }, { 242,0 }};
	garkiDatGam.garki304 = { { 18,3 },{ 32,20 },{ 46,19 },{ 60,80 },{ 74,69 },{ 88,357 },{ 102,310 },{ 116,151 },{ 130,141 },{ 144,129 },{ 158,15 },{ 172,2 },{ 186,1 },{ 197,0 },{ 211,0 } , { 225,0 }, { 239,0 }, { 253,0 }, { 267,0 } };
	garkiDatGam.garki553 = { { 12,112 },{ 26,143 },{ 40,656 },{ 54,250 },{ 68,673 },{ 82,683 },{ 96,336 },{ 110,124 },{ 124,19 },{ 138,2 },{ 153,7 },{ 165,13 },{ 179,1 },{ 193,0 } , { 207,0 }, { 221,0 } };
	
	garkiDatGam.garki802 = { { 4,9 },{ 18,69 },{ 32,285 },{ 46,66 },{ 60,20 },{ 74,103 },{ 88,67 },{ 102,6 },{ 116,8 },{ 130,3 },{ 144,2 },{ 158,0 },{ 171,2 } ,{ 185,0 },{ 199,0 },{ 212,0 },{ 227,0 },{ 241,0 },
	{ 255,0 },{ 269,0 },{ 283, 2 },{ 297,0 },{ 311,2 },{ 325,3 },{ 339,2 },{ 353,7 } ,{ 368,39 },{ 382,17 },{ 396,13 },{ 410,13 },{ 424,21 },{ 440,24 },{ 454,12 },{ 468,0 },{ 482,3 },{ 496,0 },{ 510,0 } };//, { 524,0 }, { 538,0 }, { 552,0 }
	//,{ 566,0 },{ 580,0 },{ 594,0 },{ 608,0 },{ 622,0 },{ 636,0 },{ 650,0 },{ 663,0 },{ 678,0 },{ 692,1 } ,{ 706,4 } ,{ 720,11 },{ 734,6 },{ 748,7 },{ 762,76 },{ 776,221 },{ 790,213 },{ 804,44 },{ 818,35 },{ 832,3 },{ 846,0 } };

	garkiDatGam.garki801 = { { 4,8 },{ 18,36 },{ 32,249 },{ 46,118 },{ 60,69 },{ 74,250 },{ 88,93 },{ 102,68 },{ 116,24 },{ 130,11 },{ 144,5 },{ 158,0 } ,{ 171,0 } ,{ 185,0 },{ 199,0 },{ 212,0 },{ 227,0 },{ 241,0 },
	{ 255,0 },{ 269,1 },{ 283, 1 },{ 297,1 },{ 311,0 },{ 325,3 },{ 339,5 },{ 353,6 } ,{ 368,83 },{ 382,23 },{ 396,24 },{ 410,41 },{ 424,61 },{ 440,63 },{ 454,14 },{ 468,5 },{ 482,4 },{ 496,0 },{ 510,0 } };//, { 524,0 }, { 538,1 }, { 552,0 }
	//,{ 566,0 },{ 580,0 },{ 594,0 },{ 608,0 },{ 622,0 },{ 636,0 },{ 650,0 },{ 663,0 },{ 678,1 },{ 692,1 } ,{ 706,3 } ,{ 720,8 },{ 734,24 },{ 748,25 },{ 762,175 },{ 776,221 },{ 790,777 },{ 804,329 },{ 818,167 },{ 832,21 },{ 846,1 } };
		

		//modParms initParams;
	initParms.uoE = 0.0351098;
	initParms.uoL = 0.0351108;
	initParms.uP = 0.252691;
	initParms.Y = 13.2635;
	initParms.w = 0.00278759;
	initParms.n = 8.47769;
	initParms.z1 = 1.58026;
	initParms.z4 = 2.09397;
	initParms.z5 = 3.5542;
	initParms.z6 = 3.97494;
	initParms.sf1 = 5.19574;
	initParms.sf4 = 4.26647;
	initParms.sf5 = 5.0198;
	initParms.sf6 = 3.81799;
	initParms.dE = 0.143307;
	initParms.dL = 0.252131;
	initParms.dP = 0.989246;
	initParms.o = 0.608759;
	initParms.tau = 8;
	initParms.uM = 0.0800447;
	initParms.Mg = 1.64751;
	initParms.p = 0.00524687;




	initParms.rF = txtReader("\\\\qdrive.dide.ic.ac.uk\\homes\\ALM210\\Imperial\\lModCpp\\Data\\rf03.txt", initParms.dt);

	
		//initParams.E0 = 24361;
		//initParams.L0 = 244;
		//initParams.P0 = 47;
		//initParams.M0 = 252;

		//initParams.startTime = 5/0.25;
		//initParams.endTime = 186 / 0.25;
		//initParams.z = pow(10, 2.18505);
		//initParams.sf = pow(10, 4.44227);
		//boost::mt19937 rng(std::time(0));
		//mPmod(initParams, rng,"powerClumped");

		vector<double> sdProps = {
			0.001, 0.001, 0.01,0.01,1,0.0001,
			10, 0.1, 0.1,0.1,0.1,
			0.1, 0.1, 0.1,0.1
			,0.01,0.05,0.05,1,0,0.5,0.001};

		vector<double> maxSdProps = {
			0.05, 0.05, 0.8, 0.5,6,0.01,
			10, 5, 5, 5,5,
			7, 7,7,7,
			0.1,0.1,0.1,1,0,4,0.01};

		vector<double> acptRs = {
			0.25,0.25,0.25,0.25,0.25,0.25,
			0.25,0.25,0.25,0.25,
			0.25,0.25,0.25,0.25,
			0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25};


		vector<double> sdProps2 = {
			0, 0,0,0,0,0,
			0, 0, 0,0,0,
			0, 0, 0,0
			,0,0,0,0,0,0,0};

		vector<tuple<string, double>> fitPrms = { { "uoE", initParms.uoE },{ "uoL", initParms.uoL },{ "uP", initParms.uP },{ "uM", initParms.uM },{ "Y", initParms.Y },
		{ "w", initParms.w },{ "n", initParms.n },{ "z1", initParms.z1 },{ "z4", initParms.z4 },{ "z5", initParms.z5 },{ "z6", initParms.z6 },
		{ "sf1", initParms.sf1 }   ,{ "sf4", initParms.sf4 } ,{ "sf5", initParms.sf5 } ,{ "sf6", initParms.sf6 }
		,{ "dE", initParms.dE },{ "dL", initParms.dL } ,{ "dP", initParms.dP } ,{ "o", initParms.o },{ "tau", initParms.tau },{ "Mg", initParms.Mg },{"p",initParms.p } };


		pmcmcOptions pmcmcOpt;
		pmcmcOpt = optionsReader("Q:\\Imperial\\lModCpp\\x64\\powerClumped\\paramOptions.txt", pmcmcOpt);//read in pMCMC options
		string outputFolder = pmcmcOpt.outputFolder; //output folder for results
		string dFunc = pmcmcOpt.dFunc; //Which density/egg laying functions to use: "expClumped", "linearClumped","powerClumped","exp", "linear" or "power"


		pMMHres results = pMMHSampler(
			initParms,//initial parameters
			dFunc,//density function to use: "expClumped", "linearClumped","powerClumped","exp", "linear" or "power"
			sdProps,//initial sf for param proposals
			acptRs,//acceptance ratios
			fitPrms,//tuple of initial parm values plus names - needed as no reflection...
			maxSdProps,//max sd for each parameter proposal in tuner
			pmcmcOpt.iterations,//iterations
			pmcmcOpt.particles,//particles
			pmcmcOpt.nburn,//nburn 
			pmcmcOpt.monitoring,//monitoring
			pmcmcOpt.startAdapt,//start adapt
			pmcmcOpt.tell,//tell
			garkiDatGam//observed data
		);

		//write results to csv
		resultsWriter("\\\\qdrive.dide.ic.ac.uk\\homes\\ALM210\\Imperial\\particleTestLinearN\\",outputFolder, results);
		string fileNameFit = "\\\\qdrive.dide.ic.ac.uk\\homes\\ALM210\\Imperial\\particleTestLinearN\\";
		fileNameFit.append(outputFolder);
		pFitFunc(250, results, garkiDatGam, initParms, fileNameFit, dFunc);
		cout << "End";

	cin.get();
}







