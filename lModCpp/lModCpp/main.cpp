#include<iostream>
#include"lModH.h"
#include <iostream>
#include <map>
#include <random>
#include <chrono>
#include <fstream>



int main()
{


	
	//lay eggs into empty system then simulate forwards, but only simulate one generation forward, have eggs, they mature, mating happens, they lay eggs but the eggs they lay don't go further, keep track of how many laid, run few thousand times with set of parameters will give an estimate of Rm.
	//Same can be done for continuous model - starting position you have single female, one they lay clustered eggs, larvae compete, then you get adults again...can do this for both models.run lots of simulations...

	modParms initParms;
	obsDatX garkiDatGam;

// read in garki village data anoph gamb 
	garkiDatGam.garki154 = { { 15,3 },{ 29,4 },{ 43,4 },{ 57,6 },{ 71,17 },{ 85,7 },{ 99,20 },{ 113,34 },{ 141,18 },{ 155,5 },{ 169,4 },{ 183,0 },{ 196,0 } , { 210,0 }, { 224,0 }, { 238,1 }, { 252,0 }, { 266,0 }};
	garkiDatGam.garki408 = { { 5,7 },{ 19,99 },{ 33,670 },{ 47,1266 },{ 61,876 },{ 75,1340 },{ 89,227 },{ 103,513 },{ 117,49 },{ 131,19 },{ 145,15 },{ 158,7 } ,{ 172,15 },{186,2} };
	garkiDatGam.garki202 = { { 8,18 },{ 22,15 },{ 36,20 },{ 50,167 },{ 64,155 },{ 78,370 },{ 92,178 },{ 106,127 },{ 120,75 },{ 134,25 },{ 147,5 },{ 162,3 },{ 175,0 } , { 189,0 }, { 203,1 }, { 216,0 }, { 231,0 }, { 245,0 } };
	garkiDatGam.garki218 = { { 19,99 },{ 33,266 },{ 47,211 },{ 61,301 },{ 75,475 },{ 89,309 },{ 103,279 },{ 117,52 },{ 131,14 },{ 147,6 },{ 156,10 }, { 172,0 }, { 186,0 }, { 200,0 }, { 215,0 }, { 228,0 }, { 242,0 }};
	garkiDatGam.garki304 = { { 18,3 },{ 32,20 },{ 46,19 },{ 60,80 },{ 74,69 },{ 88,357 },{ 102,310 },{ 116,151 },{ 130,141 },{ 144,129 },{ 158,15 },{ 172,2 },{ 186,1 },{ 197,0 },{ 211,0 } , { 225,0 }, { 239,0 }, { 253,0 }, { 267,0 } };
	garkiDatGam.garki553 = { { 12,112 },{ 26,143 },{ 40,656 },{ 54,250 },{ 68,673 },{ 82,683 },{ 96,336 },{ 110,124 },{ 124,19 },{ 138,2 },{ 153,7 },{ 165,13 },{ 179,1 },{ 193,0 } , { 207,0 }, { 221,0 },
	{ 235,0 },{ 249,1 },{ 263,2 },{ 276,1 },{ 291,0 },{ 304,3 },{ 319,0 },{ 333,4 },{ 347,17 },{ 362,33 },{ 376,110 },{ 390,122 },{ 404,44 },{ 418,48 } ,{ 434,144 },{ 448,139 },{ 462,20 },{ 476,13 } ,{ 490,0 },{ 504,1 },
	{ 518, 0 }, { 532,0 }, { 546,0 }, { 560,1 }, { 588,0 }, { 602,0 }, { 616,0 }, { 630,1 }, { 644,1 }, { 658,0 }, { 672,2 }, { 686,0 }, { 700,1 }, { 714,1 }, { 728,3 }, { 742,12 }, { 756,39 }, { 770,171 }, { 784,740 }, { 798,544 } , { 812,305 }, { 826,60 }, { 840,3 }, { 854,1 } };
	garkiDatGam.garki802 = { { 4,9 },{ 18,69 },{ 32,285 },{ 46,66 },{ 60,20 },{ 74,103 },{ 88,67 },{ 102,6 },{ 116,8 },{ 130,3 },{ 144,2 },{ 158,0 },{ 171,2 } ,{ 185,0 },{ 199,0 },{ 212,0 },{ 227,0 },{ 241,0 },
	{ 255,0 },{ 269,0 },{ 283, 2 },{ 297,0 },{ 311,2 },{ 325,3 },{ 339,2 },{ 353,7 } ,{ 368,39 },{ 382,17 },{ 396,13 },{ 410,13 },{ 424,21 },{ 440,24 },{ 454,12 },{ 468,0 },{ 482,3 },{ 496,0 },{ 510,0 }, { 524,0 }, { 538,0 }, { 552,0 }
	,{ 566,0 },{ 580,0 },{ 594,0 },{ 608,0 },{ 622,0 },{ 636,0 },{ 650,0 },{663,0},{678,0},{ 692,1 } ,{ 706,4 } ,{ 720,11 },{ 734,6 },{ 748,7 },{ 762,76 },{ 776,221 },{ 790,213 },{ 804,44 },{ 818,35 },{ 832,3 },{ 846,0 } };
	garkiDatGam.garki801 = { { 4,8 },{ 18,36 },{ 32,249 },{ 46,118 },{ 60,69 },{ 74,250 },{ 88,93 },{ 102,68 },{ 116,24 },{ 130,11 },{ 144,5 },{ 158,0 } ,{ 171,0 } ,{ 185,0 },{ 199,0 },{ 212,0 },{ 227,0 },{ 241,0 },
	{ 255,0 },{ 269,1 },{ 283, 1 },{ 297,1 },{ 311,0 },{ 325,3 },{ 339,5 },{ 353,6 } ,{ 368,83 },{ 382,23 },{ 396,24 },{ 410,41 },{ 424,61 },{ 440,63 },{ 454,14 },{ 468,5 },{ 482,4 },{ 496,0 },{ 510,0 }, { 524,0 }, { 538,1 }, { 552,0 }
	,{ 566,0 },{ 580,0 },{ 594,0 },{ 608,0 },{ 622,0 },{ 636,0 },{ 650,0 },{ 663,0 },{ 678,1 },{ 692,1 } ,{ 706,3 } ,{ 720,8 },{ 734,24 },{ 748,25 },{ 762,175 },{ 776,221 },{ 790,777 },{ 804,329 },{ 818,167 },{ 832,21 },{ 846,1 } };
		

		vector<double> sdProps = {
			0.001, 0.001, 0.01,0.01,1,0.0001,
			10, 0.1, 0.1,0.1,0.1,0.1,0.1,
			0.1, 0.1, 0.1,0.1,0.1,0.1
			,0.01,0.05,0.05,1,1,0.5,0.001};

		vector<double> maxSdProps = {
			0.05, 0.05, 0.8, 0.5,6,0.01,
			10, 0, 0, 0,5,5,5,
			0, 0,7,7,7,7,
			0.1,0.1,0.1,1,1,4,0.01};

		vector<double> acptRs = {
			0.25,0.25,0.25,0.25,0.25,0.25,
			0.25,0.25,0.25,0.25,0.25,0.25,
			0.25,0.25,0.25,0.25,0.25,0.25,
			0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25};

		vector<double> sdProps2 = {
			0, 0,0,0,0,0,
			0, 0, 0,0,0,
			0, 0, 0,0
			,0,0,0,0,0,0,0};


		//Read in pMCMC options from text file to enable multiple instances to run on the cluster with differing setups
		pmcmcOptions pmcmcOpt;
		pmcmcOpt = optionsReader("\\\\qdrive.dide.ic.ac.uk\\homes\\ALM210\\Imperial\\lModCpp\\x64\\expClumpedLong\\paramOptions.txt", pmcmcOpt);//read in pMCMC options
		string outputFolder = pmcmcOpt.outputFolder; //output folder for results
		string dFunc = pmcmcOpt.dFunc; //Which density/egg laying functions to use: "expClumped", "linearClumped","powerClumped","expNoClumped", "linearNoClumped" or "powerNoClumped"

		string initParamsLoc = pmcmcOpt.initParamsLoc; //location of initial parameter values - taken from previous pMCMC run and found in pMCMC options text file
		initParms = initParamsReader(initParamsLoc, initParms);//initial parameters, taken from previous runs

		vector<tuple<string, double>> fitPrms = { { "uoE", initParms.uoE },{ "uoL", initParms.uoL },{ "uP", initParms.uP },{ "uM", initParms.uM },{ "Y", initParms.Y },
		{ "w", initParms.w },{ "n", initParms.n },{ "z1", initParms.z1 },{ "z2", initParms.z2 },{ "z3", initParms.z3 },{ "z4", initParms.z4 },{ "z5", initParms.z5 },{ "z6", initParms.z6 },
		{ "sf1", initParms.sf1 } ,{ "sf2", initParms.sf2 } ,{ "sf3", initParms.sf3 }   ,{ "sf4", initParms.sf4 } ,{ "sf5", initParms.sf5 } ,{ "sf6", initParms.sf6 }
		,{ "dE", initParms.dE },{ "dL", initParms.dL } ,{ "dP", initParms.dP } ,{ "o", initParms.o },{ "tau", initParms.tau },{ "Mg", initParms.Mg },{"p",initParms.p } };

		pMMHres results = pMMHSampler(
			initParms,//initial parameters
			dFunc,//density function to use: "expClumped", "linearClumped","powerClumped","expNoClumped", "linearNoClumped" or "powerNoClumped"
			sdProps,//initial sd for param proposals
			acptRs,//acceptance ratios
			fitPrms,//tuple of initial parm values plus names - needed as no reflection...
			maxSdProps,//max sd for each parameter proposal in tuner
		500000,//iterations
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







