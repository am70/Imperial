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

	
	//vector<tuple<int, int, int, int, double>> test = { {1,2,1,2,-4},{ 2,2,1,2,-4.5 } ,{ 3,2,1,2,-5 } ,{ 4,2,1,2,-3 } };

	//vector<tuple<int, int, int, int, double>> test2 = normalise(test, 100.00);

	//for (auto i = 0; i != size(test2); ++i) {
	//	cout << get<4>(test2[i])<<endl;
	//}


	//vector<tuple<int, int, int, int, double>> test3 = rSamp(test2);

	//for (auto i = 0; i != size(test3); ++i) {
	//	double ff = get<0>(test3[i]);
	//	cout << ff;
	//}

	//cin.get();

		modParms initParams;
	obsDatX garkiDatGam;

// read in garki village data anoph gamb 
	garkiDatGam.garki154 = { { 15,3 },{ 29,4 },{ 43,4 },{ 57,6 },{ 71,17 },{ 85,7 },{ 99,20 },{ 113,34 },{ 141,18 },{ 155,5 },{ 169,4 },{ 183,0 },{ 196,0 } , { 210,0 }, { 224,0 }, { 238,1 }, { 252,0 }, { 266,0 }};

	garkiDatGam.garki408 = { { 5,7 },{ 19,99 },{ 33,670 },{ 47,1266 },{ 61,876 },{ 75,1340 },{ 89,227 },{ 103,513 },{ 117,49 },{ 131,19 },{ 145,15 },{ 158,7 } ,{ 172,15 },{186,2} };

	garkiDatGam.garki202 = { { 8,18 },{ 22,15 },{ 36,20 },{ 50,167 },{ 64,155 },{ 78,370 },{ 92,178 },{ 106,127 },{ 120,75 },{ 134,25 },{ 147,5 },{ 162,3 },{ 175,0 } , { 189,0 }, { 203,1 }, { 216,0 }, { 231,0 }, { 245,0 } };
	garkiDatGam.garki218 = { { 19,99 },{ 33,266 },{ 47,211 },{ 61,301 },{ 75,475 },{ 89,309 },{ 103,279 },{ 117,52 },{ 131,14 },{ 147,6 },{ 156,10 }, { 172,0 }, { 186,0 }, { 200,0 }, { 215,0 }, { 228,0 }, { 242,0 }};
	garkiDatGam.garki304 = { { 18,3 },{ 32,20 },{ 46,19 },{ 60,80 },{ 74,69 },{ 88,357 },{ 102,310 },{ 116,151 },{ 130,141 },{ 144,129 },{ 158,15 },{ 172,2 },{ 186,1 },{ 197,0 },{ 211,0 } , { 225,0 }, { 239,0 }, { 253,0 }, { 267,0 } };
	garkiDatGam.garki553 = { { 12,112 },{ 26,143 },{ 40,656 },{ 54,250 },{ 68,673 },{ 82,683 },{ 96,336 },{ 110,124 },{ 124,19 },{ 138,2 },{ 153,7 },{ 165,13 },{ 179,1 },{ 193,0 } , { 207,0 }, { 221,0 } };
	
	garkiDatGam.garki802 = { { 4,9 },{ 18,69 },{ 32,285 },{ 46,66 },{ 60,20 },{ 74,103 },{ 88,67 },{ 102,6 },{ 116,8 },{ 130,3 },{ 144,2 },{ 158,0 },{ 171,2 } ,{ 185,0 },{ 199,0 },{ 212,0 },{ 227,0 },{ 241,0 },
	{ 255,0 },{ 269,0 },{ 283, 2 },{ 297,0 },{ 311,2 },{ 325,3 },{ 339,2 },{ 353,7 } ,{ 368,39 },{ 382,17 },{ 396,13 },{ 410,13 },{ 424,21 },{ 440,24 },{ 454,12 },{ 468,0 },{ 482,3 },{ 496,0 },{ 510,0 },{ 524,0 } };//, { 538,0 }, { 552,0 }
	//,{ 566,0 },{ 580,0 },{ 594,0 },{ 608,0 },{ 622,0 },{ 636,0 },{ 650,0 },{ 663,0 },{ 678,0 },{ 692,1 } ,{ 706,4 } ,{ 720,11 },{ 734,6 },{ 748,7 },{ 762,76 },{ 776,221 },{ 790,213 },{ 804,44 },{ 818,35 },{ 832,3 },{ 846,0 } };

		garkiDatGam.garki801 = { { 4,8 },{ 18,36 },{ 32,249 },{ 46,118 },{ 60,69 },{ 74,250 },{ 88,93 },{ 102,68 },{ 116,24 },{ 130,11 },{ 144,5 },{ 158,0 } ,{ 171,0 } ,{ 185,0 },{ 199,0 },{ 212,0 },{ 227,0 },{ 241,0 },
		{ 255,0 },{ 269,1 },{ 283, 1 },{ 297,1 },{ 311,0 },{ 325,3 },{ 339,5 },{ 353,6 } ,{ 368,83 },{ 382,23 },{ 396,24 },{ 410,41 },{ 424,61 },{ 440,63 },{ 454,14 },{ 468,5 },{ 482,4 },{ 496,0 },{ 510,0 },{ 524,0 } };//, { 538,1 }, { 552,0 }
	//,{ 566,0 },{ 580,0 },{ 594,0 },{ 608,0 },{ 622,0 },{ 636,0 },{ 650,0 },{ 663,0 },{ 678,1 },{ 692,1 } ,{ 706,3 } ,{ 720,8 },{ 734,24 },{ 748,25 },{ 762,175 },{ 776,221 },{ 790,777 },{ 804,329 },{ 818,167 },{ 832,21 },{ 846,1 } };
	
	

	int ff = 35;
	int tau;

	while (ff <= 35) {
		tau = ff;
	

		//modParms initParams;
		initParams.uoE = 0.0361393;
		initParams.uoL = 0.0295857;
		initParams.uP = 0.279611;
		initParams.uM = 0.0842671;
		initParams.Y = 12.038;
		initParams.w = 0.0110063;
		initParams.n = 7;

		initParams.z1 = 2.11392;
		initParams.z2 = 3.54433;
		initParams.z3 = 4.03068;
		initParams.z4 = 2.77054;
		initParams.z5 = 3.42122;
		initParams.z6 = 3.14363;

		initParams.sf1 = 3.77229;
		initParams.sf2 = 2.81513;
		initParams.sf3 = 5.95307;
		initParams.sf4 = 2.81513;
		initParams.sf5 = 3.62676;
		initParams.sf6 = 2.26616;

		initParams.dE = 0.14424;
		initParams.dL = 0.213182;
		initParams.dP = 0.972388;
		initParams.o = 0.579851;
		initParams.dt = 0.25;
		initParams.Mg = 2.5;
		initParams.tau = 7;

		initParams.p = 0.0196771;

		initParams.rF = txtReader("Q:\\Imperial\\lModCpp\\Data\\rf01.txt", initParams.dt);


		vector<double> sdProps = {
			0.001, 0.001, 0.01,0.01,1,0.0001,
			1, 0.1, 0.1,0.1,0.1,
			0.1, 0.1, 0.1,0.1
			,0.01,0.05,0.05,0,1,0.5,0.001 };

		vector<double> maxSdProps = {
			0.05, 0.05, 0.8, 0.5,6,0.001,
			1, 5, 5, 5,5,
			5, 5,5,5,
			0.1,0.1,0.1,0,3,0.5,0.001 };
		vector<double> acptRs = {
			0.25,0.25,0.25,0.25,0.25,0.25,
			0.25,0.25,0.25,0.25,
			0.25,0.25,0.25,0.25,
			0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25 };


		//vector<double> sdProps2 = {
		//	0, 0, 0,0,0,0,
		//	0, 0, 0,0,0,0,
		//	0, 0, 0, 0, 0,0,
		//	0,0,0,0,0,0,0 };

		vector<tuple<string, double>> fitPrms = { { "uoE", initParams.uoE },{ "uoL", initParams.uoL },{ "uP", initParams.uP },{ "uM", initParams.uM },{ "Y", initParams.Y },
		{ "w", initParams.w },{ "n", initParams.n },{ "z1", initParams.z1 },{ "z4", initParams.z4 },{ "z5", initParams.z5 },{ "z6", initParams.z6 },
		{ "sf1", initParams.sf1 }  ,{ "sf4", initParams.sf4 } ,{ "sf5", initParams.sf5 } ,{ "sf6", initParams.sf6 }
		,{ "dE", initParams.dE },{ "dL", initParams.dL } ,{ "dP", initParams.dP } ,{ "o", initParams.o },{ "tau", initParams.tau },{ "Mg", initParams.Mg },{"p",initParams.p } };



		pMMHres results = pMMHSampler(
			initParams,//initial parameters
			6,//fixed parameter(s) if any - obsolete? 
			sdProps,//initial sf for param proposals
			acptRs,//acceptance ratios
			fitPrms,//tuple of initial parm values plus names - needed as no reflection - maybe can be coded better
			maxSdProps,//max sd for each parameter proposal in tuner
			1000000,//iterations
			75,//particles
			50000,//nburn 
			1,//monitoring
			1000,//start adapt
			25,//tell
			garkiDatGam//observed data
		);

		//std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
		//std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << std::endl;



		//write results to csv
		string fileName = "Q:\\Imperial\\particleTestLinearN\\";
		string loc = to_string(tau);
		fileName.append(loc);
		//boost::filesystem::create_directory(fileName);
		fileName.append("\\results.txt");
		
		ofstream myfile;
		myfile.open(fileName);

		for (auto iter = 0; iter != size(results.ll); ++iter) {
			myfile << results.uoE.at(iter) << " " << results.uoL.at(iter) << " " << results.uP.at(iter) << " " << results.Y.at(iter) << " " << results.w.at(iter)
				<< " " << results.n.at(iter) << " " << results.z1.at(iter) << " " << results.z4.at(iter) << " " << results.z5.at(iter) << " " << results.z6.at(iter)
				<< " " << results.sf1.at(iter) << " " << results.sf4.at(iter) << " " << results.sf5.at(iter) << " " << results.sf6.at(iter) <<
				" " << results.dE.at(iter) << " " << results.dL.at(iter) << " " << results.dP.at(iter) << " " << results.o.at(iter) << " " << results.tau.at(iter)<<" " << results.uM.at(iter) << " " << results.Mg.at(iter) << " " << results.p.at(iter) << " " << results.ll.at(iter) << endl;
		}

		string fileNameFit = "Q:\\Imperial\\particleTestLinearN\\";
		fileNameFit.append(loc);
		pFitFunc(250, results, garkiDatGam, 20, initParams, fileNameFit);
		cout << "End";

		ff = ff+2;
			
	}
	cin.get();

	//std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

}







