#include"lModH.h"



//text file reader
//@param file text file location
//@return vector of ints from text file
vector<double> txtReader(string file, float dt) {
	vector<double>dat;
	// open file    
	ifstream inputFile(file);
	if (inputFile) {
		double value;
		// read the elements in the file into a vector  
		while (inputFile >> value) {
			for (int j = 0; j < (1 / dt); j++) {
				dat.emplace_back(value*dt);

			}
		}
	}
	return (dat);
}

//text file reader for mosquito data
//@param file text file location
//@return vector of ints from text file
//vector<tuple<double,double>> txtReaderMos(string file, float dt) {
//	vector<double>dat;
//	// open file    
//	ifstream inputFile(file);
//	if (inputFile) {
//		double value;
//		// read the elements in the file into a vector  
//		while (inputFile >> value) {
//			for (int j = 0; j < (1 / dt); j++) {
//				dat.emplace_back(value*dt);
//
//			}
//		}
//	}
//	return (dat);
//}

double medianFnc(vector<double> vec) {
	double medPos;
	double medPos2;
	double res;
	int n = boost::size(vec);
	//cout << "vec size = " << n << endl;
	std::sort(vec.begin(), vec.end());

	if (n % 2 == 0) {
		medPos = (n / 2) + 1;
		res = vec.at(medPos);
	}
	else {
		medPos = (n / 2);
		medPos2 = (n / 2) + 1;
		res = (vec.at(medPos)+ vec.at(medPos2))/2.0;
	}
	return res;
}



////plot fit function
double pFitFunc(int particles, pMMHres results, obsDatX obsDat, int fixedParam, modParms prms) {
	string outputFile;

	vector<double> rainfall_05 = txtReader("Q:\\Imperial\\lModCpp\\Data\\rf05.txt", 0.25);
	vector<double> rainfall_07 = txtReader("Q:\\Imperial\\lModCpp\\Data\\rf07.txt", 0.25);
	vector<double> rainfall_08 = txtReader("Q:\\Imperial\\lModCpp\\Data\\rf08.txt", 0.25);
	vector<double> rainfall_04 = txtReader("Q:\\Imperial\\lModCpp\\Data\\rf04.txt", 0.25);
	vector<double> rainfall_02 = txtReader("Q:\\Imperial\\lModCpp\\Data\\rf02.txt", 0.25);
	vector<double> rainfall_01 = txtReader("Q:\\Imperial\\lModCpp\\Data\\rf01.txt", 0.25);

	prms.uoE = medianFnc(results.uoE);
	prms.uoL = medianFnc(results.uoL);
	prms.uM = medianFnc(results.uM);

	prms.uP = medianFnc(results.uP);
	prms.Y = medianFnc(results.Y);
	prms.z1 = medianFnc(results.z1);
	prms.z2 = medianFnc(results.z2);
	prms.z3 = medianFnc(results.z3);
	prms.z4 = medianFnc(results.z4);
	prms.z5 = medianFnc(results.z5);
	prms.z6 = medianFnc(results.z6);
	prms.w = medianFnc(results.w);
	prms.sf1 = medianFnc(results.sf1);
	prms.sf2 = medianFnc(results.sf2);
	prms.sf3 = medianFnc(results.sf3);
	prms.sf4 = medianFnc(results.sf4);
	prms.sf5 = medianFnc(results.sf5);
	prms.sf6 = medianFnc(results.sf6);
	prms.dE = medianFnc(results.dE);
	prms.dL = medianFnc(results.dL);
	prms.dP = medianFnc(results.dP);
	prms.o = medianFnc(results.o);
	prms.n = medianFnc(results.n);


	for (auto j = 0; j != 6; ++j) {
		vector<tuple<int, int>> oDat;
		if (j == 0) {
			oDat = obsDat.garki154;
			prms.sf = pow(10, prms.sf1);
			prms.z = pow(10, prms.z1);
			prms.rF = rainfall_05;
			outputFile = "Q:\\Imperial\\fitPlots\\garki154.txt";
		}
		else if (j == 1) {
			oDat = obsDat.garki202;
			prms.sf = pow(10, prms.sf2);
			prms.z = pow(10, prms.z2);
			prms.rF = rainfall_04;
			outputFile = "Q:\\Imperial\\fitPlots\\garki202.txt";
		}
		else if (j == 2) {
			oDat = obsDat.garki218;
			prms.sf = pow(10, prms.sf3);
			prms.z = pow(10, prms.z3);
			prms.rF = rainfall_07;
			outputFile = "Q:\\Imperial\\fitPlots\\garki218.txt";
		}
		else if (j == 3) {
			oDat = obsDat.garki304;
			prms.sf = pow(10, prms.sf4);
			prms.z = pow(10, prms.z4);
			prms.rF = rainfall_08;
			outputFile = "Q:\\Imperial\\fitPlots\\garki304.txt";
		}
		else if (j == 4) {
			oDat = obsDat.garki553;
			prms.sf = pow(10, prms.sf5);
			prms.z = pow(10, prms.z5);
			prms.rF = rainfall_02;
			outputFile = "Q:\\Imperial\\fitPlots\\garki553.txt";
		}
		else if (j == 5) {
			oDat = obsDat.garki802;
			prms.sf = pow(10, prms.sf6);
			prms.z = pow(10, prms.z6);
			prms.rF = rainfall_01;
			outputFile = "Q:\\Imperial\\fitPlots\\garki802.txt";
		}

		//run particle filter
		pFilt(particles,
			oDat,//garki data
			prms,//parameters
			true,//full output or just likelihood
			fixedParam,//fixed parameters
			outputFile//name of output file for plots
		);
	}

}