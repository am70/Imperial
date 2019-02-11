#include"lModH.h"



//text file reader for rainfall
//@param file text file location
//@param dt discrete time period for division of daily rainfall
//@return vector of ints from text file
vector<double> txtReader(string file, double dt) {
	vector<double>dat;
	// open file    
	ifstream inputFile(file);
	if (inputFile.good() != TRUE) { cout << "Incorrect rainfall file name"; cin.get(); }

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




//text file reader for pmcmc options
//@param file text file location
//@param pmcmcOpt
//@pmcmcOpt
pmcmcOptions optionsReader(string optionsTextFile, pmcmcOptions pmcmcOpt) {

	try {
		std::ifstream in(optionsTextFile);
		if (in.good() != TRUE) { cout << "Incorrect pmcmcOptions file name"; cin.get(); }

		for (std::string parameter, value;
			std::getline(in, parameter, ',') && std::getline(in, value);
			)

		{
			cout << "Parameter " << parameter << " = " <<value<<endl;
			if (parameter == "outputFolder") pmcmcOpt.outputFolder = value;
			if (parameter == "dFunc") pmcmcOpt.dFunc = value;
			if (parameter == "iterations") pmcmcOpt.iterations = stoi(value);
			if (parameter == "particles") pmcmcOpt.particles = stoi(value);
			if (parameter == "nburn") pmcmcOpt.nburn = stoi(value);
			if (parameter == "monitoring") pmcmcOpt.monitoring = stoi(value);
			if (parameter == "startAdapt") pmcmcOpt.startAdapt = stoi(value);
			if (parameter == "tell") pmcmcOpt.tell = stoi(value);
			if (parameter == "initParams") pmcmcOpt.initParamsLoc = value;

		}
	}
	catch(...)  { cerr << "problem with pMCMC options file, check input file"<<endl; cin.get();}
	return (pmcmcOpt);
}




//text file reader for initial parameters
//@param file text file location
//@param modParms
//@return modParms
modParms initParamsReader(string initParamsTxtFile, modParms initParams) {
	try {
		std::ifstream in(initParamsTxtFile);
		if (in.good() != TRUE) { cout << "Incorrect parameter file name"; cin.get(); }

		for (std::string parameter, value;
			std::getline(in, parameter, '=') && std::getline(in, value);
			)
		{
			value.erase(remove_if(value.begin(), value.end(), isspace), value.end());
			parameter.erase(remove_if(parameter.begin(), parameter.end(), isspace), parameter.end());

			cout << parameter << " = " <<  value << endl;
			if (parameter == "initParms.uoE") initParams.uoE = stod(value);
			if (parameter == "initParms.uoL") initParams.uoL = stod(value);
			if (parameter == "initParms.uP") initParams.uP = stod(value);
			if (parameter == "initParms.Y") initParams.Y = stod(value);
			if (parameter == "initParms.w") initParams.w = stod(value);
			if (parameter == "initParms.n") initParams.n = stod(value);
			if (parameter == "initParms.z1") initParams.z1 = stod(value);
			if (parameter == "initParms.z2") initParams.z2 = stod(value);
			if (parameter == "initParms.z3") initParams.z3 = stod(value);
			if (parameter == "initParms.z4") initParams.z4 = stod(value);
			if (parameter == "initParms.z5") initParams.z5 = stod(value);
			if (parameter == "initParms.z6") initParams.z6 = stod(value);
			if (parameter == "initParms.z7") initParams.z7 = stod(value);
			if (parameter == "initParms.z8") initParams.z8 = stod(value);

			if (parameter == "initParms.sf1") initParams.sf1 = stod(value);
			if (parameter == "initParms.sf2") initParams.sf2 = stod(value);
			if (parameter == "initParms.sf3") initParams.sf3 = stod(value);
			if (parameter == "initParms.sf4") initParams.sf4 = stod(value);
			if (parameter == "initParms.sf5") initParams.sf5 = stod(value);
			if (parameter == "initParms.sf6") initParams.sf6 = stod(value);
			if (parameter == "initParms.sf7") initParams.sf7 = stod(value);
			if (parameter == "initParms.sf8") initParams.sf8 = stod(value);

			if (parameter == "initParms.dE") initParams.dE = stod(value);
			if (parameter == "initParms.dL") initParams.dL = stod(value);
			if (parameter == "initParms.dP") initParams.dP = stod(value);
			if (parameter == "initParms.o") initParams.o = stod(value);
			if (parameter == "initParms.tau") initParams.tau = stod(value);
			if (parameter == "initParms.uM") initParams.uM = stod(value);
			if (parameter == "initParms.Mg") initParams.Mg = stod(value);
			if (parameter == "initParms.p") initParams.p = stod(value);
			if (parameter == "initParms.lK") initParams.lK = stod(value);
			if (parameter == "initParms.lKs") initParams.lKs = stod(value);
			if (parameter == "initParms.lKm") initParams.lKm = stod(value);


		}
		return (initParams);

	}
	catch (...) { cerr << "problem with pMCMC options file, check input file/spellings etc." << endl; cin.get(); }
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
void pFitFunc(int particles, pMMHres results, obsDatX obsDat, modParms prms, string outputFile, string dFunc) {

	string orgFile = outputFile;

	vector<double> rainfall_05 = txtReader("\\\\qdrive.dide.ic.ac.uk\\homes\\ALM210\\Imperial\\lModCpp\\Data\\rf05.txt", prms.dt);
	vector<double> rainfall_07 = txtReader("\\\\qdrive.dide.ic.ac.uk\\homes\\ALM210\\Imperial\\lModCpp\\Data\\rf07.txt", prms.dt);
	vector<double> rainfall_08 = txtReader("\\\\qdrive.dide.ic.ac.uk\\homes\\ALM210\\Imperial\\lModCpp\\Data\\rf08.txt", prms.dt);
	vector<double> rainfall_04 = txtReader("\\\\qdrive.dide.ic.ac.uk\\homes\\ALM210\\Imperial\\lModCpp\\Data\\rf04.txt", prms.dt);
	vector<double> rainfall_02 = txtReader("\\\\qdrive.dide.ic.ac.uk\\homes\\ALM210\\Imperial\\lModCpp\\Data\\rf02.txt", prms.dt);
	vector<double> rainfall_01 = txtReader("\\\\qdrive.dide.ic.ac.uk\\homes\\ALM210\\Imperial\\lModCpp\\Data\\rf01.txt", prms.dt);
	vector<double> rainfall_03 = txtReader("\\\\qdrive.dide.ic.ac.uk\\homes\\ALM210\\Imperial\\lModCpp\\Data\\rf03.txt", prms.dt);
	//vector<double> rainfall_01_2 = txtReader("\\\\qdrive.dide.ic.ac.uk\\homes\\ALM210\\Imperial\\lModCpp\\Data\\rf01_2.txt", prms.dt);
	//vector<double> rainfall_02_2 = txtReader("\\\\qdrive.dide.ic.ac.uk\\homes\\ALM210\\Imperial\\lModCpp\\Data\\rf02_2.txt", prms.dt);



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
	prms.z7 = medianFnc(results.z7);
	prms.z8 = medianFnc(results.z8);
	prms.w = medianFnc(results.w);
	prms.sf1 = medianFnc(results.sf1);
	prms.sf2 = medianFnc(results.sf2);
	prms.sf3 = medianFnc(results.sf3);
	prms.sf4 = medianFnc(results.sf4);
	prms.sf5 = medianFnc(results.sf5);
	prms.sf6 = medianFnc(results.sf6);
	prms.sf7 = medianFnc(results.sf7);
	prms.sf8 = medianFnc(results.sf8);
	prms.dE = medianFnc(results.dE);
	prms.dL = medianFnc(results.dL);
	prms.dP = medianFnc(results.dP);
	prms.o = medianFnc(results.o);
	prms.n = medianFnc(results.n);
	prms.p = medianFnc(results.p);
	prms.tau = medianFnc(results.tau);
	prms.lK = medianFnc(results.lK);
	prms.lKs = medianFnc(results.lKs);
	prms.lKm = medianFnc(results.lKm);



	for (auto j = 0; j != 8; ++j) {
		vector<tuple<int, int>> oDat;
		if (j == 0) {
			oDat = obsDat.garki408;
			prms.sf = pow(10, prms.sf1);
			prms.z = pow(10, prms.z1);
			prms.rF = rainfall_03;
			outputFile.append("\\garki408");
		}
		else if (j == 1) {
			outputFile = orgFile;
			oDat = obsDat.garki154;
			prms.sf = pow(10, prms.sf2);
			prms.z = pow(10, prms.z2);
			prms.rF = rainfall_05;
			outputFile.append("\\garki154");
		}
		 if (j == 2) {
			outputFile = orgFile;
			oDat = obsDat.garki801;
			prms.sf = pow(10, prms.sf3);
			prms.z = pow(10, prms.z3);
			prms.rF = rainfall_01;
			outputFile.append("\\garki801");
		}
		  if (j == 3) {
			outputFile = orgFile;
			oDat = obsDat.garki802;
			prms.sf = pow(10, prms.sf4);
			prms.z = pow(10, prms.z4);
			prms.rF = rainfall_01;
			outputFile.append("\\garki802");
		}
		else if (j == 4) {
			outputFile = orgFile;
			oDat = obsDat.garki553;
			prms.sf = pow(10, prms.sf5);
			prms.z = pow(10, prms.z5);
			prms.rF = rainfall_02;
			outputFile.append("\\garki553");
		}
		
		else if (j == 5) {
			outputFile = orgFile;
			oDat = obsDat.garki801_2;
			prms.sf = pow(10, prms.sf6);
			prms.z = pow(10, prms.z6);
			prms.rF = rainfall_01_2;
			outputFile.append("\\garki801");
		}
		else if (j == 6) {
			outputFile = orgFile;
			oDat = obsDat.garki802_2;
			prms.sf = pow(10, prms.sf7);
			prms.z = pow(10, prms.z7);
			prms.rF = rainfall_01_2;
			outputFile.append("\\garki802");
		}
		else if (j == 7) {
			outputFile = orgFile;
			oDat = obsDat.garki553_2;
			prms.sf = pow(10, prms.sf8);
			prms.z = pow(10, prms.z8);
			prms.rF = rainfall_02_2;
			outputFile.append("\\garki553");
		}

		//run particle filter
		pFilt(particles,
			oDat,//garki data
			prms,//parameters
			true,//full output or just likelihood
			outputFile,//name of output file for plots
			false,
			dFunc
		);
	}

	string medFile = orgFile;
	medFile.append("\\medianParms.txt");
	ofstream myfile;
	myfile.open(medFile);
	myfile << "initParms.uoE = " << prms.uoE << endl << "initParms.uoL  = " << prms.uoL << endl << " initParms.uP = " << prms.uP << endl << "initParms.Y = " << prms.Y << endl << "initParms.w = " << prms.w << endl
		<< "initParms.n = " << prms.n << endl << "initParms.z1 = " << prms.z1 << endl << "initParms.z2 = " << prms.z2 << endl << "initParms.z3 = " << prms.z3 << endl << "initParms.z4 = " << prms.z4 << endl << "initParms.z5 = " << prms.z5 << endl << "initParms.z6 = " << prms.z6 << endl << "initParms.z7 = " << prms.z7 << endl << "initParms.z8 = " << prms.z8 << endl
		<< "initParms.sf1 = " << prms.sf1 << endl << "initParms.sf2 = " << prms.sf2 << endl << "initParms.sf3 = " << prms.sf3 << endl << "initParms.sf4 = " << prms.sf4 << endl << "initParms.sf5 = " << prms.sf5 << endl << "initParms.sf6 = " << prms.sf6 << endl << endl << "initParms.sf7 = " << prms.sf7 << endl << endl << "initParms.sf8 = " << prms.sf8 << endl <<
		"initParms.dE = " << prms.dE << endl << "initParms.dL = " << prms.dL << endl << "initParms.dP = " << prms.dP << endl << "initParms.o = " << prms.o << endl << "initParms.tau = " << prms.tau << endl << "initParms.uM = " << prms.uM << endl << "initParms.Mg = " << prms.Mg << endl << "initParms.p = " << prms.p << endl << "initParms.lK = " << prms.lK << endl
		<< "initParms.lKs = " << prms.lKs << endl << "initParms.lKm = " << prms.lKm << endl;

	orgFile.append("\\endParms.txt");
	ofstream myfile2;
	myfile2.open(orgFile);
	myfile2 << "initParms.uoE = " << prms.uoE << endl << "initParms.uoL  = " << prms.uoL << endl << " initParms.uP = " << prms.uP << endl << "initParms.Y = " << prms.Y << endl << "initParms.w = " << prms.w << endl
		<< "initParms.n = " << prms.n << endl << "initParms.z1 = " << prms.z1 << endl << "initParms.z2 = " << prms.z2 << endl << "initParms.z3 = " << prms.z3 << endl << "initParms.z4 = " << prms.z4 << endl << "initParms.z5 = " << prms.z5 << endl << "initParms.z6 = " << prms.z6 << endl << "initParms.z7 = " << prms.z7 << endl << "initParms.z8 = " << prms.z8 << endl
		<< "initParms.sf1 = " << prms.sf1 << endl << "initParms.sf2 = " << prms.sf2 << endl << "initParms.sf3 = " << prms.sf3 << endl << "initParms.sf4 = " << prms.sf4 << endl << "initParms.sf5 = " << prms.sf5 << endl << "initParms.sf6 = " << prms.sf6 << endl << "initParms.sf7 = " << prms.sf7 << endl << endl << "initParms.sf8 = " << prms.sf8 << endl <<
		"initParms.dE = " << prms.dE << endl << "initParms.dL = " << prms.dL << endl << "initParms.dP = " << prms.dP << endl << "initParms.o = " << prms.o << endl << "initParms.tau = " << prms.tau << endl << "initParms.uM = " << prms.uM << endl << "initParms.Mg = " << prms.Mg << endl << "initParms.p = " << prms.p << endl << "initParms.lK = " << prms.lK << endl
		<< "initParms.lKs = " << prms.lKs << endl << "initParms.lKm = " << prms.lKm << endl;
}


double resultsWriter(string fileName, string outputFolder, pMMHres results) {

	fileName.append(outputFolder);
	fileName.append("\\results.txt");

	ofstream myfile;
	myfile.open(fileName);

	for (auto iter = 0; iter != size(results.ll); ++iter) {
		myfile << results.uoE.at(iter) << " " << results.uoL.at(iter) << " " << results.uP.at(iter) << " " << results.Y.at(iter) << " " << results.w.at(iter)
			<< " " << results.n.at(iter) << " " << results.z1.at(iter) << " " << results.z2.at(iter) << " " << results.z3.at(iter) << " " << results.z4.at(iter) << " " << results.z5.at(iter) << " " << results.z6.at(iter) << " " << results.z7.at(iter) << " " << results.z8.at(iter)
			<< " " << results.sf1.at(iter) << " " << results.sf2.at(iter) << " " << results.sf3.at(iter) << " " << results.sf4.at(iter) << " " << results.sf5.at(iter) << " " << results.sf6.at(iter) << " " << results.sf7.at(iter) << " " << results.sf8.at(iter) <<
			" " << results.dE.at(iter) << " " << results.dL.at(iter) << " " << results.dP.at(iter) << " " << results.o.at(iter) << " " << results.tau.at(iter) << " " << results.uM.at(iter) << " " << results.Mg.at(iter) << " " << results.p.at(iter) << " " << results.lK.at(iter) <<
			" " << results.lKs.at(iter) << " " << results.lKm.at(iter) << " " << results.ll.at(iter) << endl;
	}

	return 0.0;

}