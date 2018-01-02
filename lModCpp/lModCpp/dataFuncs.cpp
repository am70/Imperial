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