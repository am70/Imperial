#include"lModH.h"

//text file reader
//@param file text file location
//@return vector of ints from text file
vector<int> txtReader(string file) {
	vector<int>dat;
	// open file    
	ifstream inputFile(file);
	if (inputFile) {
		int value;
		// read the elements in the file into a vector  
		while (inputFile >> value) {
			dat.emplace_back(value);
		}
	}
	return (dat);
}

