#include<iostream>
#include"lModH.h"

#include <iostream>
#include <map>
#include <random>


int main()
{




	vector<int> rainfall1 = txtReader("C:\\Imperial\\rf1.txt");
	vector<int> rainfall2 = txtReader("C:\\Imperial\\rf2.txt");


	modParms prms;

	prms.rF = rainfall1;
	//vector<int> rainfall1 = txtReader("C:\\Imperial\\rf1.txt");
	//prms.rF = rainfall1;
	vector<tuple<int, int>>garki101 = { { 353,1 },{ 368,18 },{ 382,3 },{ 396,3 },{ 410,7 },{ 424,24 },{ 440,40 },{ 454,1 },{ 468,2 },{ 482,0 },{ 496,0 },{ 510,0 },{ 524,0 },{ 538,0 },{ 552,0 } };
	//vector<int> times;

	//vector<tuple<int, double>> foo = { {1,0.4},{2,0.8},{ 3,0.01 },{ 4,0.3 },{ 5,0.5 },{ 6,0.9 } };


	//	std::random_device rd;
	//	std::mt19937 gen(rd());
	//	std::discrete_distribution<> d({ .40, .10, .10, .40 });
	//	std::map<int, int> m;
	//	for (int n = 0; n<10000; ++n) {
	//		++m[d(gen)];
	//	}
	//	for (auto p : m) {
	//		std::cout << p.first << " generated " << p.second << " times\n";
	//	}


	//	vector<double> lltemp = { 1.0,22.0,34.5,23.536345 };

	//	double llMean = (boost::accumulate(lltemp, 0.0))/lltemp.size();// / lltemp.size;

	//	cout << llMean;
	//	cin.get();
	//for (auto i = begin(garki101); i != end(garki101); ++i) { //calculate discreet time steps for observed data
	//	times.push_back (get<0>(*i));// prms.dt
	//};

	//cout << size(times);
	//cin.get();
	int t = 0;
	int ff = 100;
	while (t < ff) {
		cout << pFilt(50,
			garki101,//garki data
			prms,//parameters
			false,//full output or just likelihood
			1,//rainfall cluster
			7//fixed parameters
		);
		t++;
	}
		cin.get();
	//wpStruct weights;
	//weights.E0 = 1000;
	//weights.L0 = 200;
	//weights.P0 = 50;
	//weights.M0 = 100;
	//weights.startTime = 341;
	//weights.endTime = 552;

	//tuple<int, int, int, int> res = modStepFnc(weights);

	//vector<tuple<int, int, int, int>> res = { { 10, 12, 3, 4 } };

	// res.resize(6, { 10, 12, 3, 4 });



	//std::vector<tuple<int>> z = { {4},{9} };
	//int f = z.back;
	//cout << res.size();
	//cin.get();

	//for (std::vector<int>::const_iterator i = times.begin(); i != times.end(); ++i)
	//	std::cout << *i << ',';
	//cin.get();



	}





		
	
