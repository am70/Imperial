// mosquito population model
#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <string>
#include<fstream>


using namespace std;
std::random_device rd;
std::uniform_real_distribution<> dis(0, 1);
std::mt19937 mrand(std::time(0));
typedef long unsigned int luint;


//bernouli trial function
int bernTrial(int n, float p) {
	int res = 0;
	int t = 1;
	while (t<n) {
		if (dis(mrand)<p) res++;
		t++;
	}
	return res;
}

//random poisson draw
luint rpois(luint lambda)
{
	if (lambda > 0) {
		std::poisson_distribution<luint> d(lambda);
		return d(mrand);
	}
	else return (0);
}


//mosquito population model function

std::vector<int> mPmod(float dE, float dL, float dP, float uoE, float uoL, float uP, float uM, float Y, float S,
	float tr, float sf, float dt, float n, float Emax, float E0, float L0, float P0, float M0,
	float time) {
	int t = 0;
	int Be = 0;
	int Bl = 0;
	int Bp = 0;
	int Bm = 0;
	float nt = 0.0;
	float E = E0;
	float L = L0;
	float P = P0;
	float M = M0;
	int K = 100;
	float uE;
	float uL;
	float trx = tr / dt;

	std::vector<int> r;

	while (t < time) {


		K = 100;//if (time <= trx) ? ? ? if (t<=trx) (1+(sf*((1/trx)*(sum(rF[0:(timeX-1)]))))) else (1+(sf*((1/trx)*(sum(rF[(timeX-trx):timeX-1])))))

		uE = uoE*dt*(1 + ((E + L) / (K)));
		uL = uoL*dt*(1 + (Y*(E + L) / (K)));

		if ((dE + uE)*dt<1) Be = bernTrial(E, (dE + uE)*dt); else Be = bernTrial(E, 1);
		if ((dL + uL)*dt<1) Bl = bernTrial(L, (dL + uL)*dt); else Bl = bernTrial(L, 1);

		Bp = bernTrial(P, (dP + uP)*dt);
		Bm = bernTrial(M, uM*dt);
		nt = bernTrial(M, (dt / S));
		E = E - Be + rpois(n*nt);
		L = L - Bl + bernTrial(Be, (dE / (uE + dE)));
		P = P - Bp + bernTrial(Bl, (dL / (uL + dL)));
		M = M + (0.5*(bernTrial(Bp, (dP / (uP + dP))))) - Bm;

		t++;
		r.push_back(M);

	}

	return r;

}





int main()
{
	vector<double>rainfall;
	// open file    
	ifstream inputFile("Q:\Imperial\rainfall.txt");

	// test file open   
	if (inputFile) {
		double value;
		
		// read the elements in the file into a vector  
		while (inputFile >> value) {
		//	getline(inputFile, value, ' ');
			rainfall.push_back(value);
		}
	}



		std::vector<int> z;
		z = mPmod(0.150, 0.269, 1.563,
			0.034, 0.035, 0.25,
			0.096, 13.25, 3, 14, 20, 0.25, 60.0, 93.6, 177, 8, 1, 7, 133);

		for (std::vector<double>::const_iterator i = rainfall.begin(); i != rainfall.end(); ++i)
			std::cout << *i << ' ';
		cin.get();

	}





		
	
