// mosquito population model
#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <string>
#include<fstream>
#include <numeric>
#include <boost/regex.hpp>




using namespace std;

std::uniform_real_distribution<> dis(0, 1);
std::mt19937 mrand(std::time(0));
typedef long unsigned int luint;

//binomial function
int binom(int n, double p) {
	std::binomial_distribution<int> distribution(n, p);
	return  distribution(mrand);
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

std::vector<int> mPmod(double dE, double dL, double dP, double uoE, double uoL, double uP, double uM, double Y, double S,
	double tr, double sf, double dt, double n, double Emax, int E0, int L0, int P0, int M0,
	double time, vector<int> rF) {
	int t = 0;
	int Be = 0;
	int Bl = 0;
	int Bp = 0;
	int Bm = 0;
	int nt = 0;
	int E = E0;
	int L = L0;
	int P = P0;
	int M = M0;
	int K = 100;
	double uE;
	double uL;
	double trx = tr / dt;
	int rFsum;

	

	std::vector<int> r;

	while (t < time) {

		if (t <= trx) { K = (1 + (sf*((1 / trx)))); }
	else {
		rFsum = std::accumulate(rF.begin() + (t - trx), rF.end() - t, 0);
			K =  (1 + (sf*((1 / trx)*rFsum)));
		}


	uE = uoE*dt*(1 + (E + L) / (K));
	uL = uoL*dt*(1 + (Y*(E + L) / (K)));


		if ((dE + uE)*dt<1) Be = binom(E, (dE + uE)*dt); else Be = binom(E, 1);
		if ((dL + uL)*dt<1) Bl = binom(L, (dL + uL)*dt); else Bl = binom(L, 1);

		Bp = binom(P, (dP + uP)*dt);
		Bm = binom(M, uM*dt);
		nt = binom(M, (dt / S));
		E = E - Be + rpois(n*nt);
		L = L - Bl + binom(Be, (dE / (uE + dE)));
		P = P - Bp + binom(Bl, (dL / (uL + dL)));
		M = M + (0.5*(binom(Bp, (dP / (uP + dP))))) - Bm;

		t++;
		r.push_back(M);

	}

	return r;

}





int main()
{
	clock_t tStart = clock();
	/* Do your stuff here */

	vector<int>rainfall;
	// open file    
	ifstream inputFile("Q:\\Imperial\\rf2.txt");

	

	// open file
	if (inputFile) {
		int value;
		
		// read the elements in the file into a vector  
		while (inputFile >> value) {
		//	getline(inputFile, value, ' ');
			rainfall.push_back(value);
		}
	
	}
	//std::vector<int>   sub(&rainfall[50], &rainfall[60]);
	int st = 0;
	int stime = 1000;
	   
	while (st < stime) {
		std::vector<int> z;
		z = mPmod(
			0.150, //dE
			0.269, //dL
			1.563, //dP
			0.034, //uoE
			0.035, //uoL
			0.25, //uP
			0.096, //uM
			13.25, //Y
			3, //S
			7, //tr
			60000.0, //sf
			0.25, //dt
			76.16875, //n
			93.6, //Emax
			39250, //E0
			216, //L0
			106, //P0
			862, //M0
			1056, //time
			rainfall);
		st++;
	}
	printf("Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
	cin.get();
	return 0;

	//	for (std::vector<double>::const_iterator i = z.begin(); i != z.end(); ++i)
	//		std::cout << *i << ',';
	//	cin.get();

		

	}





		
	
