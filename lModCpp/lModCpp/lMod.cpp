// mosquito population model
#include <iostream>
#include <vector>
#include <random>
#include <ctime>
#include <string>
#include<fstream>
#include <numeric>



using namespace std;
std::random_device rd;
std::uniform_real_distribution<> dis(0, 1);
std::mt19937 mrand(std::time(0));
typedef long unsigned int luint;



//binomial function
int binom(int n, double p) {
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

std::vector<double> mPmod(double dE, double dL, double dP, double uoE, double uoL, double uP, double uM, double Y, double S,
	double tr, double sf, double dt, double n, double Emax, double E0, double L0, double P0, double M0,
	double time, vector<int> rF) {
	double t = 0;
	double Be = 0;
	double Bl = 0;
	double Bp = 0;
	double Bm = 0;
	double nt = 0.0;
	double E = E0;
	double L = L0;
	double P = P0;
	double M = M0;
	double K = 100;
	double uE;
	double uL;
	double trx = tr / dt;
	double rFsum;

	

	std::vector<double> r;

	while (t < time) {

		if (t <= trx) { K = (1 + (sf*((1 / trx)))); }
	else {
			vector<int>::const_iterator first = rF.begin() + t - trx;
			vector<int>::const_iterator last = rF.begin() + t;
			vector<int> rFx(first, last);
			rFsum = std::accumulate(rFx.begin(), rFx.end(), 0);
			K = (1 + (sf*((1 / trx)*rFsum)));
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
		r.push_back(Be);

	}

	return r;

}





int main()
{


	vector<int>rainfall;
	// open file    
	ifstream inputFile("C:\\Imperial\\rf1.txt");

	

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
	int stime = 500;
	vector<int> gg;


		std::vector<double> z;
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
			50, //sf
			0.25, //dt
			60, //n
			93.6, //Emax
			177, //E0
			8, //L0
			1, //P0
			7, //M0
			600, //time
			rainfall);

		for (std::vector<double>::const_iterator i = z.begin(); i != z.end(); ++i)
			std::cout << *i << ',';
		cin.get();

		

	}





		
	
