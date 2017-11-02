#include "lModH.h"
#include<iostream>



vector<tuple<int,int,int,int>> mPmod(modParms parmsx) {
	int t = parmsx.startTime;
	int time = parmsx.endTime;
	int Be = 0;
	int Bl = 0;
	int Bp = 0;
	int Bm = 0;
	int nt = 0;
	int E = parmsx.E0;
	int L = parmsx.L0;
	int P = parmsx.P0;
	int M = parmsx.M0;
	int K = 100;
	double trx = parmsx.tr / parmsx.dt;
	double dt = parmsx.dt;
	int rFsum;
	double  uoE = parmsx.uoE;
	double  uoL = parmsx.uoL;
	double  Y = parmsx.Y;
	double  n = parmsx.n;
	double sf = parmsx.sf;
	double dP = parmsx.dP;
	double uP = parmsx.uP;
	int S = parmsx.S;
	double dE = parmsx.dE;
	double dL = parmsx.dL;
	double uM = parmsx.uM;
	double uE;
	double uL;

	vector<int> rF = parmsx.rF;


	vector<tuple<int, int, int, int>> r;

	while (t < time) {
		//cout << "M = " << M << endl;


		if (t <= trx) { K = (1 + (sf*((1 / trx)))); }
		else {
			rFsum = std::accumulate(rF.begin() + (t - trx), rF.begin() + t, 0);
			K = (1 + (sf*((1 / trx)*rFsum)));
		}

		uE = uoE*exp((E + L) / (K));
		uL = uoL*exp((Y*(E + L) / (K)));

	//	cout << (E) << " " << (L) << " " << (P) << " " << (M) << " " << (Be) << " " << (Bl) << " " << (Bp) << " " << (sf) << " " << (K) << endl;
		
		
		Be = binom(E, (dE + uE)*dt);
		 Bl = binom(L, (dL + uL)*dt); 

		// cout << (E) << " " << (dE + uE)*dt << " " << (Be) << endl;


		Bp = binom(P, (dP + uP)*dt);
		Bm = binom(M, uM*dt);
		nt = binom(M, (dt / S));
		E = E - Be + rpois(n*nt);

		L = L - Bl + binom(Be, (dE / (uE + dE)));

		P = P - Bp + binom(Bl, (dL / (uL + dL)));

		if (M + (0.5*(binom(Bp, (dP / (uP + dP))))) - Bm > 0)
			M = M + (0.5*(binom(Bp, (dP / (uP + dP))))) - Bm;
		else M = 1;

		t++;
		r.emplace_back(make_tuple(E,L,P,M));
	}
	return r;

}

