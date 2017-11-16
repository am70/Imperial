#include "lModH.h"

vector<tuple<int, int, int, int>> mPmod(modParms parmsx, boost::mt19937 rd) {

	int t = parmsx.startTime;
	int time = parmsx.endTime;
	double Be = 0;
	double Bl = 0;
	double Bp = 0;
	double Bm = 0;
	double nt = 0;
	double E = parmsx.E0;
	double L = parmsx.L0;
	double P = parmsx.P0;
	double M = parmsx.M0;
	double K = 100;
	double trx = parmsx.tr / parmsx.dt;
	double dt = parmsx.dt;
	double rFsum;
	double  uoE = parmsx.uoE;
	double  uoL = parmsx.uoL;
	double  Y = parmsx.Y;
	double  n = parmsx.n;
	double sf = parmsx.sf;
	double dP = parmsx.dP;
	double uP = parmsx.uP;
	double S = parmsx.S;
	double dE = parmsx.dE;
	double dL = parmsx.dL;
	double uM = parmsx.uM;
	double uE;
	double uL;
	double mRan;

	vector<int> rF = parmsx.rF;
	vector<tuple<int, int, int, int>> r;

	while (t < time) {

		if (t <= trx) { K = (1 + (sf*((1 / trx)))); }
		else {
			rFsum = std::accumulate(rF.begin() + (t - trx), rF.begin() + t, 0);
			K = (1 + (sf*((1 / trx)*rFsum)));
		}

		uE = uoE*exp((E + L) / (K));
		uL = uoL*exp((Y*(E + L) / (K)));

		if ((dE + uE)*dt < 1) {
			boost::binomial_distribution<int> distributionBe(E, (dE + uE)*dt);
			Be = distributionBe(rd);
		}
		else {
			boost::binomial_distribution<int> distributionBe(E,1);
			Be = distributionBe(rd);
		}

		if ((dL + uL)*dt < 1) {
			boost::binomial_distribution<int> distributionBl(L, (dL + uL)*dt);
			Bl = distributionBl(rd);
		}
		else {
			boost::binomial_distribution<int> distributionBl(L, 1);
			Bl = distributionBl(rd);
		}

		boost::binomial_distribution<int> distributionBp(P, (dP + uP)*dt);
		Bp = distributionBp(rd);

		boost::binomial_distribution<int> distributionBm(M, uM*dt);
		Bm = distributionBm(rd);

		boost::binomial_distribution<int> distributionNt(M, (dt / S));
		nt = distributionNt(rd);

		if (n*nt > 0) {
			boost::poisson_distribution<long unsigned int> distributionRp(n*nt);
			E = E - Be + distributionRp(rd);
		}
		else E = E - Be;

		boost::binomial_distribution<int> distributionL(Be, (dE / (uE + dE)));
		L = L - Bl + distributionL(rd);

		boost::binomial_distribution<int> distributionP(Bl, (dL / (uL + dL)));
		P = P - Bp + distributionP(rd);

		boost::binomial_distribution<int> distributionM(Bp, (dP / (uP + dP)));
		mRan = distributionM(rd);


		if (M + mRan - Bm > 0)
			M = M + mRan - Bm;
		else M = 0;

		t++;
		r.emplace_back(make_tuple(rint(E), rint(L), rint(P), rint(M)));
	}
	return r;
}

