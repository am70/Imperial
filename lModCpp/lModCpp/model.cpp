#include "lModH.h"

vector<tuple<int, int, int, int,double>> mPmod(modParms parmsx, boost::mt19937 rd) {

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
	double trx = parmsx.fxdPrm / parmsx.dt;
	double dt = parmsx.dt;
	double rFsum;
	double  uoE = parmsx.uoE;
	double  uoL = parmsx.uoL;
	double  Y = parmsx.Y;
	double  n = parmsx.n;
	double sf = parmsx.sf;
	double o = parmsx.o;
	double uP = parmsx.uP;
	double S = parmsx.S;
	double dE = parmsx.dE;// parmsx.dE;
	double dL = parmsx.dL;
	double uM = parmsx.uM;
	double uE;
	double uL;
	double uN;
	//double uoN =0.25;
	double rEff;

	double mRan;
	double dP = parmsx.dP;
	double Mg = parmsx.Mg;

	vector<double> rF = parmsx.rF;
	vector<tuple<int, int, int, int,double>> r;

	while (t < time) {

		if (t <= trx) {
			rFsum = std::accumulate(rF.begin(), rF.begin() + t, 0.0);
			K = ((sf*((1 / t)*rFsum)));
		}
		else {
			rFsum = std::accumulate(rF.begin() + (t - trx), rF.begin() + t, 0.0);
			K =  ((sf*((1 / trx)*rFsum)));
		}

	uE = uoE*exp(((E + L) / (K)));
		uL = uoL*exp(((Y*(E + L) / (K))));

		if (uL < 0)
			uL = 0;
		if (uE < 0)
			uE = 0;

		if ((dE + uE)*dt < 1) {
			boost::binomial_distribution<int> distributionBe(E, (dE + uE)*dt);
			Be = distributionBe(rd);
		}
		else {
			boost::binomial_distribution<int> distributionBe(E, 1);
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
			E = rint(E - Be + distributionRp(rd));
		}
		else E = rint(E - Be);

		boost::binomial_distribution<int> distributionL(Be, (dE / (uE + dE)));
		L = rint(L - Bl + distributionL(rd));

		boost::binomial_distribution<int> distributionP(Bl, (dL / (uL + dL)));
		P = rint(P - Bp + distributionP(rd));

		boost::binomial_distribution<int> distributionM(Bp, (dP / (uP + dP)));
		mRan = distributionM(rd);

		if (M + mRan - Bm > 1)
			M = rint(M + mRan - Bm);
		else M = 1;

		boost::poisson_distribution<long unsigned int> distributionRp2(rint(Mg+1));
		M = M + distributionRp2(rd);


		rEff = 0.5*((93.6*dt) / (exp(uM*S) - 1))*(1 / (1 + uE / dE))*(1 / (1 + uL / dL))*(1 / (1 + (uP*dt) / dP));

		t++;
		r.emplace_back(make_tuple(E, L, P, M,rEff));


	}

	//ofstream myfile2;
	//myfile2.open("Q:\\Imperial\\fitPlots\\sTest3.txt");
	//for (auto iter = 0; iter != size(r); ++iter) {
	//	cout << get<3>(r[iter]) << endl;
	//	myfile2 << get<3>(r[iter]) << ","<<endl;
	//}
	//cin.get();
	return r;
}

