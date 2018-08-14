#include "lModH.h"

vector<tuple<int, int, int, int,double>> mPmod(modParms parmsx, boost::mt19937_64 rd, string dFunc) {

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
	double K;
	double trx = round(parmsx.tau / parmsx.dt);
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
	vector<double> sortrF;
	double mRan;
	double dP = parmsx.dP;
	double Mg = parmsx.Mg;
	double uEn;
	double uLn;
	double lK = parmsx.lK;
	vector<double> rF = parmsx.rF;
	vector<tuple<int, int, int, int,double>> r;


	while (t < time) {

		if (t <= trx) {
			rFsum = std::accumulate(rF.begin(), rF.begin() + (t-1) , 0.0);
			K = ((sf*((1 / trx)*rFsum)));
		}
		else {
			rFsum = std::accumulate(rF.begin() + ((t-1) - trx), rF.begin() + (t-1), 0.0);//t-1 at begining as c++ starts on 0
			K =  ((sf*((1 / trx)*rFsum)));
		}

		// K = K+-lK*pow(K,2);


		//((sf*(1 / (trx*(1 - exp(-time / trx))))*rFsum)); exp rainfall carrying cap

	
	
		if (dFunc == "powerNoClumped" || dFunc == "linearNoClumped" || dFunc == "powerClumped" || dFunc == "linearClumped" || dFunc == "logisticClumped"|| dFunc == "logisticNoClumped") {
			uE = uoE * (1 + pow(((E + L) / (K)), o));
			uL = uoL * (1 + (Y*pow(((E + L) / (K)), o)));
			 uEn = uoE * (1 + pow(((n) / (K)), o));
			 uLn = uoL * (1 + (Y*pow(((n) / (K)), o)));
		}
		else {
			uE = uoE * exp(((E + L) / (K)));
			uL = uoL * exp(Y*(((E + L) / (K))));
			 uEn = uoE * exp(((n) / (K)));
			 uLn = uoL * exp(Y*(((n) / (K))));
		}



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
				boost::binomial_distribution<int> distributionBl(L,1);
				Bl = distributionBl(rd);
			}


			if ((dP + uP)*dt < 1) {
				boost::binomial_distribution<int> distributionBp(P, (dP + uP)*dt);
				Bp = distributionBp(rd);
			}
			else {
				boost::binomial_distribution<int> distributionBp(P, 1);
				Bp = distributionBp(rd);
			}
		
	
			if (M >= 1) {
				boost::binomial_distribution<int> distributionBm(M, uM*dt);
				Bm = distributionBm(rd);
			}
			else Bm = 0;


			if (dFunc == "expClumped"|| dFunc == "linearClumped"|| dFunc == "powerClumped" || dFunc == "logisticClumped") {
				if (M >= 1) {
					boost::binomial_distribution<int> distributionNt(M, (dt / S));
					nt = distributionNt(rd);
				
				}
				else nt = 0;

				if (n*nt > 0) {
					boost::poisson_distribution<long unsigned int> distributionRp(n*nt);
					E = round(E - Be + distributionRp(rd));
				}
				else E = round(E - Be);
			}
			else E = round(E - Be + M * (n*dt));


		if (Be >= 1) {
			boost::binomial_distribution<int> distributionL(Be, (dE / (uE + dE)));
			L = round(L - Bl + distributionL(rd));
		}
		else 
			L = round(L - Bl);

		if (Bl >= 1) {
			boost::binomial_distribution<int> distributionP(Bl, (dL / (uL + dL)));
			P = round(P - Bp + distributionP(rd));
		}
		else 
			P = round(P - Bp);

		if (Bp >= 1) {
			boost::binomial_distribution<int> distributionM(Bp, (dP / (uP + dP)));
			mRan = distributionM(rd);
		}
		else
			mRan = 0; 

		if (Mg > 0) {
			boost::poisson_distribution<long unsigned int> distributionRp2(Mg);
			M = round(M + (0.5*mRan) - Bm + distributionRp2(rd));
		}
		
		else 
			M = round(M + (0.5*mRan) - Bm);

		if (M < 1)
			M = 1;

		if (dFunc == "expClumped" || dFunc == "linearClumped" || dFunc == "powerClumped" || dFunc == "logisticClumped") {
			rEff = 0.5*((n) / (exp(uM*S) - 1))*(1 / (1 + uE / dE))*(1 / (1 + uL / dL))*(1 / (1 + (uP) / dP));

		} else { 
			double Es = n * (exp(S*uM) - 1) / uM;//calculate value for oviposition cycle from daily rate
			rEff = 0.5*((Es) / (exp(uM*S) - 1))*(1 / (1 + uE / dE))*(1 / (1 + uL / dL))*(1 / (1 + (uP) / dP));}

		t++;
		r.emplace_back(make_tuple(E, L, P, M, rEff));

	}

	//ofstream myfile2;
	//myfile2.open("Q:\\Imperial\\fitPlots\\sTest3.txt");
	//for (int iter = 0; iter != size(r); ++iter) {
	//	myfile2 << get<4>(r[iter]) <<endl;
	//}
	//cin.get();
	return r;
}

