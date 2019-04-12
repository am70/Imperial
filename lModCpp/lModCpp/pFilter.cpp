#include "lModH.h"
boost::mt19937 mrand(std::time(0));
double m_pi = 3.14159265358979323846;


/* initial state sampler, samples random initial states # need to fit initial E and use this to inform L, P & M
@param N number of particles for particle filter
@param t starting time period for rainfall
@param prms model parameters
@fxdParams fixed parameter (currently just for n in mosParamsP) - maybe update for multiple fixed parameters
@return conditions for E, L, P and M, double is empty for addition of weight later. */
vector<tuple<double, double, double, double, double>> iState(int N, int time, modParms iParms, string dFunc) {
	try {
		double z = iParms.z;
		double dE = iParms.dE;
		double dL = iParms.dL;
		double dP = iParms.dP;
		double y = iParms.Y;
		double S = iParms.S;
		double sf = iParms.sf;
		double UoE = iParms.uoE;
		double UoL = iParms.uoL;
		vector<double> rF = iParms.rF;
		double uM = iParms.uM;
		double uP = iParms.uP;
		double B = iParms.B;
		double K;
		double uE;
		double uL;
		double E;
		double L;
		double P;
		double M;
		double t = iParms.startTime;
		double a;
		double Mg = iParms.Mg;
		//int t = iParms.startTime;
		double trx = rint(iParms.tau / iParms.dt);
		double rFsum;
		double n = iParms.n;
		


		if (dFunc == "powerNoClumped" || dFunc == "expNoClumped" || dFunc == "linearNoClumped" || dFunc == "logisticNoClumped") n = iParms.n * (exp(S*uM) - 1) / uM;
		


	

		//exp carrying cap
		if (dFunc == "expNoClumped" || dFunc == "expClumped") {
			double calc;

	/*		for (int r = 0; r <= t; r++) {
				calc = (-(t - r)) / trx;
				rFsum = rFsum + exp(calc) * rF[r];
			}

			K = sf * (1.0 / (trx* (1 - exp(-(t + 1) / trx)))) * rFsum;*/

			if (t <= trx) {
				rFsum = std::accumulate(rF.begin(), rF.begin() + (t - 1), 0.0);
				K = ((sf*((1 / trx)*rFsum)));
			}
			else {
				rFsum = std::accumulate(rF.begin() + ((t - 1) - trx), rF.begin() + (t - 1), 0.0);//t-1 at begining as c++ starts on 0
				K = ((sf*((1 / trx)*rFsum)));
			}


			a = (0.5 * dL*dP) / (uM*(dP + uP));

			uE = UoE * exp((z / K));
			uL = UoL * exp((y*z / K));
			E = round(((B*a)*z) / ((dE + (B*a) + uE)));
			L = round((dE*z) / ((dE + dL + uL)));
			M = round((0.5 * dL*dP*L) / (uM*(dP + uP)));
			P = round((2 * uM*M) / dP);
		}
		//else if (dFunc == "linearNoClumped" || dFunc == "linearClumped") {
		//	//linear carrying cap
		//	if (dFunc == "linearClumped") n = n  / (exp(S*uM) - 1)*uM;
		//	dE = 1 / dE;
		//	dL = 1 / dL;
		//	dP = 1 / dP;
		//	M = z;
		//	double W = -0.5*(y*(UoL / UoE) - (dE / dL) + (y - 1)*UoL*dE) + sqrt(0.25*pow((y*(UoL / UoE) - (dE / dL) + (y - 1)*UoL*dE), 2) + y * ((n*UoL*dE) / (2 * UoE*uM*dL*(1 + dP * uP))));
		//	E = rint(2 * W*uM*dL*(1 + dP * uP)*M);
		//	L = rint(2 * uM*dL*(1 + dP * uP)*M);
		//	P = rint(2 * dP*uM*M);
		//}
		else {
			//power carrying cap
			a = ((n / S) * dP * dL) / ((2 * uM) * (uP + dP));
			double b = (UoE / (y * UoL)) * (dL + UoL) - dE - UoE;
			double c = -(UoE * dE) / (UoL * y);
			double x = (-b + sqrt(pow(b, 2) - 4 * a * c)) / (2 * a);

			L = z;
			E = round(L / x);
			P = round((dL * L) / (uP + dP));
			M = round((dP * P) / (2 * uM));
		}



		if (Mg > 1) {
			boost::poisson_distribution<long unsigned int> distributionRp2(round(Mg));
			M = M + distributionRp2(mrand);
		}
		if (M < 1)
			M = 1;
		vector<tuple<double, double, double, double, double>> states = { { E, L, P, M ,0.0 } };
		states.resize(N, { E, L, P, M, 0.0 });//return initial states, repeated to number of particles
		return states;
	}
	catch (...) { cerr << "error in initial states function, check input" << endl;  cin.get();
	}
}


/*binomial function
@param n number of trials
@param p probability of success
return number of successes*/
int binom(int n, double p, boost::mt19937 rd) {
	boost::binomial_distribution<int> distribution(n, p);
	int res = distribution(rd);
	return  res;
}


/*random poisson draw
@param lambda
@return random poisson draw*/
luint rpois(luint lambda, boost::mt19937 rd) {
	if (lambda > 0) {
		boost::poisson_distribution<luint> d(lambda);
		return d(rd);
	}
	else return (0);
}

/*lbeta function
@param a
@param b
@return log beta*/
double lbeta(double a, double b) {
	return double(lgamma(a) + lgamma(b)) - lgamma(a + b);
}

//double lchoose(double n, double k) {
//	double 	nk = log(n - k);
//	double ln = log(n);
//	double lk = log(k);
//	double res = boost::math::tgamma<double>(1+ln) - boost::math::tgamma<double>(1+lk)- boost::math::tgamma<double>(1 + nk);
//	cout << res <<endl;
//	return res;
//}

/*negBin function
@param k Observed data point
@param n Simulated data point
@param p fraction of population*/
double nBgP(double k, double n, double p) {
	if (p > 1)
		p = 1;
	if(k==0)
		k=1;
		double res =  (lgamma(n + k) - (lgamma(k) + lgamma(1 + n))) + (k*log(p)) + n*log(1-p);
		return res;
}



/*negBin function
@param k Observed data point
@param n Simulated data point
@param r overdispersion parameter
@param p fraction of population*/
double nB(double k, double n, double r, double p) {
	double m = p*n;
	double res = (lgamma(r + k) - (lgamma(k + 1) + lgamma(r))) + k*(log(m) - (log(r + m))) + r*(log(r) - log(r + m));
	return res;
}

/*binomial function
@param k Observed data point
@param n Simulated data point
@param p fraction of population*/
double dbinom(double k, double n, double p) {
	if (n >= k) {
		if (p > 1)
			p = 1;
		if (k == 0)
			k = 1;
		double res = (boost::math::lgamma(n + 1) - (boost::math::lgamma(k + 1) + boost::math::lgamma(n - k + 1))) + (k*log(p)) + ((n - k)*log(1 - p));
		return res;
	}
	else return -1000;

}

/* Beta binomial likelihood function
@param k Observed data point
@param n Simulated data point
@param p fraction of population
@param w overdisperion parameter
@return log likelihood*/
double betaBinom(double k, double n, double p, double w) {
	try {
		n = n;
		if (n >= k) {
			if (n <= 0) return -1000;
			else {
				if (k <= 0) k = 1;
				if (k == n) n = n + 1;
				double a = p * ((1 / w) - 1);
				double b = (1 - p) * ((1 / w) - 1);
				double logn = log(n);
				double logk = log(k);
				double lognk = log(n - k);
				double res = lbeta(k + a, n - k + b) - lbeta(a, b) + boost::math::lgamma(n + 1) - (boost::math::lgamma(k + 1) + boost::math::lgamma(n - k + 1));
				return res;
			}
		}
		else return -1000;
	}
	catch (...) { cerr << "betaBinom Func error, check input, currently:" << endl << "k = "<<k
		<<endl << "n = " << n << endl << "p = " << p << endl << "w = " << w << endl; cin.get(); }

}


/*model step function
@param wp a structure containing start times and states for model step
@param obsData, observed data for calculating likelihood value
@return tuple containing state ints for E, L, P, M and likelihood (non-log), at the end of 
a model run with designated start and end times*/
tuple<double, double, double, double, double> modStepFnc(modParms wp, int obsData, boost::mt19937_64 rd, string dFunc) {
	vector<tuple<double, double, double, double,double>> modRun;
	modRun.reserve(wp.endTime - wp.startTime);
	double weight; 
	try {
		modRun = mPmod(wp, rd, dFunc);
	}
	catch (...) {
		cerr << "error with model run in modStep func" << endl; cin.get();
	}
	double sim = get<3>(modRun.back());
	weight = betaBinom(obsData, sim, wp.p,wp.w); //add weights to tuple
	tuple<double, double, double, double, double> res = { get<0>(modRun.back()),get<1>(modRun.back()),get<2>(modRun.back()),get<3>(modRun.back()), weight};
	return res;
}

/*model step function for plotting
@param wp a structure containing start times and states for model step
@param obsData, observed data for calculating likelihood value
@return tuple containing state ints for E, L, P, M and likelihood, at the end of
a model run with designated start and end times*/
vector<tuple<double, double, double, double, double,double>> modStepFncPlot(modParms wp, int obsData, boost::mt19937_64 rd, string dFunc) {
	vector<tuple<double, double, double, double,double>> modRun;
	modRun.reserve(wp.endTime - wp.startTime);
	double weight;
	
	try {
		modRun = mPmod(wp, rd, dFunc);
	}
	catch (...) {
		cerr << "error with model run in modStepFncPlot func" << endl; cin.get();
	}	
	
	double sim = get<3>(modRun.back());
	weight = betaBinom(obsData, sim,wp.p,wp.w);//add weights to tuple
	tuple<double, double, double, double, double,double> res;
	vector<tuple<double, double, double, double, double,double>> res2;

	for (int j = 0; j < boost::size(modRun); j++) {
		res = { get<0>(modRun[j]),get<1>(modRun[j]),get<2>(modRun[j]),get<3>(modRun[j]), weight,get<4>(modRun[j])};
		res2.emplace_back(res);
	}

	return res2;
}


/*rand particle sample function
@param samp vector of tuples containing particles E,L,M,P
@return resamped particles*/
vector<tuple<double, double, double, double, double>> rSamp(vector<std::tuple<double, double, double, double, double>>& samp) {
	vector<std::tuple<double, double, double, double, double>> temp;
	vector<double> w;// weights
	for (auto&& i : samp)
		w.push_back(get<4>(i)); // get weights
	std::discrete_distribution<int> dd{ w.begin(), w.end() }; // create distribution
	for (size_t i{}; i < samp.size(); ++i)
		temp.push_back(samp[dd(mrand)]); // build return vector by selecting according to weights
	return temp;
}


/*rand result matrix sample function
@param samp vector of vectors containing trajectories of model runs
@return id number for random model run based on final likelihood
vec<particles<vec<times>>*/
double rSampMat(vector<vector<double>> vP) {
	vector<tuple<double, double>>temp;
	tuple<double, double>tempRes;
	vector<double> w;
	for (int i = 0; i != boost::size(vP[0]); i++) { //for i in number of particles
		temp.emplace_back(vP.at(boost::size(vP[0]))[i],i);//tuple 0 = end of particle vec likelihood
		w.emplace_back(exp(vP.at(boost::size(vP[0]))[i]));
	}
		std::discrete_distribution<int> dd{ w.begin(), w.end()}; // create distribution
		tempRes=(temp[dd(mrand)]);
	return accumulate(vP[get<1>(tempRes)].begin(), vP[get<1>(tempRes)].end(), 0.0);
}

/*normalise weights
@param particles vector of tuples for particles
@param llsum current sum of log likelihoods in particles*/
vector<tuple<double, double, double, double, double>> normalise(vector<tuple<double, double, double, double, double>> particles, double llSum) {
	const auto pComp = [](const auto& lhs, const auto& rhs)
	{ return get<4>(lhs) < get<4>(rhs);};
	const double maxVal = get<4>(*std::max_element(particles.begin(), particles.end(), pComp));
	double sum = 0;
	for (int i = 0; i != boost::size(particles); i++) {

		get<4>(particles[i]) = (get<4>(particles[i])) - maxVal;
		sum = sum + exp((get<4>(particles[i])));

	}

	for (int i = 0; i != boost::size(particles); i++) {
		get<4>(particles[i]) = exp(get<4>(particles[i])) / sum;
	}

	return particles;
}




/* Particle filter function - runs model in steps between each observed data point
before using weighted re - sampling of model state at each step to start next
@param n number of particles 
@param iState function for calculating initial state of particles
@param obsData observed data to fit model to
@param prms model parameters
@param resM True/False whether to output likelihood values or results of simulation for plotting figures - not currently implemented
@param fxedParams fixed parameters - currently just accepting tr but may expand to include more than one parameter
@return log likelihood value*/
double pFilt(int n,
	vector<tuple<int, int>> obsData,
	modParms prms,
	bool resM,
	string outputFile,
	bool reff,
	string dFunc
	)
{

		vector<int> times;
		double ll; //log likelihood value
		modParms wp;
		vector<tuple<double, double, double, double, double>>  particles;
		particles.reserve(n);
		vector<double> modPlotRes;
		vector<double> modPlot;
		vector<double> modPlotReff;
		vector<double> plotResTemp;
		vector<double> plotResTempReff;


		//calculate discreet time steps for observed data
		for (auto i = begin(obsData); i != end(obsData); ++i) {
			times.emplace_back((get<0>(*i)) / prms.dt);
		};
		prms.startTime = times.at(0);
		//get initial state of particles
		particles = iState(n, times.front(), prms, dFunc);

		ll = (betaBinom(get<1>(obsData[0]), get<3>(particles[0]), prms.p, wp.w));

		vector<vector<double>> resMat(size(particles), vector<double>(size(obsData)));

		for (auto i = 0; i != times.size() - 1; ++i) {
			modParms wp;
			wp.dE = prms.dE;
			wp.dL = prms.dL;
			wp.dP = prms.dP;

			wp.uoE = prms.uoE;
			wp.uoL = prms.uoL;
			wp.uP = prms.uP;
			wp.uM = prms.uM;

			wp.Y = prms.Y;
			wp.sf = prms.sf;
			wp.z = prms.z;
			wp.n = prms.n;
			wp.w = prms.w;
			wp.rF = prms.rF;
			wp.Mg = prms.Mg;
			wp.p = prms.p;
			wp.lK = prms.lK;
			wp.lKs = prms.lKs;
			wp.lKm = prms.lKm;


			if (dFunc == "expNoClumped" || dFunc == "linearNoClumped" || dFunc == "linearClumped" || dFunc == "expClumped")
				wp.o = 1;
			else
				wp.o = prms.o;

			//if (dFunc != "logisticClumped" || "logisticNoClumped")
				//wp.lK == 0;

			wp.tau = prms.tau;

			wp.startTime = times.at(i);
			wp.endTime = times.at(i + 1);
			//run model in step and update particles in parallel
			double lltemp = 0;


			if (resM == false) {


#pragma omp parallel for schedule(static) reduction(+:lltemp) //reduction needed due to problem with thread racing
				for (int j = 0; j < boost::size(particles); j++) {

					std::random_device rdev;
					uint64_t seed = (uint64_t(rdev()) << 32) | rdev();

					boost::mt19937_64 mrandThread(seed);
					wp.E0 = get<0>(particles[j]);
					wp.L0 = get<1>(particles[j]);
					wp.P0 = get<2>(particles[j]);
					wp.M0 = get<3>(particles[j]);

					int obsDatPoint = get<1>(obsData[i]);
					particles.at(j) = modStepFnc(wp, obsDatPoint, mrandThread, dFunc);

					resMat[j][i] = get<4>(particles.at(j));

					lltemp = lltemp + get<4>(particles.at(j));
					//if (std::isnan(get<4>(particles.at(j))) == 0)  get<4>(particles.at(j)) = 0;
				}

				double llMean = lltemp / boost::size(particles);
				ll = ll + llMean;//add start val for frst obs point
				//normalise particle probabilities
				particles = normalise(particles, lltemp);
				//re-sample particles
				particles = rSamp(particles);

			}

			//alternative loop for getting plot outputs from particle filter - could be coded better...
			if (resM == true) {
				vector < tuple<double, double, double, double, double, double>> plotRes;
				plotResTemp.clear();
				plotResTemp.resize(wp.endTime - wp.startTime);
				plotResTempReff.clear();
				plotResTempReff.resize(wp.endTime - wp.startTime);
				for (int j = 0; j < boost::size(particles); j++) {
					std::random_device rdev;
					uint64_t seed = (uint64_t(rdev()) << 32) | rdev();
					boost::mt19937_64 mrandThread(seed);
					wp.E0 = get<0>(particles[j]);
					wp.L0 = get<1>(particles[j]);
					wp.P0 = get<2>(particles[j]);
					wp.M0 = get<3>(particles[j]);

					int obsDatPoint = get<1>(obsData[i]);
					plotRes = modStepFncPlot(wp, obsDatPoint, mrandThread, dFunc);
					for (int x = 0; x < boost::size(plotRes); x++) {
						plotResTempReff.at(x) = plotResTempReff.at(x) + (get<5>(plotRes[x]));
						plotResTemp.at(x) = plotResTemp.at(x) + (get<3>(plotRes[x]));
					}
					particles.at(j) = { get<0>(plotRes.back()),get<1>(plotRes.back()),get<2>(plotRes.back()),get<3>(plotRes.back()), get<4>(plotRes.back()) };
					lltemp = lltemp + get<4>(particles.at(j));
				}
				for (int g = 0; g < boost::size(plotResTemp); g++) {
					modPlot.emplace_back(plotResTemp.at(g) / boost::size(particles));
					modPlotReff.emplace_back(plotResTempReff.at(g) / boost::size(particles));
				}
				//take mean of likelihoods from particles
				double llMean = lltemp / boost::size(particles);
				//normalise particle probabilities
				particles = normalise(particles, lltemp);
				//re-sample particles
				particles = rSamp(particles);
				ll = ll + llMean;

			}


		}
		if (resM == true) {
			string opF = outputFile;
			outputFile.append(".txt");
			ofstream myfile;
			myfile.open(outputFile);
			for (int j = 0; j < boost::size(modPlot); j++) {
				myfile << modPlot.at(j) << endl;
			}

			opF.append("Reff.txt");
			ofstream myfileReff;
			myfileReff.open(opF);
			for (int j = 0; j < boost::size(modPlotReff); j++) {
				myfileReff << modPlotReff.at(j) << endl;
			}
		}

		//take random sample at end of resMat
		double resNum = ll;// rSampMat(resMat);

		return resNum;
	
}