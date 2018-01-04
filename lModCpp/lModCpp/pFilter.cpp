#include "lModH.h"
boost::mt19937 mrand(std::time(0));
double inf = std::numeric_limits<double>::infinity();
double m_pi = 3.14159265358979323846;


/* initial state sampler, samples random initial states # need to fit initial E and use this to inform L, P & M
@param N number of particles for particle filter
@param t starting time period for rainfall
@param prms model parameters
@fxdParams fixed parameter (currently just for n in mosParamsP) - maybe update for multiple fixed parameters
@return conditions for E, L, P and M, double is empty for addition of weight later. */
vector<tuple<int, int, int, int, double>> iState(int N, int time, modParms iParms, int fxdParm) {

	int z = iParms.z;//fitted(E - L)
	double dE = iParms.dE;
	double dL = iParms.dL;
	double dP = iParms.dP;
	double y = iParms.Y;
	double S = iParms.S;
	double sf = iParms.sf;
	double n = iParms.n;
	double UoE = iParms.uoE;
	double UoL = iParms.uoL;
	vector<double> rF = iParms.rF;
	double uM = iParms.uM;
	double uP = iParms.uP;
	double B = iParms.B;
	double K;
	double uE;
	double uL;
	int E;
	double L;
	int P;
	int M;
	double t = iParms.startTime;
	double a;
	//int t = iParms.startTime;
	double trx = iParms.tr / iParms.dt;
	double rFsum;

	if (t<= trx) {
		rFsum = std::accumulate(rF.begin(), rF.begin() + t, 0);
		K = ((sf*((1 / t)*rFsum)));
	}
	else {
		rFsum = std::accumulate(rF.begin() + (t - trx), rF.begin() + t, 0);
		K = ((sf*((1 / trx)*rFsum)));
	}

	a = (0.5 * dL*dP) / (uM*(dP + uP));

	uE = UoE*exp((z / K));
	uL = UoL*exp((y*z / K));

	E = rint(((B*a)*z) / ((dE + (B*a) + uE)));
	L = rint((dE*z) / ((dE + dL + uL)));
	M = rint((0.5 * dL*dP*L) / (uM*(dP + uP)));
	P = rint((2 * uM*M) / dP);

	/*a = ((n / S) * dP * dL) / ((2 * uM) * (uP * dP));
	double	b = (UoE / (y * UoL)) * (dL + UoL) - dE - UoE;
	double	c = -(UoE * dE) / (UoL * y);
	double	x = (-b + sqrt(pow(b,2) * -4 * a * c)) / (2 * a);

	L = z;
	E =rint(L / x);
	P =rint((dL * L) / (uP + dP));
	M =rint((dP * P) / (2 * uM));
*/

	if (M < 20) M = 20;

	vector<tuple<int, int, int, int, double>> states = { { E, L, P, M ,0.0 } };
	states.resize(N, { E, L, P, M, 0.0 });//return initial states, repeated to number of particles
	return states;
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

//negBin function
double nBgP(double k, double n, double p) {
	if (n >= k) {
		double res = (lgamma(n + k) - (lgamma(k) + lgamma(1 + n))) + (log(pow(p, k)) + log(pow((1 - p), n)));
	//	cout << " k = " << k << " n =" << n << " res = " << res << endl;
		return res;
	}
	else return -10000;
}

/* Beta binomial likelihood function
@param k Observed data point
@param n Simulated data point
@param p fraction of population
@param w overdisperion parameter
@return log likelihood*/
double betaBinom(double k, double n, double p, double w) {
	//n = n + 100;
	if (n >= k) {
		if (k <= 0) k = 1;
		if (k == n) n = n + 1;
		double a = p * ((1 / w) - 1);
		double b = (1 - p) * ((1 / w) - 1);
		double logn = log(n);
		double logk = log(k);
		double lognk = log(n - k);
		double res = lbeta(k + a, n - k + b) - lbeta(a, b) + (n *logn - k*logk - (n - k)* lognk + 0.5*(logn - logk - lognk - log(2 * m_pi)));//using stirling approximation for bionomial coefficient

		return res;
	}

	else
		//cout << " k = " << k << " n =" << n << endl;

		return -1000;
}


/*model step function
@param wp a structure containing start times and states for model step
@param obsData, observed data for calculating likelihood value
@return tuple containing state ints for E, L, P, M and likelihood (non-log), at the end of 
a model run with designated start and end times*/
tuple<int, int, int, int, double> modStepFnc(modParms wp, int obsData, boost::mt19937 rd) {
	vector<tuple<int, int, int, int>> modRun;
	modRun.reserve(wp.endTime - wp.startTime);
	double weight; 
	try
	{
		modRun = mPmod(wp, rd);
	}
	catch (...) { cout << "mod run exception"; }
	double sim = get<3>(modRun.back());
	weight = betaBinom(obsData, sim, 0.01,wp.w); //add weights to tuple
	tuple<int, int, int, int, double> res = { get<0>(modRun.back()),get<1>(modRun.back()),get<2>(modRun.back()),get<3>(modRun.back()), weight };
	return res;
}

/*model step function for plotting
@param wp a structure containing start times and states for model step
@param obsData, observed data for calculating likelihood value
@return tuple containing state ints for E, L, P, M and likelihood (non-log), at the end of
a model run with designated start and end times*/
vector<tuple<int, int, int, int, double>> modStepFncPlot(modParms wp, int obsData, boost::mt19937 rd) {
	vector<tuple<int, int, int, int>> modRun;
	modRun.reserve(wp.endTime - wp.startTime);
	double weight;
	try
	{

		modRun = mPmod(wp, rd);
	}
	catch (...) { cout << "mod run exception"; }

	double sim = get<3>(modRun.back());
	weight = betaBinom(obsData, sim, 0.01, wp.w); //add weights to tuple
	tuple<int, int, int, int, double> res;
	vector<tuple<int, int, int, int, double>> res2;
	for (int j = 0; j < boost::size(modRun); j++) {
		res = { get<0>(modRun[j]),get<1>(modRun[j]),get<2>(modRun[j]),get<3>(modRun[j]), weight };
		res2.emplace_back(res);
	}
	return res2;
}


/*rand particle sample function
@param samp vector of tuples containing particles E,L,M,P
@return resamped particles*/
auto rSamp(vector<std::tuple<int, int, int, int, double>>& samp) {
	vector<std::tuple<int, int, int, int, double>> temp;
	vector<double> w;// weights
	for (auto&& i : samp)
		w.push_back(get<4>(i)); // get weights
	std::discrete_distribution<int> dd{ w.begin(), w.end() }; // create distribution
	for (size_t i{}; i < samp.size(); ++i)
		temp.push_back(samp[dd(mrand)]); // build return vector by selecting according to weights
	return temp;
}

/*normalise weights
@param particles vector of tuples for particles
@param llsum current sum of log likelihoods in particles*/
vector<tuple<int, int, int, int, double>> normalise(vector<tuple<int, int, int, int, double>> particles, double llSum) {
	const auto pComp = [](const auto& lhs, const auto& rhs)
	{ return get<4>(lhs) < get<4>(rhs); };
	const double maxVal = get<4>(*std::max_element(particles.begin(), particles.end(), pComp));

	int sum = 0;
	for (int i = 0; i != boost::size(particles); i++) {
		get<4>(particles[i]) = (get<4>(particles[i])) - maxVal;
		sum = sum + exp((get<4>(particles[i])) - maxVal);
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
	vector< tuple<int, int> > obsData,
	modParms prms,
	bool resM,
	int fxdParams,
	string outputFile) {
	vector<int> times;
	double ll; //log likelihood value
	modParms wp;
	vector<tuple<int, int, int, int, double>>  particles;
	particles.reserve(n);
	vector<double> modPlotRes;
	vector<double> modPlot;
	vector<double> plotResTemp;


	//calculate discreet time steps for observed data
	for (auto i = begin(obsData); i != end(obsData); ++i) {
		times.emplace_back((get<0>(*i)) / prms.dt);
	};
	prms.startTime = times.at(0);
	//get initial state of particles
	particles = iState(n, times.front(), prms, fxdParams);

		ll = (betaBinom(get<1>(obsData[0]), get<3>(particles[0]), 0.01,prms.w));
		

		for (auto i = 0; i != times.size() - 1; ++i) {
			modParms wp;
			wp.dE = prms.dE;
			wp.dL = prms.dL;
			wp.dP = prms.dP;

			wp.uoE = prms.uoE;
			wp.uoL = prms.uoL;
			wp.uP = prms.uP;
			wp.Y = prms.Y;
			wp.sf = prms.sf;
			wp.z = prms.z;
			wp.n = prms.n;
			wp.w = prms.w;
			wp.rF = prms.rF;
			wp.fxdPrm = fxdParams;
			wp.startTime = times.at(i);
			wp.endTime = times.at(i + 1);
			//run model in step and update particles in parallel
			double lltemp = 0;


			if (resM == false) {
#pragma omp parallel for schedule(static) reduction(+:lltemp)  //reduction needed due to problem with thread racing
				for (int j = 0; j < boost::size(particles); j++) {
					boost::mt19937 mrandThread(std::random_device{}());
					wp.E0 = get<0>(particles[j]);
					wp.L0 = get<1>(particles[j]);
					wp.P0 = get<2>(particles[j]);
					wp.M0 = get<3>(particles[j]);

					int obsDatPoint = get<1>(obsData[i]);
					particles.at(j) = modStepFnc(wp, obsDatPoint, mrandThread);
					lltemp = lltemp + get<4>(particles.at(j));
					//if (std::isnan(get<4>(particles.at(j))) == 0)  get<4>(particles.at(j)) = 0;
				}
				//normalise particle probabilities
				particles = normalise(particles, lltemp);
				//re-sample particles
				particles = rSamp(particles);
				//take mean of likelihoods from particles
				double llMean = lltemp / boost::size(particles);
				ll = ll + llMean;//add start val for frst obs point
			}

			//alternative loop for getting plot outputs from particle filter - could be coded better...
			if (resM == true) {
				vector < tuple<int, int, int, int, double>> plotRes;
				plotResTemp.clear();
				plotResTemp.resize(wp.endTime - wp.startTime);
				for (int j = 0; j < boost::size(particles); j++) {
					boost::mt19937 mrandThread(std::random_device{}());
					wp.E0 = get<0>(particles[j]);
					wp.L0 = get<1>(particles[j]);
					wp.P0 = get<2>(particles[j]);
					wp.M0 = get<3>(particles[j]);

					int obsDatPoint = get<1>(obsData[i]);
					plotRes = modStepFncPlot(wp, obsDatPoint, mrandThread);
					for (int x = 0; x < boost::size(plotRes); x++) {
						plotResTemp.at(x)= plotResTemp.at(x)+(get<3>(plotRes[x]));
					}
					particles.at(j) = { get<0>(plotRes.back()),get<1>(plotRes.back()),get<2>(plotRes.back()),get<3>(plotRes.back()), get<4>(plotRes.back()) };;
					lltemp = lltemp + get<4>(particles.at(j));
					//if (std::isnan(get<4>(particles.at(j))) == 0)  get<4>(particles.at(j)) = 0;
				}
				for (int g = 0; g < boost::size(plotResTemp); g++) {
					modPlot.emplace_back(plotResTemp.at(g) / boost::size(particles));
				}
				//normalise particle probabilities
				particles = normalise(particles, lltemp);
				//re-sample particles
				particles = rSamp(particles);
				//take mean of likelihoods from particles
				double llMean = lltemp / boost::size(particles);
				ll = ll + llMean;//add start val for frst obs point

			}

	
	}
		if (resM == true) {

			for (int j = 0; j < boost::size(modPlot); j++) {
				cout <<j << ",";
			}

			ofstream myfile;
			myfile.open(outputFile);
			for (int j = 0; j < boost::size(modPlot); j++) {
				myfile << modPlot.at(j) <<endl;
			}
		}
	return ll;
}
