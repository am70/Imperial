#include "lModH.h"
boost::mt19937 mrand(std::time(0));
double inf = std::numeric_limits<double>::infinity();


// initial state sampler, samples random initial states # need to fit initial E and use this to inform L, P & M
// @param N number of particles for particle filter
// @param t starting time period for rainfall
// @param prms model parameters
// @fxdParams fixed parameter (currently just for n in mosParamsP) - maybe update for multiple fixed parameters
// @return conditions for E, L, P and M, double is empty for addition of weight later. 
vector<tuple<int, int, int, int, double>> iState(int N, int time, modParms iParms, int fxdParm) {

	double z = iParms.z;//fitted(E - L)
	double dE = iParms.dE;
	double dL = iParms.dL;
	double dP = iParms.dP;
	double y = iParms.Y;
	double S = iParms.S;
	double sf = iParms.sf;
	double n = iParms.n;
	double UoE = iParms.uoE;
	double UoL = iParms.uoL;
	vector<int> rF = iParms.rF;
	double uM = iParms.uM;
	double uP = iParms.uP;
	double tr = fxdParm / iParms.dt;
	double B = iParms.B;
	double K;
	double uE;
	double uL;
	double E;
	double L;
	double P;
	double M;
	double a;
	//int t = iParms.startTime;
	double trx = iParms.tr / iParms.dt;

	int rFsum = std::accumulate(rF.begin() + (time - trx), rF.begin() + time-1, 0);
	K = (1 + (sf*((1 / trx)*rFsum)));
	a = (0.5 * dL*dP) / (uM*(dP + uP));

	uE = UoE*exp(z / K);
	uL = UoL*exp(y*z / K);

	E = ((B*a)*z) / (dE + (B*a) + uE);
	L = (dE*z) / (dE + dL + uL);
	M = (0.5 * dL*dP*L) / uM*(dP + uP);
	P = 2 * uM*M / dP;

	vector<tuple<int, int, int, int, double>> states = { { rint(E), rint(L), rint(P), rint(M) ,0.0} };
	states.resize(N, { E, L, P, M, 0.0 });//return initial states, repeated to number of particles
	return states;
}


//binomial function
//@param n number of trials
//@param p probability of success
//@return number of successes
int binom(int n, double p, boost::mt19937 rd) {
	boost::binomial_distribution<int> distribution(n, p);
	int res = distribution(rd);
	return  res;
}


//random poisson draw
//@param lambda
//@return random poisson draw
luint rpois(luint lambda, boost::mt19937 rd){
	if (lambda > 0) {
		boost::poisson_distribution<luint> d(lambda);
		return d(rd);
	}
	else return (0);
}

//lbeta function
//@param a
//@param b
//@return log beta
double lbeta(double a, double b) {
	return double(lgamma(a) + lgamma(b)) - lgamma(a + b);
}

// Beta binomial likelihood function
// @param k Observed data point
// @param n Simulated data point
// @param p fraction of population
// @param w overdisperion parameter
// @return log likelihood
double betaBinom(double k, double n, double p, double w) {
	 if (n > k) {
		double a = p * ((1 / w) - 1);
		double b = (1 - p) * ((1 / w) - 1);
		return lbeta(k + a, n - k + b) - lbeta(a, b) + log(boost::math::binomial_coefficient<double>(n, k));
	} else
		 return -inf;
}


//model step function
//@param wp a structure containing start times and states for model step
//@param obsData, observed data for calculating likelihood value
//@return tuple containing state ints for E, L, P, M and likelihood (non-log), at the end of 
//a model run with designated start and end times
tuple<int, int, int, int, double> modStepFnc(modParms wp, int obsData, boost::mt19937 rd) {
	vector<tuple<int, int, int, int>> modRun;
	modRun.reserve(wp.endTime - wp.startTime);
	double weight;
	modRun = mPmod(wp, rd);
	double sim = 10 + get<3>(modRun.back());
	weight = betaBinom(obsData, sim, 0.01, wp.w); //add weights to tuple
	tuple<int, int, int, int, double> res = { get<0>(modRun.back()),get<1>(modRun.back()),get<2>(modRun.back()),get<3>(modRun.back()), weight };
	return res;
}

//rand particle sample function
//@param samp vector of tuples containing particles E,L,M,P
//@return resamped particles
auto rSamp(vector<std::tuple<int, int, int, int, double>>& samp){
	vector<std::tuple<int, int, int, int, double>> temp;
	vector<double> w;// weights
	for (auto&& i : samp)
		w.push_back(get<4>(i)); // get weights
	std::discrete_distribution<int> dd{ w.begin(), w.end() }; // create distribution
	for (size_t i{}; i < samp.size(); ++i)
		temp.push_back(samp[dd(mrand)]); // build return vector by selecting according to weights
	return temp;
}

//normalise weights
//@param particles vector of tuples for particles
//@param llsum current sum of log likelihoods in particles
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

// Particle filter function - runs model in steps between each observed data point
// before using weighted re - sampling of model state at each step to start next
// @param n number of particles 
// @param iState function for calculating initial state of particles
// @param obsData observed data to fit model to
// @param prms model parameters
// @param resM True/False whether to output likelihood values or results of simulation for plotting figures - not currently implemented
// @param fxedParams fixed parameters - currently just accepting tr but may expand to include more than one parameter
// @return log likelihood value
double pFilt(int n,
	vector< tuple<int, int> > obsData,
	modParms prms,
	bool resM,
	int fxdParams) {
	vector<int> times;
	double ll; //log likelihood value
	modParms wp;
	vector<tuple<int, int, int, int, double>>  particles;
	particles.reserve(n);

	//calculate discreet time steps for observed data
	for (auto i = begin(obsData); i != end(obsData); ++i) {
		times.emplace_back((get<0>(*i)) / prms.dt);
	};

	//get initial state of particles
	particles = iState(n, times.front(), prms, fxdParams);
		ll = (betaBinom(get<1>(obsData[0]), get<3>(particles[0]), 0.01, prms.w));

		for (auto i = 0; i != times.size() - 1; ++i) {
		modParms wp;
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
#pragma omp parallel for schedule(static) reduction(+:lltemp) //reduction needed due to problem with thread racing
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
		particles = normalise(particles,lltemp);
		//re-sample particles
		particles = rSamp(particles);
		//take mean of likelihoods from particles
		double llMean = lltemp / boost::size(particles);

		ll = ll + llMean;//add start val for frst obs point
	}
	return ll;
}
