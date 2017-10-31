#include "lModH.h"
boost::mt19937 mrand(std::time(0));

modParms tParms;//parameter struct to store new param values


// initial state sampler, samples random initial states # need to fit initial E and use this to inform L, P & M
// @param N number of particles for particle filter
// @param t starting time period for rainfall
// @param prms model parameters
// @fxdParams fixed parameter (currently just for n in mosParamsP) - maybe update for multiple fixed parameters
// @return conditions for E, L, P and M, double is empty for addition of weight later. 
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
	vector<int> rF = iParms.rF;
	double uM = iParms.uM;
	double uP = iParms.uP;
	double tr = iParms.tr / iParms.dt;
	double B = iParms.B;
	double K;
	double uE;
	double uL; 
	int E;
	int L;
	int P;
	int M;
	double a;
	//int t = iParms.startTime;
	double trx = iParms.tr / iParms.dt;

	

	int rFsum = std::accumulate(rF.begin() + (time - trx), rF.begin() + time, 0);
	K = (1 + (sf*((1 / trx)*rFsum)));

	a = (0.5 * dL*dP) / (uM*(dP + uP));

	uE = UoE*exp(z / K);
	uL = UoL*exp(y*z / K);


	E = ((B*a)*z) / (dE + (B*a) + uE);
	L = (dE*z) / (dE + dL + uL);
	M = (1 / 2 * dL*dP*L) / uM*(dP + uP);
	P = 2 * uM*M / dP;
	
	vector<tuple<int, int, int, int, double>> states = { { E, L, P, M ,0.0} };
	states.resize(N, { E, L, P, M, 0.0 });
	return states;
}

//binomial function
//@param n number of trials
//@param p probability of success
//@return number of successes
int binom(int n, double p) {
	boost::binomial_distribution<int> distribution(n, p);
	return  distribution(mrand);
}

//random poisson draw
//@param lambda
//@return random poisson draw
luint rpois(luint lambda)
{
	if (lambda > 0) {
		boost::poisson_distribution<luint> d(lambda);
		return d(mrand);
	}
	else return (0);
}

//lbeta function
//@param a
//@param b
//@return log beta
double lbeta(double a, double b) {
	return double(lgamma(a)+lgamma(b)) - lgamma(a + b);
}

// Beta binomial likelihood function
// @param k Observed data point
// @param n Simulated data point
// @param p fraction of population
// @param w overdisperion parameter
// @return log likelihood
double betaBinom(double k, double n, double p, double w) {
	double a = p * ((1 / w) - 1);
	double b = (1 - p) * ((1 / w) - 1);
	return lbeta(k + a, n - k + b) - lbeta(a, b) + log(boost::math::binomial_coefficient<double>(n, k));
}


//model step function
//@param wp a structure containing start times and states for model step
//@param obsData, observed data for calculating likelihood value
//@return tuple containing state ints for E, L, P, M and likelihood (non-log), at the end of 
//a model run with designated start and end times
tuple<int, int, int, int, double> modStepFnc(wpStruct wp, int obsData) {
	vector<int> rFc;
	vector<tuple<int, int, int, int>> modRun;
	double weight;

	rFc = wp.rFclust;

	
	tParms.E0 = wp.E0;
	tParms.L0 = wp.L0;
	tParms.P0 = wp.P0;
	tParms.M0 = wp.M0;
	tParms.startTime = wp.startTime;
	tParms.endTime = wp.endTime;
	tParms.rF = rFc;
	modRun = mPmod(tParms);



	weight = betaBinom(obsData, 1+get<3>(modRun.back()), 0.01, wp.w); //add weights to tuple
	tuple<int, int, int, int, double> res = {get<0>(modRun.back()),get<1>(modRun.back()),get<2>(modRun.back()),get<3>(modRun.back()), weight};


	return res;
}

//rand sample function
auto rSamp(std::vector<std::tuple<int, int, int, int, double>>& samp)
{
	std::vector<std::tuple<int,int,int,int, double>> temp;           
	std::vector<double> v;                                  // weights
	for (auto&& i : samp)
		v.push_back(exp(std::get<4>(i)));                        // get exponential weights
	std::discrete_distribution<int> dd{ v.begin(), v.end() }; // create distribution
	//static std::random_device rd;
	for (size_t i{}; i < samp.size(); ++i)
		temp.push_back(samp[dd(mrand)]);                        // build return vector by selecting according to weights
	return temp;
}


// Particle filter function - runs model in steps between each observed data point
// before using weighted re - sampling of model state at each step to start next
// @param n number of particles 
// @param iState function for calculating initial state of particles
// @param obsData observed data to fit model to
// @param prms model parameters
// @param resM True/False whether to output likelihood values or results of simulation for plotting figures
// @param fxedParams fixed parameters - currently just accepting tr but may expand to include more than one parameter
// @return if resM = F log likelihood value, if resM = T mean results of simulation
double pFilt(int n,
	vector< tuple<int, int> > obsData,
	modParms prms,
	bool resM,
	int rFclust,
	int fxdParams) {

	
	vector<int> times;
	vector<double> ll; //vector of log likelihood values
	vector<double> lltemp;
	vector<tuple<int, int, int, int, double>>  particles;
	particles.reserve(n);
	for (auto i = begin(obsData); i != end(obsData); ++i) {
		//calculate discreet time steps for observed data
		times.emplace_back((get<0>(*i))/prms.dt);
	};

	

	particles = iState(n,times.front(),prms,fxdParams);

	for (auto i = 0; i != times.size()- 1; ++i) {
		int loc = 0;
		int run = 0;

		wpStruct wp;
		wp.uoE = prms.uoE;
		wp.uoL = prms.uoL;
		wp.uP = prms.uP;
		wp.Y = prms.Y;
		wp.sf = prms.sf;
		wp.n = prms.n;
		wp.w = prms.w;
		wp.rFclust = prms.rF;
		wp.fxdPrm = fxdParams;
		wp.startTime = times.at(i);
		wp.endTime = times.at(i+1);
		

		//run model in step and update particles - should be in parallel
		lltemp.clear();
		for (auto j = begin(particles); j != end(particles); ++j) {
			wp.E0 = get<0>((*j));
			wp.L0 = get<1>((*j));
			wp.P0 = get<2>((*j));
			wp.M0 = get<3>((*j));
			int obsDatPoint = get<1>(obsData[i]);
			particles.at(loc) = modStepFnc(wp, obsDatPoint);
			lltemp.emplace_back(get<4>(particles.at(loc)));
			loc++;
		}
		//re-sample particles
		particles = rSamp(particles);

		//take mean of likelihoods from particles
		double llMean = (boost::accumulate(lltemp, 0.0))/ lltemp.size();

		ll.emplace_back(llMean);
	
	}
return (boost::accumulate(ll, 0.0));
}
