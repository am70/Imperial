#include "lModH.h"
boost::mt19937 mrand(std::time(0));
vector<int> rainfall1 = txtReader("Q:\\Imperial\\rf1.txt");
vector<int> rainfall2 = txtReader("Q:\\Imperial\\rf2.txt");
modParms tParms;//parameter struct to store new param values


// initial state sampler, samples random initial states # need to fit initial E and use this to inform L, P & M
// @param N number of particles for particle filter
// @param t starting time period for rainfall
// @param prms model parameters
// @fxdParams fixed parameter (currently just for n in mosParamsP) - maybe update for multiple fixed parameters
// @return conditions for E, L, P and M
vector<tuple<int, int, int, int>> iState(int N, t, modParms iParms, int fxdParm,) {
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
	double a;
	double uE;
	double uL; 
	int E;
	int L;
	int P;
	int M;
	int t = iParms.startTime;
	double trx = iParms.tr / iParms.dt;

	int rFsum = std::accumulate(rF.begin() + (t - trx), rF.end() - t, 0);
	K = (1 + (sf*((1 / trx)*rFsum)));

	a = (1 / 2 * dL*dP) / (uM*(dP + uP));

	uE = UoE*exp(z / K);
	uL = UoL*exp(y*z / K);


	E = ((B*a)*z) / (dE + (B*a) + uE);
	L = (dE*z) / (dE + dL + uL);
	M = (1 / 2 * dL*dP*L) / uM*(dP + uP);
	P = 2 * uM*M / dP;

	vector<tuple<int, int, int, int>> states = { { E, L, P, M } };
	states.resize(N, { E, L, P, M });
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
	return double(lgamma(a)*lgamma(b)) / lgamma(a + b);
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
	lbeta(k + a, n - k + b) - lbeta(a, b) + log(boost::math::binomial_coefficient<double>(n, k));
}


//model step function
//@param wp a structure containing start times and states for model step
//@return tuple containing state ints for E, L, P & M at the end of 
//a model run with designated start and end times
tuple<int, int, int, int> modStepFnc(wpStruct wp) {
	vector<int> rFc;
	vector<tuple<int, int, int, int>> modRun;
	if (wp.rFclust == 1)
		rFc = rainfall1;
	else
		rFc = rainfall2;


	tParms.E0 = wp.E0;
	tParms.L0 = wp.L0;
	tParms.P0 = wp.P0;
	tParms.M0 = wp.M0;
	tParms.startTime = wp.startTime;
	tParms.endTime = wp.endTime;
	tParms.rF = rFc;

	modRun = mPmod(tParms);

	tuple<int,int,int,int> res = modRun.back();
	return res;
}

//particle filter
vector<int> pFilt(int n,
	vector<int> iState,
	vector< tuple<int, int> > obsData,
	modParms prms,
	string resM,
	string rFclust,
	int fxdParams) {
	//read times data
	vector<int> times;
	vector<tuple<int,int,int,int>> particles = iState(n, times.front, prms, fxdParams) //initial state


	//for (auto i = begin(obsData); i != end(obsData); ++i) {
	//	times.push_back << get<0>(*i) / delta;
	//}

	//particles = iState(n, t = times[1], prms = prms, fxdParams) #initial state

}
