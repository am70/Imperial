#include"lModH.h"


// proposal sd tuning function
// @param current s.d.
// @param target acceptance ratio
// @param current acceptance ratio
// @param maximum proposal s.d.
// @return tuned s.d.
double tuner(double curSd, double acptR, double curAcptR, double maxSd) {
	boost::math::normal dist(0.0, 1.0);
	if (curAcptR == 1);
	curAcptR = 0.99;
	if (curAcptR == 0);
	curAcptR = 0.01;
	double curSd = (curSd *quantile(dist, (acptR / 2))) / quantile(dist, (curAcptR / 2));
}


vector<pMMHres> mcmcSampler(
	vector<double> initParams,
	int fixedParam,
	vector<double> sdProps,
	vector<double> maxSddProps,
	int niter,
	int particles,
	int nburn,
	int monitoring,
	int startAdapt,
	int adptBurn,
	int	tell,
	vector<tuple<int, int>>	oDat) 
{
	

}