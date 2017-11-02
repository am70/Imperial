#include"lModH.h"

vector<int> rainfall1 = txtReader("Q:\\Imperial\\rf1.txt");
vector<int> rainfall2 = txtReader("Q:\\Imperial\\rf2.txt");

// proposal sd tuning function
// @param current s.d.
// @param target acceptance ratio
// @param current acceptance ratio
// @param maximum proposal s.d.
// @return tuned s.d.
double tuner(double curSd, double acptR, double curAcptR, double maxSd) {
	boost::math::normal dist(0.0, 1.0);
	if (curAcptR == 1)
	curAcptR = 0.99;
	if (curAcptR == 0)
	curAcptR = 0.01;
	return (curSd *quantile(dist, (acptR / 2))) / quantile(dist, (curAcptR / 2));
}


//random num 0-1
double rn01(void)
{
	boost::mt19937 rng(std::time(0));
	static boost::uniform_01<boost::mt19937> zeroone(rng);
	return zeroone();
}


//horrible structure updating function that can probably be better somehow...
modParms parmUpdt(modParms prms, string prmName, double propPrm) {
	if (prmName == "uoE")
		prms.uoE = propPrm;
	if (prmName == "uoL")
		prms.uoL = propPrm;
	if (prmName == "uP")
		prms.uP = propPrm;
	if (prmName == "Y")
		prms.Y = propPrm;
	if (prmName == "z1")
		prms.z1 = propPrm;
	if (prmName == "z2")
		prms.z2 = propPrm;
	if (prmName == "z3")
		prms.z3 = propPrm;
	if (prmName == "z4")
		prms.z4 = propPrm;
	if (prmName == "w")
		prms.w = propPrm;
	if (prmName == "sf1")
		prms.sf1 = propPrm;
	if (prmName == "sf2")
		prms.sf2 = propPrm;
	if (prmName == "sf3")
		prms.sf3 = propPrm;
	if (prmName == "sf4")
		prms.sf4 = propPrm;
	if (prmName == "n")
		prms.n = propPrm;
	return prms;
}


double llFunc(int particles, modParms prms, obsDatX obsDat) {
	for (auto j = 0; j != 4; ++j) {
		vector<tuple<int, int>> oDat;

		if (j == 0) {
			oDat = obsDat.garki101;
			prms.sf = prms.sf1;
			prms.z = prms.z1;
			prms.rF = rainfall1;
		}
		else if (j == 1) {
			oDat = obsDat.garki104;
			prms.sf = prms.sf2;
			prms.z = prms.z2;
			prms.rF = rainfall1;
		}
		else if (j == 2) {
			oDat = obsDat.garki219;
			prms.sf = prms.sf3;
			prms.z = prms.z3;
			prms.rF = rainfall2;
		}
		else oDat = obsDat.garki220;
		prms.sf = prms.sf4;
		prms.z = prms.z4;
		prms.rF = rainfall2;

	
		//run particle filter
		return pFilt(particles,
			oDat,//garki data
			prms,//parameters
			false,//full output or just likelihood
			7//fixed parameters
		);
	}
}


double propPrmFunc(double sd, double parm) {
	boost::mt19937 rng(std::time(0));
	static boost::normal_distribution<> nd(0.0, sd);
	boost::variate_generator<boost::mt19937&,
    boost::normal_distribution<> > var_nor(rng, nd);
	return var_nor()+parm;
}


pMMHres pMMHSampler(
	modParms initParams,
	int fixedParam,
	vector<double> sdProps,
	vector<double> acptRs,
	vector<tuple<string,double>> fitPrms,
	vector<double> maxSdProps,
	int niter,
	int particles,
	int nburn,
	int monitoring,
	int startAdapt,
	int	tell,
	obsDatX	oDat)
{
	int iter = 0;
	unsigned int parmNum = 0;
	double llProp;//proposed ll
	double llCur;//current ll
	vector<double> llRes;

	double llRatio; //ratio between proposed and current ll
	modParms prms = initParams;

	double propPrm; //proposed new parameter
	llCur = llFunc(particles, prms, oDat); //get value for initial ll

	vector<double> acptRcur = acptRs;
	pMMHres results;

	
	while (iter < niter) {

		propPrm = propPrmFunc(sdProps[parmNum], get<1>(fitPrms[parmNum]));//propose new parameter
		prms = parmUpdt(prms, get<0>(fitPrms[parmNum]), propPrm);//update parameter struct

		llProp = llFunc(particles,prms, oDat);//find log likelihood from particle filter

		llRatio = llProp - llCur;//ratio between previous ll and proposed ll

		if (llRatio >= 0 || rn01() <= exp(llRatio)) { 
			llCur = llProp; //update current ll if ll is accepted
		} 
		else { 
			prms = parmUpdt(prms, get<0>(fitPrms[parmNum]), get<1>(fitPrms[parmNum])); //change parameter back if not accepted
		} 


		if (iter > startAdapt) {
			//start adapting with tuner 
			sdProps[parmNum] = tuner(sdProps[parmNum], acptRs[parmNum], acptRcur[parmNum], maxSdProps[parmNum]);
		}

		if (iter > nburn) { //if passed burnin, start adding to results
			results.uoE.emplace_back(prms.uoE);
			results.uoL.emplace_back(prms.uoL);
			results.uP.emplace_back(prms.uP);
			results.Y.emplace_back(prms.Y);
			results.w.emplace_back(prms.w);
			results.n.emplace_back(prms.n);
			results.z1.emplace_back(prms.z1);
			results.z2.emplace_back(prms.z2);
			results.z3.emplace_back(prms.z3);
			results.z4.emplace_back(prms.z4);
			results.sf1.emplace_back(prms.sf1);
			results.sf2.emplace_back(prms.sf2);
			results.sf3.emplace_back(prms.sf3);
			results.sf4.emplace_back(prms.sf4);
			results.ll.emplace_back(llCur);
		}

		parmNum++; //parameter number, if greater than number of parms, revert back to 0
		if (parmNum >= boost::size(sdProps)) {
			parmNum = 0;
		}
		
		acptRcur[parmNum] = (acptRcur[parmNum] / iter); //calc current acceptance ratio for parameter

		if ((monitoring > 1 && iter % tell == 0)) {
			cout << "iteration " << iter << endl;
			cout << "current = " << llCur << endl;//print proposed ll
			cout << "proposed = " << llProp << endl;//print proposed ll
		}

		iter++;//increase iteration
	}
	return results;
}