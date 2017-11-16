#include"lModH.h"
boost::mt19937 rng(std::time(0));
//boost::random::mt19937 engine{ boost::random::random_device{}() };



vector<int> rainfall1 = txtReader("Q:\\Imperial\\rf1.txt");
vector<int> rainfall2 = txtReader("Q:\\Imperial\\rf2.txt");
vector<int> rainfall3 = txtReader("Q:\\Imperial\\rf3.txt");




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
	double res = (curSd *quantile(dist, (acptR / 2))) / quantile(dist, (curAcptR / 2));
	if (res > maxSd) res = maxSd;
	return res;
}


//random num function 0-1
double rn01(void)
{
	static boost::uniform_01<boost::mt19937> zeroone(rng);
	return zeroone();
}


//structure updating function - could be improved?
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
	if (prmName == "z5")
		prms.z5 = propPrm;
	if (prmName == "z6")
		prms.z6 = propPrm;
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
	if (prmName == "sf5")
		prms.sf5 = propPrm;
	if (prmName == "sf6")
		prms.sf6 = propPrm;
	if (prmName == "n")
		prms.n = propPrm;
	return prms;
}

//log likelihood function 
//@param particles number of particles in pMMH
//@param prms current parameters
//@param obsDatX observed data
//@param fixedParam fixed parameters
//@return log likelihood value
double llFunc(int particles, modParms prms, obsDatX obsDat,int fixedParam) {
	vector<double> pfiltRes;

	for (auto j = 0; j != 6; ++j) {
		vector<tuple<int, int>> oDat;
		if (j == 0) {
			oDat = obsDat.garki101;
			prms.sf = pow(10,prms.sf1);
			prms.z = pow(10,prms.z1);
			prms.rF = rainfall1;
		}
		else if (j == 1) {
			oDat = obsDat.garki104;
			prms.sf = pow(10, prms.sf2);
			prms.z = pow(10, prms.z2);
			prms.rF = rainfall1;
		}
		else if (j == 2) {
			oDat = obsDat.garki219;
			prms.sf = pow(10, prms.sf3);
			prms.z = pow(10, prms.z3);
			prms.rF = rainfall2;
		}
		else if (j==3) {
			oDat = obsDat.garki220;
			prms.sf = pow(10, prms.sf4);
			prms.z = pow(10, prms.z4);
			prms.rF = rainfall2;
		}
		else if (j == 4) {
			oDat = obsDat.garki782;
			prms.sf = pow(10, prms.sf5);
			prms.z = pow(10, prms.z5);
			prms.rF = rainfall3;
		}
		else  {
			oDat = obsDat.garki783;
			prms.sf = pow(10, prms.sf6);
			prms.z = pow(10, prms.z6);
			prms.rF = rainfall3;
		}

		//run particle filter
		pfiltRes.emplace_back(pFilt(particles,
			oDat,//garki data
			prms,//parameters
			false,//full output or just likelihood
			fixedParam//fixed parameters
		));
	}
	return boost::accumulate(pfiltRes, 0.0);
}

//Parameter proposal function
//@param sd parameter standard deviation
//@param parameter
//@return proposed parameter
double propPrmFunc(double sd, double parm) {
	boost::normal_distribution<> nd(0.0, sd);
	double ran = nd(rng);
	double prop = ran + parm;
	if (prop < 0) prop = 0.0;
	return prop;
}

//log prior function
//@param current parameters
//@return sum loglikelihood for each parameter
double lprior(modParms prms) {
	double res = 0;
	boost::math::normal_distribution<double> d1(0.035, 0.009);//uoE
	res = res + (log(pdf(d1, prms.uoE)));
	boost::math::uniform_distribution<double> u1(0.001, 0.99);//uoE unif
	res = res + (log(pdf(u1, prms.uoE)));

	boost::math::normal_distribution<double> d2(0.035, 0.009);//uoL
	res = res + (log(pdf(d2, prms.uoL)));
	boost::math::uniform_distribution<double> u2(0.001, 0.99);//uoL unif
	res = res + (log(pdf(u2, prms.uoL)));

	boost::math::normal_distribution<double> d3(0.25, 0.07);//uP
	res = res + (log(pdf(d3, prms.uP)));
	boost::math::uniform_distribution<double> u3(0.001, 0.99);//uP unif
	res = res + (log(pdf(u3, prms.uP)));

	boost::math::normal_distribution<double> d4(13.06, 8);//Y
	res = res + (log(pdf(d4, prms.Y)));

	boost::math::uniform_distribution<double> u5(0.0001, 1);//w unif
	res = res + (log(pdf(u5, prms.w)));
	
	boost::math::uniform_distribution<double> u6(1, 7);//z1:4 unif
	res = res + (log(pdf(u6, prms.z1)));
	res = res + (log(pdf(u6, prms.z2)));
	res = res + (log(pdf(u6, prms.z3)));
	res = res + (log(pdf(u6, prms.z4)));
	res = res + (log(pdf(u6, prms.z5)));
	res = res + (log(pdf(u6, prms.z6)));


	boost::math::uniform_distribution<double> u7(1, 7);//sf1:4 unif
	res = res + (log(pdf(u7, prms.sf1)));
	res = res + (log(pdf(u7, prms.sf2)));
	res = res + (log(pdf(u7, prms.sf3)));
	res = res + (log(pdf(u7, prms.sf4)));
	res = res + (log(pdf(u7, prms.sf5)));
	res = res + (log(pdf(u7, prms.sf6)));


	boost::math::uniform_distribution<double> u8(10, 93.6);//n unif
	res = res + (log(pdf(u8, prms.n)));

	return(res);
}


//pMMH sampler function
//@param initParams initial parameter values
//@param fixedParam fixed parameter(s)
//@param sdProps starting s.d. for proposals for each parameter
//@param acptRs target acceptance ratios for each parameter
//@param fitPrms tuple containing name and values for each parameter 
//@param maxSdProps maximum s.d. for each parameter proposal to stop it going too high during tuning
//@param niter number of iterations in MMH algorithm
//@param nburn number of iterations to automatically burn
//@param monitoring if true print output every x number of iterations according to @param tell
//@param startAdapt when to start adapting proposal s.d.'s
//@param tell how often to print outputs
//@param oDat struct containing observed data
//@return struct containing results of pMMH
pMMHres pMMHSampler(
	modParms initParams,
	int fixedParam,
	vector<double> sdProps,
	vector<double> acptRs,
	vector<tuple<string, double>> fitPrms,
	vector<double> maxSdProps,
	int niter,
	int particles,
	int nburn,
	bool monitoring,
	int startAdapt,
	int	tell,
	obsDatX	oDat){
	int iter = 0;
	unsigned int parmNum = 0;
	double llProp;//proposed ll
	double llCur;//current ll
	vector<double> llRes;
	double llRatio; //ratio between proposed and current ll
	modParms prms = initParams;
	double propPrm; //proposed new parameter
	llCur = llFunc(particles, prms, oDat, fixedParam) + lprior(prms); //get value for initial ll
	vector<double> acptRcur = acptRs;//current acceptance ratio
	vector<double> acpts(acptRs.size(), 0.0);//number of acceptances (use acptRs to get correct vector length)
	vector<int> parmIter(acptRs.size(), 0);//iteration number for specific parameters
	pMMHres results;


	for (auto iter = 0; iter != niter; ++iter) {
		propPrm = propPrmFunc(sdProps[parmNum], get<1>(fitPrms[parmNum]));//propose new parameter
		prms = parmUpdt(prms, get<0>(fitPrms[parmNum]), propPrm);//update parameter struct
		llProp = llFunc(particles, prms, oDat, fixedParam) + lprior(prms);//find log likelihood from particle filter

		//print outputs
		if ((monitoring = true && iter % tell == 0)) {
			cout << "current = " << llCur << endl;
			cout << "proposed = " << llProp << endl;
			cout << "iteration " << iter << " of " << niter << endl;
			cout << " uoE = " << prms.uoE << " uoL = " << prms.uoL << " uoP = " << prms.uP << " Y = " << prms.Y << " w = " << prms.w << " n = " << prms.n << " z1 = " << prms.z1 << endl
				<< " z2 = " << prms.z2 << " z3 = " << prms.z3 << " z4 = " << prms.z4 << " z5 = " << prms.z5 << " z6 = " << prms.z6 << " sf1 = " << prms.sf1 << " sf2 = " << prms.sf2 << " sf3 = " << prms.sf3 << " sf4 = " << prms.sf4 << " sf5 = " << prms.sf5 << " sf6 = " << prms.sf6 << endl;
			cout << "||---------aratio--------||" << endl;
			cout << " uoE = " << acptRcur[0] << " uoL = " << acptRcur[1] << " uoP = " << acptRcur[2] << " Y = " << acptRcur[3] << " w = " << acptRcur[4] << " n = " << acptRcur[5] << " z1 = " << acptRcur[6] << endl
				<< " z2 = " << acptRcur[7] << " z3 = " << acptRcur[8] << " z4 = " << acptRcur[9] << " z5 = " << acptRcur[10] << " z6 = " << acptRcur[11] << " sf1 = " << acptRcur[12] << " sf2 = " << acptRcur[13] << " sf3 = " << acptRcur[14] << " sf4 = " << acptRcur[15] << " sf5 = " << acptRcur[16] << " sf6 = " << acptRcur[17] << endl;
			cout << "||---------aNum--------||" << endl;
			cout << " uoE = " << acpts[0] << " uoL = " << acpts[1] << " uoP = " << acpts[2] << " Y = " << acpts[3] << " w = " << acpts[4] << " n = " << acpts[5] << " z1 = " << acpts[6] << endl
				<< " z2 = " << acpts[7] << " z3 = " << acpts[8] << " z4 = " << acpts[9] << " z5 = " << acpts[10] << " z6 = " << acpts[11] << " sf1 = " << acpts[12] << " sf2 = " << acpts[13] << " sf3 = " << acpts[14] << " sf4 = " << acpts[15] << " sf5 = " << acpts[16] << " sf6 = " << acpts[17] << endl;
			cout << "||---------sd--------||" << endl;
			
			cout << " uoE = " << sdProps[0] << " uoL = " << sdProps[1] << " uoP = " << sdProps[2] << " Y = " << sdProps[3] << " w = " << sdProps[4] << " n = " << sdProps[5] << " z1 = " << sdProps[6] << endl
				<< " z2 = " << sdProps[7] << " z3 = " << sdProps[8] << " z4 = " << sdProps[9] <<" z5 = " << sdProps[10] << " z6 = " << sdProps[11] << " sf1 = " << sdProps[12] << " sf2 = " << sdProps[13] << " sf3 = " << sdProps[14] << " sf4 = " << sdProps[15] << " sf5 = " << sdProps[16] << " sf6 = " << sdProps[17] << endl;
			cout << "||-----------------------||" << endl;
			cout << "||-----------------------||" << endl;
		}

		llRatio = llProp - llCur;//ratio between previous ll and proposed ll

		if (llRatio >= 0 || rn01() <= exp(llRatio)) {
			llCur = llProp; //update current ll if ll is accepted
			acpts[parmNum] = acpts[parmNum] + 1;
			get<1>(fitPrms[parmNum]) = propPrm;
		}
		else {
			prms = parmUpdt(prms, get<0>(fitPrms[parmNum]), get<1>(fitPrms[parmNum])); //change parameter back if not accepted
		}

		if (iter > startAdapt) {
			//start adapting with tuner 
			acptRcur[parmNum] = (acpts[parmNum] / parmIter[parmNum]); //calc current acceptance ratio for parameter
			sdProps[parmNum] = tuner(sdProps[parmNum], acptRs[parmNum], acptRcur[parmNum], maxSdProps[parmNum]);
		}

		//if passed burn-in, start adding to results
		if (iter > nburn) {
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
			results.z5.emplace_back(prms.z5);
			results.z6.emplace_back(prms.z6);
			results.sf1.emplace_back(prms.sf1);
			results.sf2.emplace_back(prms.sf2);
			results.sf3.emplace_back(prms.sf3);
			results.sf4.emplace_back(prms.sf4);
			results.sf5.emplace_back(prms.sf5);
			results.sf6.emplace_back(prms.sf6);
			results.ll.emplace_back(llCur);
		}

		parmNum++; //parameter number, if greater than number of parms, revert back to 0
		if (parmNum >= boost::size(sdProps)) {
			parmNum = 0;
		}



		parmIter[parmNum]++;
	}

	return results;
}