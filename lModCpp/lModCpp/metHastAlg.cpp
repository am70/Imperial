#include"lModH.h"
#include <algorithm>
boost::mt19937 rng(std::time(0));


vector<double> rainfall_05 = txtReader("\\\\qdrive.dide.ic.ac.uk\\homes\\ALM210\\Imperial\\lModCpp\\Data\\rf05.txt", 0.25);
vector<double> rainfall_07 = txtReader("\\\\qdrive.dide.ic.ac.uk\\homes\\ALM210\\Imperial\\lModCpp\\Data\\rf07.txt", 0.25);
vector<double> rainfall_08 = txtReader("\\\\qdrive.dide.ic.ac.uk\\homes\\ALM210\\Imperial\\lModCpp\\Data\\rf08.txt", 0.25);
vector<double> rainfall_04 = txtReader("\\\\qdrive.dide.ic.ac.uk\\homes\\ALM210\\Imperial\\lModCpp\\Data\\rf04.txt", 0.25);
vector<double> rainfall_02 = txtReader("\\\\qdrive.dide.ic.ac.uk\\homes\\ALM210\\Imperial\\lModCpp\\Data\\rf02.txt", 0.25);
vector<double> rainfall_01 = txtReader("\\\\qdrive.dide.ic.ac.uk\\homes\\ALM210\\Imperial\\lModCpp\\Data\\rf01.txt", 0.25);
vector<double> rainfall_03 = txtReader("\\\\qdrive.dide.ic.ac.uk\\homes\\ALM210\\Imperial\\lModCpp\\Data\\rf03.txt", 0.25);



/*proposal sd tuning function
@param current s.d.
@param target acceptance ratio
@param current acceptance ratio
@param maximum proposal s.d.
@return tuned s.d.*/
double tuner(double curSd, double acptR, double curAcptR, double maxSd){
	boost::math::normal dist(0.0, 1.0);
	if (curAcptR == 1)
		curAcptR = 0.99;
	if (curAcptR == 0)
		curAcptR = 0.01;
	double res = (curSd *quantile(dist, (acptR / 2))) / quantile(dist, (curAcptR / 2));
	if (res > maxSd) res = maxSd;
	if (res <= 0) res = 0;
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
	if (prmName == "dE")
		prms.dE = propPrm;
	if (prmName == "dL")
		prms.dL = propPrm;
	if (prmName == "dP")
		prms.dP = propPrm;
	if (prmName == "o")
		prms.o = propPrm;
	if (prmName == "uM")
		prms.uM = propPrm;
	if (prmName == "Mg")
		prms.Mg = propPrm;
	if (prmName == "p")
		prms.p = propPrm;
	if (prmName == "tau")
		prms.tau = propPrm;
	return prms;
}

/*log likelihood function
@param particles number of particles in pMMH
@param prms current parameters
@param obsDatX observed data
@param fixedParam fixed parameters
@return log likelihood value*/
double llFunc(int particles, modParms prms, obsDatX obsDat, string dFunc, std::vector<int> seeds) {
	vector<double> pfiltRes;

	for (auto j = 0; j != 3; ++j) {

		vector<tuple<int, int>> oDat;
		//if (j == 0) {
		//	oDat = obsDat.garki408;
		//	prms.sf = pow(10, prms.sf1);
		//	prms.z = pow(10, prms.z1);
		//	prms.rF = rainfall_03;

		//}
		//else if (j == 1) {
		//	oDat = obsDat.garki154;
		//	prms.sf = pow(10, prms.sf2);
		//	prms.z = pow(10, prms.z2);
		//	prms.rF = rainfall_05;
		//}
		// if (j == 0) {
		//	oDat = obsDat.garki408;
		//	prms.sf = pow(10, prms.sf3);
		//	prms.z = pow(10, prms.z3);
		//	prms.rF = rainfall_03;
		//}
		 if (j == 0) {
			oDat = obsDat.garki801;
			prms.sf = pow(10, prms.sf4);
			prms.z = pow(10, prms.z4);
			prms.rF = rainfall_01;
		}
		else if (j == 1) {
			oDat = obsDat.garki553;
			prms.sf = pow(10, prms.sf5);
			prms.z = pow(10, prms.z5);
			prms.rF = rainfall_02;
		}
		else if (j == 2) {
			oDat = obsDat.garki802;
			prms.sf = pow(10, prms.sf6);
			prms.z = pow(10, prms.z6);
			prms.rF = rainfall_01;
		}

		//run particle filter

		pfiltRes.emplace_back(pFilt(particles,
			oDat,//garki data
			prms,//parameters
			false,//full output or just likelihood
			"\\\\qdrive.dide.ic.ac.uk\\homes\\ALM210\\Imperial\\fitPlots\\test.txt",
			false,
			dFunc//density function 
			,seeds
		));
	}

	return boost::accumulate(pfiltRes, 0.0);
}

/*Parameter proposal function
@param sd parameter standard deviation
@param parameter
@return proposed parameter*/
double propPrmFunc(double sd, double parm) {
	boost::normal_distribution<> nd(0.0, sd);
	double ran = nd(rng);
	double prop = ran + parm;
	if (prop <= 0) prop = 0;
	return prop;
}

/*log prior function
@param current parameters
@return sum loglikelihood for each parameter*/
double lprior(modParms prms) {
	double res = 0;
	boost::math::normal_distribution<double> d1(0.035, 0.0056);//uoE
	res = res + (log(pdf(d1, prms.uoE)));
	boost::math::uniform_distribution<double> u1(0.001, 0.99);//uoE unif
	res = res + (log(pdf(u1, prms.uoE)));

	boost::math::normal_distribution<double> d2(0.035, 0.0056);//uoL
	res = res + (log(pdf(d2, prms.uoL)));
	boost::math::uniform_distribution<double> u2(0.001, 0.99);//uoL unif
	res = res + (log(pdf(u2, prms.uoL)));

	boost::math::normal_distribution<double> d3(0.25, 0.05);//uP
	res = res + (log(pdf(d3, prms.uP)));
	boost::math::uniform_distribution<double> u3(0.001, 0.99);//uP unif
	res = res + (log(pdf(u3, prms.uP)));

	boost::math::normal_distribution<double> uM1(0.091,0.005);//uM unif
	res = res + (log(pdf(uM1, prms.uM)));


	boost::math::normal_distribution<double> d4(13.06, 4);//Y
	res = res + (log(pdf(d4, prms.Y)));

	boost::math::uniform_distribution<double> d4x(0.1, 70);//Y
	res = res + (log(pdf(d4x, prms.Y)));

	//boost::math::uniform_distribution<double> u5(0.000001, 1);//w unif
	//res = res + (log(pdf(u5, prms.w)));
	
	boost::math::uniform_distribution<double> u61(1, 15);//z1 unif
	boost::math::uniform_distribution<double> u62(1, 15);//z2 unif
	boost::math::uniform_distribution<double> u63(1, 10);//z3 unif
	boost::math::uniform_distribution<double> u64(1, 15);//z4 unif
	boost::math::uniform_distribution<double> u65(1, 10);//z5 unif
	boost::math::uniform_distribution<double> u66(1, 15);//z6 unif

	res = res + (log(pdf(u61, prms.z1)));
	res = res + (log(pdf(u62, prms.z2)));
	res = res + (log(pdf(u63, prms.z3)));
	res = res + (log(pdf(u64, prms.z4)));
	res = res + (log(pdf(u65, prms.z5)));
	res = res + (log(pdf(u66, prms.z6)));



	boost::math::uniform_distribution<double> u71(1, 10);//sf1 unif
	boost::math::uniform_distribution<double> u72(1, 10);//sf2 unif
	boost::math::uniform_distribution<double> u73(1, 10);//sf3 unif
	boost::math::uniform_distribution<double> u74(1, 10);//sf4 unif
	boost::math::uniform_distribution<double> u75(1, 10);//sf5 unif
	boost::math::uniform_distribution<double> u76(1, 10);//sf6 unif

	res = res + (log(pdf(u71, prms.sf1)));
	res = res + (log(pdf(u72, prms.sf2)));
	res = res + (log(pdf(u73, prms.sf3)));
	res = res + (log(pdf(u74, prms.sf4)));
	res = res + (log(pdf(u75, prms.sf5)));
	res = res + (log(pdf(u76, prms.sf6)));

	//dE = 0.15, dL = 0.269, dP = 1.563

	boost::math::uniform_distribution<double> u8(0.01, 93);//n unif
	res = res + (log(pdf(u8, prms.n)));



	boost::math::normal_distribution<double> d5(0.150602, 0.04);//dE
	res = res + (log(pdf(d5, prms.dE)));

	boost::math::normal_distribution<double> d6(0.268812, 0.06);//dL
	res = res + (log(pdf(d6, prms.dL)));

	boost::math::normal_distribution<double> d7(1, 0.1);//dP
	res = res + (log(pdf(d7, prms.dP)));

	boost::math::normal_distribution<double> dtau(7, 1.5);//tau
	res = res + (log(pdf(dtau, rint(prms.tau))));


	boost::math::uniform_distribution<double> d8(0.001,10);//o
	res = res + (log(pdf(d8, prms.o)));


	boost::math::uniform_distribution<double> mm(1, 10);//mg
	res = res + (log(pdf(mm, prms.Mg)));

	return(res);
}



std::vector<int> gen_seeds(int num) {
	// start by stuffing them into a set to guarantee uniqueness.
	std::set<int> s;

	std::random_device g;

	while (s.size() < num)
		s.insert(g());

	// then return them in a vector to give random access:
	std::vector<int> seeds(s.begin(), s.end());
	return seeds;
}



/*pMMH sampler function
@param initParams initial parameter values
@param fixedParam fixed parameter(s)
@param sdProps starting s.d. for proposals for each parameter
@param acptRs target acceptance ratios for each parameter
@param fitPrms tuple containing name and values for each parameter 
@param maxSdProps maximum s.d. for each parameter proposal to stop it going too high during tuning
@param niter number of iterations in MMH algorithm
@param nburn number of iterations to automatically burn
@param monitoring if true print output every x number of iterations according to @param tell
@param startAdapt when to start adapting proposal s.d.'s
@param tell how often to print outputs
@param oDat struct containing observed data@return struct containing results of pMMH*/
pMMHres pMMHSampler(
	modParms initParams,
	string dFunc,
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
	int seedIter;

	if (dFunc != "expNoClumped" && dFunc != "linearNoClumped" && dFunc != "powerNoClumped" && dFunc != "expClumped" && dFunc != "linearClumped" && dFunc != "powerClumped") {
		cerr << "dFunc must equal correct value: linearNoClumped, powerNoClumped, expNoClumped, linearClumped, powerClumped or expClumped";
		cin.get();
	}

	std::vector<int> seedStore = gen_seeds(particles*niter);
	vector<int> seeds(seedStore.begin(), seedStore.begin() + particles);

	llCur = llFunc(particles, prms, oDat, dFunc, seeds) + lprior(prms); //get value for initial ll
	vector<double> acptRcur = acptRs;//current acceptance ratio
	vector<double> acpts(acptRs.size(), 0.0);//number of acceptances (use acptRs to get correct vector length)
	vector<int> parmIter(acptRs.size(), 0);//iteration number for specific parameters
	pMMHres results;



	for (auto iter = 0; iter != niter; ++iter) {

		if (maxSdProps[parmNum] != 0) {

			propPrm = propPrmFunc(sdProps[parmNum], get<1>(fitPrms[parmNum]));//propose new parameter

			prms = parmUpdt(prms, get<0>(fitPrms[parmNum]), propPrm);//update parameter struct

			if ((monitoring = true && iter % tell == 0)) {

				cout << endl << "||-----------------------||" << dFunc << "||-----------------------||" << endl;
				cout << "iteration " << iter << " of " << niter << endl;
				cout << " uoE = " << prms.uoE << " uoL = " << prms.uoL << " uP = " << prms.uP << " uM = " << prms.uM << " Y = " << prms.Y << " w = " << prms.w << " n = " << prms.n << " z1 = " << prms.z1 << endl
					<< " z2 = " << prms.z2 << " z3 = " << prms.z3
					<< " z4 = " << prms.z4 << " z5 = " << prms.z5 << " z6 = " << prms.z6 << " sf1 = " << prms.sf1 << " sf2 = " << prms.sf2 << " sf3 = " << prms.sf3 << " sf4 = " << prms.sf4 << " sf5 = " << prms.sf5 << " sf6 = " << prms.sf6
					<< "dE = " << prms.dE << " dL = " << prms.dL << " dP = " << prms.dP << " o = " << prms.o << " Mg = " << prms.Mg << " p = " << prms.p << " tau = " << prms.tau << endl;
				cout << "||---------aratio--------||" << endl;


				for (auto iter = 0; iter != size(sdProps); ++iter) {
					cout << get<0>(fitPrms[iter]) << " = " << acptRcur[iter] << " ";
				}
				cout << endl << "||---------aNum--------||" << endl;
				for (auto iter = 0; iter != size(sdProps); ++iter) {
					cout << get<0>(fitPrms[iter]) << " = " << acpts[iter] << " ";
				}
				cout << endl << "||---------sd--------||" << endl;
				for (auto iter = 0; iter != size(sdProps); ++iter) {
					cout << get<0>(fitPrms[iter]) << " = " << sdProps[iter] << " ";
				}
				cout << endl << "||-----------------------||" << endl;
			}


			/*ofstream myfile;
			myfile.open("Q:\\test.csv");
				myfile << " Mg = " << prms.Mg << endl;*/

			vector<int> seeds(seedStore.begin()+ seedIter, seedStore.begin()+ seedIter + particles);
			llProp = llFunc(particles, prms, oDat, dFunc, seeds) + lprior(prms);//find log likelihood from particle filter
			seedIter = seedIter + particles;

			//boost::math::normal_distribution<double> curr(get<1>(fitPrms[parmNum]), sdProps[parmNum]);//normal dist for metHast 
			//boost::math::normal_distribution<double> prp(propPrm, sdProps[parmNum]);//normal dist for metHast 

			//llRatio = (llProp -log(pdf(curr, propPrm))) - (llCur- log(pdf(prp, get<1>(fitPrms[parmNum]))));//ratio between previous ll and proposed ll
			double llRatio = llProp - llCur;//ratio between previous ll and proposed ll


			//print outputs
			if ((monitoring = true && iter % tell == 0)) {
				cout << "current = " << llCur << endl;
				cout << "proposed = " << llProp << endl;
				cout << "llratio = " << llRatio << endl;
				cout << "current parameter = " << get<0>(fitPrms[parmNum]) << endl;

			}


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
				results.uM.emplace_back(prms.uM);

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

				results.dE.emplace_back(prms.dE);
				results.dL.emplace_back(prms.dL);
				results.dP.emplace_back(prms.dP);
				results.o.emplace_back(prms.o);
				results.Mg.emplace_back(prms.Mg);
				results.p.emplace_back(prms.p);
				results.tau.emplace_back(prms.tau);


				results.ll.emplace_back(llCur);
			}

			parmNum++;

			; //parameter number, if greater than number of parms, revert back to 0
			if (parmNum >= boost::size(sdProps)) {
				parmNum = 0;
			}

			parmIter[parmNum]++;
		}
		else {
			parmNum++;
			parmIter[parmNum]++;
		}

	}

	return results;
}