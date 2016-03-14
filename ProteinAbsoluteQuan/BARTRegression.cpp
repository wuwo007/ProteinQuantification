#include"stdafx.h"
#include"BARTRegression.h"
#include"ProteinWorker.h"
using std::cout;
using std::endl;
void CBART::mf_xitofile(std::ofstream& fnm, xinfo xi)
{
	fnm << xi.size() << endl;
	for (uint i = 0; i< xi.size(); i++) {
		fnm << xi[i].size() << endl;
		for (uint j = 0; j<xi[i].size(); j++) fnm << xi[i][j] << "  ";
		fnm << endl;
	}
	fnm << endl;
}
// construct the train set 
bool CBART::mf_ConstructXY(std::ofstream& log, vector<CProtein> proteins )
{

	cout << "\tConstructing XY\n";
	log << "\tConstructing XY\n";
	
	m_dMiny = INFINITY; //use range of y to calibrate prior for bottom node mu's
	m_dMaxy = -INFINITY;
	vector<CProtein>::iterator  ProteinIter;
	int j = 0; //feature num
	int n = 0; //peptide num per protein
	double dytemp;

	m_iTrainSampleNumber = 0;  //Change i to m_NumberOfLimitedSamples by gzhq 20150528
	m_iNumberOfTrainProteins = 0;

	m_vTrainx.clear();
	m_vTrainy.clear();
	m_cAllys.clear();
	
	for (ProteinIter = proteins.begin(); ProteinIter != proteins.end(); ProteinIter++)
	{
		if (ProteinIter->m_iPeptidesNumber >= UniquePeptidesTrainThreshold)
		{
			m_iNumberOfTrainProteins++;
			for (n = 0; n < ProteinIter->m_iPeptidesNumber; n++)
			{
				for (j = 0; j < ProteinIter->m_vPeptidesAttributes.at(n).m_vecAttributes.size(); j++)
				{
					m_vTrainx.push_back(ProteinIter->m_vPeptidesAttributes.at(n).m_vecAttributes.at(j));
				}
				dytemp = ProteinIter->m_vAdjustPeptidesIntensityFactor.at(n);
				if (dytemp < m_dMiny) m_dMiny = dytemp;
				if (dytemp > m_dMaxy) m_dMaxy = dytemp;
				m_cAllys.sy += dytemp; // sum of y
				m_cAllys.sy2 += dytemp*dytemp; // sum of y^2
				m_vTrainy.push_back(dytemp);
				m_iTrainSampleNumber++;
			}

		} // end if PeptidesThreshold

	} // end for proteins

	if (m_iTrainSampleNumber < 1) 
	{
		cout << "error FeatureNumber<1\n";
		log << "Error:\tFeatureNumber<1 when constructing the train set\n";
		exit(0);
		return 0;
	}

	m_cAllys.n = m_iTrainSampleNumber;
	cout << "y read in:\n";
	cout << "n: " << m_iTrainSampleNumber << endl;
	cout << "y first and last:\n";
	cout << m_vTrainy[0] << ", " << m_vTrainy[m_iTrainSampleNumber - 1] << endl;
	log << "y read in:\n";
	log << "n: " << m_iTrainSampleNumber << endl;
	log << "y first and last:\t";
	log << m_vTrainy[0] << ", " << m_vTrainy[m_iTrainSampleNumber - 1] << endl;
	m_dMeany = m_cAllys.sy / m_iTrainSampleNumber;
	m_dStdy = sqrt((m_cAllys.sy2 - m_iTrainSampleNumber*m_dMeany*m_dMeany) / (m_iTrainSampleNumber - 1));
	cout << "ybar,shat: " << m_dMeany << ", " << m_dStdy << endl;
	cout << "Constructed XY using " << m_iTrainSampleNumber << " samples in " << m_iNumberOfTrainProteins<<" proteins." << endl; // added by CC 20150527
	log << "ybar,shat: " << m_dMeany << ", " << m_dStdy << endl;
	log << "Constructed XY using " << m_iTrainSampleNumber << " samples in " << m_iNumberOfTrainProteins << " proteins." << endl; // added by CC 20150527

	m_iFeatureNumber = m_vTrainx.size() / m_iTrainSampleNumber;

	if (m_vTrainx.size() != m_iTrainSampleNumber*m_iFeatureNumber)
	{
		cout << "error: input x file has wrong number of values\n";
		log << "Error:\tInput x file has wrong number of values when constructing the train set\n";
		exit(0);
	}


	cout << "x read in:\n";
	cout << "p: " << m_iFeatureNumber << endl;
	cout << "first row: " << m_vTrainx[0] << " ...  " << m_vTrainx[m_iFeatureNumber - 1] << endl;
	cout << "last row: " << m_vTrainx[(m_iTrainSampleNumber - 1)*m_iFeatureNumber] << " ...  " << m_vTrainx[m_iTrainSampleNumber*m_iFeatureNumber - 1] << endl;
	log << "x read in:\n";
	log << "p: " << m_iFeatureNumber << endl;
	log << "first row: " << m_vTrainx[0] << " ...  " << m_vTrainx[m_iFeatureNumber - 1] << endl;
	log << "last row: " << m_vTrainx[(m_iTrainSampleNumber - 1)*m_iFeatureNumber] << " ...  " << m_vTrainx[m_iTrainSampleNumber*m_iFeatureNumber - 1] << endl;

	 return true;
}

bool CBART::mf_ConstructTestX(std::ofstream& log,vector<CProtein> proteins)
{
	//data for predictions
	cout << "\tConstructing Test XY\n";
	log << "\tConstructing Test XY\n";
	vector<CProtein>::iterator  ProteinIter;
	int j = 0; //feature num
	int n = 0; //peptide num per protein
	double dytemp;

	m_iTestSampleNumber = 0;
	m_iNumberOfTestProteins = 0;
	m_vTestx.clear();
	m_vTesty.clear();

	for (ProteinIter = proteins.begin(); ProteinIter != proteins.end(); ProteinIter++)
	{
		if (ProteinIter->m_iPeptidesNumber >= UniquePeptidesTestThreshold)
		{
			m_iNumberOfTestProteins++;
			for (n = 0; n < ProteinIter->m_iPeptidesNumber; n++)
			{
				for (j = 0; j < ProteinIter->m_vPeptidesAttributes.at(n).m_vecAttributes.size(); j++)
				{
					m_vTestx.push_back(ProteinIter->m_vPeptidesAttributes.at(n).m_vecAttributes.at(j));
				}
				dytemp = ProteinIter->m_vAdjustPeptidesIntensityFactor.at(n);
				m_vTesty.push_back(dytemp);
				m_iTestSampleNumber++;
			}

		} // end if PeptidesThreshold

	} // end for proteins

	dip.n = 0;
	if (m_vTestx.size() != m_iFeatureNumber*m_iTestSampleNumber)
	{
		cout << "error, wrong number of elements in prediction data set\n";
		log << "Error:\tWrong number of elements in prediction data set\n";
		exit(0);
		return 0;
	}
	if (m_iTestSampleNumber)
	{
		dip.n = m_iTestSampleNumber; dip.p = m_iFeatureNumber; dip.x = &m_vTestx[0]; dip.y = &m_vTesty[0]; //there are no y's!
		cout << "Constructed TestXY using " << m_iTestSampleNumber << " peptides in " << m_iNumberOfTestProteins << " proteins." << endl;
		log << "Constructed TestXY using " << m_iTestSampleNumber << " peptides in " << m_iNumberOfTestProteins << " proteins." << endl;

	}
	cout << "x for prediction read in:\n";
	cout << "n: " << dip.n << endl;
	log << "x for prediction read in:\n";
	log << "n: " << dip.n << endl;

	if (dip.n) {
		cout << "first row: " << dip.x[0] << " ...  " << dip.x[m_iFeatureNumber - 1] << endl;
		cout << "last row: " << dip.x[(dip.n - 1)*m_iFeatureNumber] << " ...  " << dip.x[dip.n*m_iFeatureNumber - 1] << endl;
		log << "first row: " << dip.x[0] << " ...  " << dip.x[m_iFeatureNumber - 1] << endl;
		log << "last row: " << dip.x[(dip.n - 1)*m_iFeatureNumber] << " ...  " << dip.x[dip.n*m_iFeatureNumber - 1] << endl;

	}
}

void CBART::mf_SelectFeature(std::ofstream& log, vector<tree> trees, map<size_t, size_t> &FeatureCount)
{
	tree::cnpv npv;
	size_t inodeIndex = 0;
	map<size_t, size_t>::iterator FeatureCountIter;
	map<size_t, size_t>::iterator mapFeatureCountIter;
	vector<size_t> SortedFeatureCount;
	vector<size_t> SortedFeatureIndex;
	
	ofstream SelectedFeatureFile(m_strSelectedFeaturePath);
	if (!SelectedFeatureFile.is_open())
	{
		cout << "can not open SelectedFeature.txt" << endl;
		log << "can not open SelectedFeature.txt" << endl;
	}
	for (int i = 0; i < trees.size(); i++)
	{
		npv.clear();
		trees[i].getnodes(npv);
		for (inodeIndex = 0; inodeIndex < npv.size(); inodeIndex++)
		{
			if (npv[inodeIndex]->ntype() != 'b')
			{
				FeatureCountIter = FeatureCount.find(npv[inodeIndex]->getv());
				if (FeatureCountIter != FeatureCount.end())
				{
					FeatureCountIter->second++;
				}
				else
				{
					FeatureCount.insert(pair<size_t, size_t>(npv[inodeIndex]->getv(), 1));
				}
			}
		}
	}
	
	mapFeatureCountIter = FeatureCount.begin();
	for (; mapFeatureCountIter != FeatureCount.end(); mapFeatureCountIter++)
	{
		SelectedFeatureFile << mapFeatureCountIter->first << "\t" << mapFeatureCountIter->second << endl;
		SortedFeatureIndex.push_back(mapFeatureCountIter->first);
		SortedFeatureCount.push_back(mapFeatureCountIter->second);
	}

	quickSort(SortedFeatureCount, SortedFeatureIndex, 0, SortedFeatureIndex.size() - 1);
	for (size_t i = 0; i < SortedFeatureIndex.size(); i++)
	{
		SelectedFeatureFile << SortedFeatureIndex.at(i) << "\t" << SortedFeatureCount.at(i) << endl;
	}
	SelectedFeatureFile.close();
}

bool CBART::mf_BartRegressionRun(std::ofstream& log, vector<CProtein>&proteins, size_t burn, size_t nd, size_t m, double lambda, double nu, double kfac, double alpha, double beta)
{

	cout << "Bayes addictive regression trees\n";
	cout << "burn,nd,number of trees: " << burn << ", " << nd << ", " << m << endl;
	cout << "lambda,nu,kfac: " << lambda << ", " << nu << ", " << kfac << endl;
	
	log << "Bayes addictive regression trees\n";
	log << "burn,nd,number of trees: " << burn << ", " << nd << ", " << m << endl;
	log << "lambda,nu,kfac: " << lambda << ", " << nu << ", " << kfac << endl;

	//string strCheckTestX = "TestXY.txt";
	//mf_saveTestXY(strCheckTestX);


	double* pAllfit;//sum of fit of all trees
	double* pResidual; //y-(pAllfit-ftemp) = y-pAllfit+ftemp
	double* pFittemp; //fit of current tree
	//in sample fit
	double* pPredictTrainMean; //posterior mean of in-sample fit, sum draws,then divide
	//out of sample fit
	double* pPredictTesttemp; //temporary fit vector to compute prediction
	double* pPredictTestMean; //posterior mean for prediction

	//x cutpoints
	size_t nc = 100; //100 equally spaced cutpoints from min to max.
	makexinfo(m_iFeatureNumber, m_iTrainSampleNumber, &m_vTrainx[0], xi, nc);

	//trees
	tree cTreeTemp;
	cTreeTemp.setm(m_dMeany / m);
	for (int i = 0; i < m; i++)
	{
		m_vTrees.push_back(cTreeTemp);
	}

	//dinfo
	pAllfit = new double[m_iTrainSampleNumber]; //sum of fit of all trees
	for (size_t i = 0; i<m_iTrainSampleNumber; i++) pAllfit[i] = m_dMeany;
	pResidual = new double[m_iTrainSampleNumber]; //y-(pAllfit-ftemp) = y-pAllfit+ftemp
	pFittemp = new double[m_iTrainSampleNumber]; //fit of current tree
	di.n = m_iTrainSampleNumber; di.p = m_iFeatureNumber; di.x = &m_vTrainx[0]; di.y = pResidual; //the y for each draw will be the residual 

	//pi.pbd = 1.0; //prob of birth/death move
	pi.pb = .5; //prob of birth given  birth/death

	pi.alpha = alpha; //prior prob a bot node splits is alpha/(1+d)^beta, d is depth of node
	pi.beta = beta; //2 for bart means it is harder to build big trees.
	pi.tau = (m_dMaxy - m_dMiny) / (2 * kfac*sqrt((double)m));
	pi.sigma = m_dStdy;
	cout << "alpha, beta: " << pi.alpha << ", " << pi.beta << endl;
	cout << "sigma, tau: " << pi.sigma << ", " << pi.tau << endl;
	log << "alpha, beta: " << pi.alpha << ", " << pi.beta << endl;
	log << "sigma, tau: " << pi.sigma << ", " << pi.tau << endl;
	
	//--------------------------------------------------
	//storage for ouput
	//in sample fit
	pPredictTrainMean = new double[m_iTrainSampleNumber]; //posterior mean of in-sample fit, sum draws,then divide
	for (size_t i = 0; i<m_iTrainSampleNumber; i++) pPredictTrainMean[i] = 0.0;

	//out of sample fit
	if (dip.n) {
		pPredictTestMean = new double[dip.n];
		pPredictTesttemp = new double[dip.n];
		for (size_t i = 0; i<dip.n; i++) pPredictTestMean[i] = 0.0;
	}

	//for sigma draw
	double rss;  //residual sum of squares
	double restemp; //a residual
	//std::ofstream bsdf("bart-sd.txt"); //note that we write all burn+nd draws to this file.

	//std::ofstream btsf("bart-ave-tree-size.txt"); //note that we write all burn+nd draws to this file.
	double ats; //place for average tree size
	
	double anb; //place for average number of bottom nodes
//	std::ofstream bnbf("bart-ave-num-bots.txt"); //note that we write all burn+nd draws to this file.

	//--------------------------------------------------
	//mcmc
	//random number generation
	uint seed = 99;
	RNG gen(seed); //this one random number generator is used in all draws

	cout << "MCMC:\n";
	log << "MCMC:\n";
	clock_t tp;
	tp = clock();
	for (size_t i = 0; i<(nd + burn); i++) 
	{
		if (i % 100 == 0)
		{
			cout << "i: " << i << endl;
			log << "i: " << i << endl;
		}
			
		//draw trees
		for (size_t j = 0; j<m; j++) 
		{
			fit(m_vTrees[j], xi, di, pFittemp);
			for (size_t k = 0; k<m_iTrainSampleNumber; k++) 
			{
				pAllfit[k] = pAllfit[k] - pFittemp[k];
				pResidual[k] = m_vTrainy[k] - pAllfit[k];
			}
			bd(log,m_vTrees[j], xi, di, pi, gen);
			drmu(m_vTrees[j], xi, di, pi, gen);
			fit(m_vTrees[j], xi, di, pFittemp);
			for (size_t k = 0; k<m_iTrainSampleNumber; k++) pAllfit[k] += pFittemp[k];
		}
		//draw sigma
		rss = 0.0;
		for (size_t k = 0; k<m_iTrainSampleNumber; k++) 
		{ 
			restemp = m_vTrainy[k] - pAllfit[k];
			rss += restemp*restemp;
		}
		pi.sigma = sqrt((nu*lambda + rss) / gen.chi_square(nu + m_iTrainSampleNumber));
		//bsdf << pi.sigma << endl;
		ats = 0.0; anb = 0.0;
		for (size_t k = 0; k<m; k++) 
		{
			ats += m_vTrees[k].treesize();
			anb += m_vTrees[k].nbots();
		}
		//btsf << ats / m << endl;
		//bnbf << anb / m << endl;
		if (i >= burn) 
		{
			for (size_t k = 0; k<m_iTrainSampleNumber; k++) 
				pPredictTrainMean[k] += pAllfit[k];
			if (dip.n) 
			{
				for (size_t j = 0; j<m; j++) 
				{
					fit(m_vTrees[j], xi, dip, pPredictTesttemp);
					for (size_t k = 0; k<dip.n; k++) pPredictTestMean[k] += pPredictTesttemp[k];
				}
			}
		}
	}
	tp = clock() - tp;
	double thetime = (double)(tp) / (double)(CLOCKS_PER_SEC);
	cout << "time for loop: " << thetime << endl;
	log << "time for loop: " << thetime << endl;
	//std::ofstream timef("time.txt");
	//timef << thetime << endl;

	m_vTrainPredicty.clear();
	for (size_t k = 0; k < m_iTrainSampleNumber; k++)
	{
		pPredictTrainMean[k] /= (double)nd;
		m_vTrainPredicty.push_back(pPredictTrainMean[k]);
	}
	m_dAutoCorrelation = spearsonCorrelation(m_vTrainPredicty, m_vTrainy);
	//std::ofstream bfitf("bart-fit.txt");
	//for (size_t i = 0; i<m_iTrainSampleNumber; i++) 
	//	bfitf << pPredictTrainMean[i] << endl;

	if (dip.n) {
		for (size_t k = 0; k < dip.n; k++)
		{
			pPredictTestMean[k] /= (double)nd;
			m_vTestPredicty.push_back(pPredictTestMean[k]);
		}
		//std::ofstream bpredfitf("bart-pred-fit.txt");
		//for (size_t i = 0; i<dip.n; i++) bpredfitf << pPredictTestMean[i] << endl;
	}

	//std::ofstream treef("trees.txt");
	//treef << xi << endl; //the cutpoints
	//mf_xitofile(treef, xi);
	//treef << m << endl;  //number of trees
	//treef << m_iFeatureNumber << endl;  //dimension of x's
	//for (size_t j = 0; j<m; j++) treef << m_vTrees[j] << endl;  //all the trees

	map<size_t, size_t> mapFeatureCount;
	map<size_t, size_t>::iterator mapFeatureCountIter;
	vector<size_t> SortedFeatureCount;
	vector<size_t> SortedFeatureIndex;

	mf_SelectFeature(log, m_vTrees, mapFeatureCount);

	mapFeatureCountIter = mapFeatureCount.begin();
	for (; mapFeatureCountIter != mapFeatureCount.end(); mapFeatureCountIter++)
	{
		//treef << mapFeatureCountIter->first << "\t" << mapFeatureCountIter->second << endl;
		SortedFeatureIndex.push_back(mapFeatureCountIter->first);
		SortedFeatureCount.push_back(mapFeatureCountIter->second);
	}

	quickSort(SortedFeatureCount, SortedFeatureIndex, 0, SortedFeatureIndex.size() - 1);

	Discretization(pPredictTestMean, m_iTestSampleNumber);
	//Assertion(pPredictTestMean, m_iTestSampleNumber);


	vector<CProtein>::iterator proteinIter;
	vector<double>::iterator peptidesIter;
	int index = 0;
	for (proteinIter = proteins.begin(); proteinIter != proteins.end(); proteinIter++)
	{
		if (proteinIter->m_iPeptidesNumber >= UniquePeptidesTestThreshold)
		{
			proteinIter->m_vPeptidesCorrectionFactor.clear();
			for (peptidesIter = proteinIter->m_vPeptidesNativeIntensity.begin(); peptidesIter != proteinIter->m_vPeptidesNativeIntensity.end(); peptidesIter++) //use vPeptidesNativeIntensity instead of vPeptidesIntensity! by CC 20150527
			{
				proteinIter->m_vPeptidesCorrectionFactor.push_back(pPredictTestMean[index]);
				index++;
			}
		}
		

	}

	delete pAllfit;
	pAllfit = NULL;
	delete pResidual;
	pResidual = NULL;
	delete pFittemp;
	pFittemp = NULL;
	delete pPredictTrainMean;
	pPredictTrainMean = NULL;
	delete pPredictTestMean;
	pPredictTestMean = NULL;
	delete pPredictTesttemp;
	pPredictTesttemp = NULL;

	return 0;
}


void CBART::mf_MakeIndexPartition(size_t SampleNumber, vector<size_t> &IndexPartation)
{
	for (size_t i = 0; i < SampleNumber; i++)
	{
		IndexPartation.push_back(floor(rand()%10));
	}
}

void CBART::mf_ConstructCVData(std::ofstream& log, vector<CProtein> proteins, vector<size_t>IndexPartation, size_t fold, vector<double> &TrainxForCV, vector<double> &TrainYForCV, vector<double>&TestxForCV, vector<double>&TestyForCV)
{
	m_dMiny = INFINITY; //use range of y to calibrate prior for bottom node mu's
	m_dMaxy = -INFINITY;

	//cout << "\tConstructing CV XY\n\n";
	vector<CProtein>::iterator  ProteinIter;
	int j = 0; //feature num
	int n = 0; //peptide num per protein
	double dytemp;

	TrainxForCV.clear();
	TrainYForCV.clear();
	TestxForCV.clear();
	m_iTrainSampleNumber = 0;  //Change i by m_NumberOfLimitedSamples gzhq 20150528
	m_iTestSampleNumber = 0;
	m_cAllys.clear();
	m_iNumberOfTestProteins = 0;
	m_iNumberOfTrainProteins = 0;
	int i = 0;
	for (ProteinIter = proteins.begin(); ProteinIter != proteins.end(); ProteinIter++)
	{		
		if (ProteinIter->m_iPeptidesNumber >= UniquePeptidesTrainThreshold)
		{
			if (IndexPartation.at(i) == fold)
			{
				for (n = 0; n < ProteinIter->m_iPeptidesNumber; n++)
				{
					for (j = 0; j < ProteinIter->m_vPeptidesAttributes.at(n).m_vecAttributes.size(); j++)
					{
						TestxForCV.push_back(ProteinIter->m_vPeptidesAttributes.at(n).m_vecAttributes.at(j));
					}
					dytemp = ProteinIter->m_vAdjustPeptidesIntensityFactor.at(n);
					TestyForCV.push_back(dytemp);
					m_iTestSampleNumber++;
				}
				m_iNumberOfTestProteins++;
			}
			else
			{
				for (n = 0; n < ProteinIter->m_iPeptidesNumber; n++)
				{
					for (j = 0; j < ProteinIter->m_vPeptidesAttributes.at(n).m_vecAttributes.size(); j++)
					{
						TrainxForCV.push_back(ProteinIter->m_vPeptidesAttributes.at(n).m_vecAttributes.at(j));
					}
					dytemp = ProteinIter->m_vAdjustPeptidesIntensityFactor.at(n);
					if (dytemp < m_dMiny) m_dMiny = dytemp;
					if (dytemp > m_dMaxy) m_dMaxy = dytemp;
					m_cAllys.sy += dytemp; // sum of y
					m_cAllys.sy2 += dytemp*dytemp; // sum of y^2
					TrainYForCV.push_back(dytemp);
					m_iTrainSampleNumber++;
				}
				m_iNumberOfTrainProteins++;
			}
			i++;

		} // end if PeptidesThreshold

	} // end for proteins

	m_iTrainSampleNumber = TrainYForCV.size();

	m_cAllys.n = m_iTrainSampleNumber;
	m_dMeany = m_cAllys.sy / m_iTrainSampleNumber;
	m_dStdy = sqrt((m_cAllys.sy2 - m_iTrainSampleNumber*m_dMeany*m_dMeany) / (m_iTrainSampleNumber - 1));

	m_iFeatureNumber = TrainxForCV.size() / m_iTrainSampleNumber;

	if (TrainxForCV.size() != m_iTrainSampleNumber*m_iFeatureNumber) 
	{
		cout << "error: input x file has wrong number of values\n";
		log << "error: input x file has wrong number of values\n";
	}

	if (m_iTestSampleNumber)
	{
		dip.n = m_iTestSampleNumber; dip.p = m_iFeatureNumber; dip.x = &TestxForCV[0]; dip.y = &TestyForCV[0]; //there are no y's!
		cout << "Constructed TestXY using " << m_iTestSampleNumber << " peptides in " << m_iNumberOfTestProteins << " proteins." << endl;
		log << "Constructed TestXY using " << m_iTestSampleNumber << " peptides in " << m_iNumberOfTestProteins << " proteins." << endl;
	}

}
void CBART::mf_TrainAndPredictForCV(std::ofstream& log, vector<CProtein> &proteins, vector<double> &TrainxForCV, vector<double> &TrainYForCV, vector<double>&TestxForCV, vector<double> &PredictedY, size_t burn, size_t nd, size_t m, double lambda, double nu, double kfac, double alpha, double beta)
{

	double* pAllfit;//sum of fit of all trees
	double* pResidual; //y-(pAllfit-ftemp) = y-pAllfit+ftemp
	double* pFittemp; //fit of current tree
	//in sample fit
	double* pPredictTrainMean; //posterior mean of in-sample fit, sum draws,then divide
	//out of sample fit
	double* pPredictTesttemp; //temporary fit vector to compute prediction
	double* pPredictTestMean; //posterior mean for prediction

	//x cutpoints
	size_t nc = 100; //100 equally spaced cutpoints from min to max.
	makexinfo(m_iFeatureNumber, m_iTrainSampleNumber, &TrainxForCV[0], xi, nc);

	//trees
	tree cTreeTemp;
	cTreeTemp.setm(m_dMeany / m);
	m_vTrees.clear();
	for (int i = 0; i < m; i++)
	{
		m_vTrees.push_back(cTreeTemp);
	}

	//dinfo
	pAllfit = new double[m_iTrainSampleNumber]; //sum of fit of all trees
	for (size_t i = 0; i<m_iTrainSampleNumber; i++) pAllfit[i] = m_dMeany;
	pResidual = new double[m_iTrainSampleNumber]; //y-(pAllfit-ftemp) = y-pAllfit+ftemp
	pFittemp = new double[m_iTrainSampleNumber]; //fit of current tree
	di.n = m_iTrainSampleNumber; di.p = m_iFeatureNumber; di.x = &m_vTrainx[0]; di.y = pResidual; //the y for each draw will be the residual 

	//pi.pbd = 1.0; //prob of birth/death move
	pi.pb = .5; //prob of birth given  birth/death

	pi.alpha =alpha; //prior prob a bot node splits is alpha/(1+d)^beta, d is depth of node
	pi.beta = beta; //2 for bart means it is harder to build big trees.

	pi.tau = (m_dMaxy - m_dMiny) / (2 * kfac*sqrt((double)m));
	pi.sigma = m_dStdy;
	//cout << "\nalpha, beta: " << pi.alpha << ", " << pi.beta << endl;
	//cout << "sigma, tau: " << pi.sigma << ", " << pi.tau << endl;

	//--------------------------------------------------

	//--------------------------------------------------
	//storage for ouput
	//in sample fit
	pPredictTrainMean = new double[m_iTrainSampleNumber]; //posterior mean of in-sample fit, sum draws,then divide
	for (size_t i = 0; i<m_iTrainSampleNumber; i++) pPredictTrainMean[i] = 0.0;

	//out of sample fit
	if (dip.n) {
		pPredictTestMean = new double[dip.n];
		pPredictTesttemp = new double[dip.n];
		for (size_t i = 0; i<dip.n; i++) pPredictTestMean[i] = 0.0;
	}

	//for sigma draw
	double rss;  //residual sum of squares
	double restemp; //a residual
	//std::ofstream bsdf("bart-sd.txt"); //note that we write all burn+nd draws to this file.

	//std::ofstream btsf("bart-ave-tree-size.txt"); //note that we write all burn+nd draws to this file.
	//double ats; //place for average tree size
	//double anb; //place for average number of bottom nodes
	//std::ofstream bnbf("bart-ave-num-bots.txt"); //note that we write all burn+nd draws to this file.

	//--------------------------------------------------
	//mcmc

	//random number generation
	uint seed = 99;
	RNG gen(seed); //this one random number generator is used in all draws

	cout << "\nMCMC:\n";
	log << "\nMCMC:\n";
	clock_t tp;
	tp = clock();
	for (size_t i = 0; i<(nd + burn); i++)
	{
		//if (i % 100 == 0) cout << "i: " << i << endl;
		//draw trees
		for (size_t j = 0; j<m; j++)
		{
			fit(m_vTrees[j], xi, di, pFittemp);
			for (size_t k = 0; k<m_iTrainSampleNumber; k++)
			{
				pAllfit[k] = pAllfit[k] - pFittemp[k];
				pResidual[k] = TrainYForCV[k] - pAllfit[k];
			}
			bd(log,m_vTrees[j], xi, di, pi, gen);
			drmu(m_vTrees[j], xi, di, pi, gen);
			fit(m_vTrees[j], xi, di, pFittemp);
			for (size_t k = 0; k<m_iTrainSampleNumber; k++) pAllfit[k] += pFittemp[k];
		}
		//draw sigma
		rss = 0.0;
		for (size_t k = 0; k<m_iTrainSampleNumber; k++)
		{
			restemp = TrainYForCV[k] - pAllfit[k];
			rss += restemp*restemp;
		}
		pi.sigma = sqrt((nu*lambda + rss) / gen.chi_square(nu + m_iTrainSampleNumber));
		//bsdf << pi.sigma << endl;
		//ats = 0.0; anb = 0.0;
		//for (size_t k = 0; k<m; k++)
		//{
		//	ats += m_vTrees[k].treesize();
		//	anb += m_vTrees[k].nbots();
		//}
		//btsf << ats / m << endl;
		//bnbf << anb / m << endl;
		if (i >= burn)
		{
			for (size_t k = 0; k<m_iTrainSampleNumber; k++)
				pPredictTrainMean[k] += pAllfit[k];
			if (dip.n)
			{
				for (size_t j = 0; j<m; j++)
				{
					fit(m_vTrees[j], xi, dip, pPredictTesttemp);
					for (size_t k = 0; k<dip.n; k++) pPredictTestMean[k] += pPredictTesttemp[k];
				}
			}
		}
	}
	tp = clock() - tp;
	double thetime = (double)(tp) / (double)(CLOCKS_PER_SEC);
	//cout << "time for loop: " << thetime << endl;
	//std::ofstream timef("time.txt");
	//timef << thetime << endl;

 	for (size_t k = 0; k<m_iTrainSampleNumber; k++) pPredictTrainMean[k] /= (double)nd;
	//std::ofstream bfitf("bart-fit.txt");
	//for (size_t i = 0; i<m_iTrainSampleNumber; i++) bfitf << pPredictTrainMean[i] << endl;

	if (dip.n) 
	{
		for (size_t k = 0; k < dip.n; k++)
		{
			pPredictTestMean[k] /= (double)nd;
			PredictedY.push_back(pPredictTestMean[k]);
		}
		
		//std::ofstream bpredfitf("bart-pred-fit.txt");
		//for (size_t i = 0; i<dip.n; i++) bpredfitf << pPredictTestMean[i] << endl;
	}
	
	delete pAllfit;
	pAllfit = NULL;
	delete pResidual;
	pResidual = NULL;
	delete pFittemp;
	pFittemp = NULL;
	delete pPredictTrainMean;
	pPredictTrainMean = NULL;
	delete pPredictTestMean;
	pPredictTestMean = NULL;
	delete pPredictTesttemp;
	pPredictTesttemp = NULL;
	
}
void CBART::mf_AdjustPeptidesIntensity(vector<CProtein>& proteins)
{
	vector<CProtein>::iterator proteinIter;
	vector<double>::iterator peptidesIter;
	double dPeptideIntensityTemp;
	int i = 0;
	proteinIter = proteins.begin();
	for (; proteinIter != proteins.end(); proteinIter++)
	{
		i = 0;
		if (proteinIter->m_iPeptidesNumber >= UniquePeptidesTrainThreshold)
		{
			proteinIter->m_bIfCalculateSPEAQInPeptides = true;
			if (proteinIter->m_iPeptidesNumber > UniquePeptidesCorrectThreshold)
			{			
				for (peptidesIter = proteinIter->m_vPeptidesAdjustIntensity.begin(); peptidesIter != proteinIter->m_vPeptidesAdjustIntensity.end(); peptidesIter++) //use vPeptidesNativeIntensity instead of vPeptidesIntensity! by CC 20150527
				{
					dPeptideIntensityTemp = *peptidesIter;
					if (proteinIter->m_vPeptidesCorrectionFactor.at(i) != 0.0)
						*peptidesIter = dPeptideIntensityTemp / proteinIter->m_vPeptidesCorrectionFactor.at(i);
					else
						*peptidesIter = 0.0;
					dPeptideIntensityTemp = *peptidesIter;
					i++;
					proteinIter->m_dProteinSPEAQInPeptides += dPeptideIntensityTemp;
				}
			}
			else
			{ //not enough peptides for correction
				for (peptidesIter = proteinIter->m_vPeptidesAdjustIntensity.begin(); peptidesIter != proteinIter->m_vPeptidesAdjustIntensity.end(); peptidesIter++) //use vPeptidesNativeIntensity instead of vPeptidesIntensity! by CC 20150527
				{
					//cout << index << endl;			
					dPeptideIntensityTemp = *peptidesIter;
					proteinIter->m_dProteinSPEAQInPeptides += dPeptideIntensityTemp;
				}
			}
			if (proteinIter->m_iPeptidesNumber > 0)
			{
				proteinIter->m_dProteinSPEAQInPeptides = proteinIter->m_dProteinSPEAQInPeptides / proteinIter->m_iPeptidesNumber;
			}

			proteinIter->m_dPeptidesIntensityMax = GetMedianIntensity(proteinIter->m_vPeptidesAdjustIntensity);
			if (proteinIter->m_dPeptidesIntensityMax != 0.0)
			{
				for (size_t i = 0; i < proteinIter->m_vAdjustPeptidesIntensityFactor.size(); i++)
					proteinIter->m_vAdjustPeptidesIntensityFactor.at(i) = proteinIter->m_vPeptidesAdjustIntensity.at(i) / proteinIter->m_dPeptidesIntensityMax;
			}
		}
	}

}
void CBART::mf_ChooseBestParameter(std::ofstream& log, vector<CProtein> &proteins, CTrainParam trainparam)
{// now we need to select the parameters: lambda nu pi.alpha pi.beta  m 
	//using ten-fold crossvalidation for parameters selection

	size_t burn = 100;
	size_t nd = 100;
	size_t m = 200;
	double lambda = 1;
	double nu = 3;
	double kfac = 2;
	double alpha = 0.95;
	double beta = 1;
	vector<size_t> vecIndexPartation;
	

	ofstream ParameterSelectByCVResult("ParameterSelectByCVResult.txt");
	ParameterSelectByCVResult << "lambda\tkfac\tnu\tm\t";
	ParameterSelectByCVResult << "alpha\tbeta\tTestCorrelation";
	ParameterSelectByCVResult<<"\tSpearsonWithReal\n";

	vector<double> vecTrainxForCV;
	vector<double> vecTrainYForCV;
	vector<double> vecTestXForCV;
	vector<double> vecTestYForCV;
	vector<double> vecPredictedY;
	int iFold = 0;
	vector<double> vecSpearsonCorrelation;
	vector<double> vecUPS2SpearsonCorr;
	double dBestSpearsonCorrelation = 0.0;
	double dBestUPS2SpearsonCorr=0.0;
	double dSpearsonCorrelationTemp;
	double dUPS2PearsonCorr;
	double dbestLamda;
	double dBestNu;
	double dBestalpha;
	double dBestBeta;
	double dBestm;
	
	char Buffer[BUFFERLENGTH];
	char *pstr;
	char *pstr1;
	string strTemp;
	FILE * pFile;

	map<string, double> mapUPS2Id2Mols;
	map<string, double> mapUPS2Id2MW;
	map<string, double> mapUPS2iBAQs;
	map<string, double>::iterator Ups2MolsIter;
	map<string, double>::iterator Ups2WMIter;
	map<string, double>::iterator UPS2IBAQIter;

	double meanReal = 0.0;
	double meanTest = 0.0;
	double sd1 = 0.0, sd2 = 0.0;
	vector<string> vecUPS2IDIdentified;
	vector<double> vecUPS2logPredictIdentified;
	vector<double> vecUPS2logMolsIdentified;

	int iNumberofUPS2Identified;
	string::size_type stPosition;
	double dMolsTemp;
	double dWMTemp;

	vector<CProtein>::iterator proteinIter;
	vector<double>::iterator peptidesIter;
	vector<double>::iterator PeptidesCorrectionFactorIter;

	vector<double> AdjustPeptidesIntensity;
	vector<double> NativePeptidesIntensity;
	double dPeptideIntensityTemp;
	int index = 0;
	int i = 0;

	if (trainparam.m_bIfCotainStardProtein)
	{
		CProteinWorker proteinworker;
		//string UPS2path = "UPS2_mols.ini";
		proteinworker.mf_LoadUPS2Mols(log, trainparam.m_strAbundanceOfStandProteinsPath, mapUPS2Id2Mols, mapUPS2Id2MW);
	}

	ofstream pOut("TestCorrelation.txt");

	mf_MakeIndexPartition(m_iNumberOfTrainProteins, vecIndexPartation);
	mf_ConstructXY(log, proteins);
	mf_ConstructTestX(log, proteins);

	for (alpha = 0.8; alpha < 1; alpha = alpha + 0.05)
		for (beta = 0.1; beta <= 2; beta = beta + 0.2)
		{
			iFold = 0;
			vecTestYForCV.clear();
			vecPredictedY.clear();
			//for 10-fold CrossValidation
			while (iFold <= 9)
			{
				mf_ConstructCVData(log,proteins, vecIndexPartation, iFold, vecTrainxForCV, vecTrainYForCV, vecTestXForCV, vecTestYForCV);
				mf_TrainAndPredictForCV(log,proteins, vecTrainxForCV, vecTrainYForCV, vecTestXForCV, vecPredictedY, burn, nd, m, lambda, nu, kfac, alpha, beta);
				iFold++;
			}
			dSpearsonCorrelationTemp = spearsonCorrelation(vecPredictedY, vecTestYForCV);
			mf_ConstructXY(log,proteins);
			mf_ConstructTestX(log, proteins);
			vecPredictedY.clear();
			mf_TrainAndPredictForCV(log, proteins, m_vTrainx, m_vTrainy, m_vTestx, vecPredictedY, burn, nd, m, lambda, nu, kfac, alpha, beta);

			Discretization(vecPredictedY);
			//Assertion(pPredictTestMean, m_iTestSampleNumber);
			for (proteinIter = proteins.begin(); proteinIter != proteins.end(); proteinIter++)
			{
				if (proteinIter->m_iPeptidesNumber >= UniquePeptidesTestThreshold)
				{
					proteinIter->m_vPeptidesCorrectionFactor.clear();
					proteinIter->m_dProteinSPEAQInPeptides = 0.0;
					proteinIter->m_dProteinSPEAQInProtein = 0.0;
					proteinIter->m_bIfCalculateSPEAQInPeptides = false;
					proteinIter->m_vPeptidesAdjustIntensity.clear();

				}
			}

			index = 0;
			for (proteinIter = proteins.begin(); proteinIter != proteins.end(); proteinIter++)
			{
				if (proteinIter->m_iPeptidesNumber >= UniquePeptidesTestThreshold)
				{
					proteinIter->m_vPeptidesCorrectionFactor.clear();
					for (peptidesIter = proteinIter->m_vPeptidesNativeIntensity.begin(); peptidesIter != proteinIter->m_vPeptidesNativeIntensity.end(); peptidesIter++) //use vPeptidesNativeIntensity instead of vPeptidesIntensity! by CC 20150527
					{
						proteinIter->m_vPeptidesCorrectionFactor.push_back(vecPredictedY[index]);
						index++;
					}
					proteinIter->m_dProteinSPEAQInPeptides = 0.0;
					proteinIter->m_dProteinSPEAQInProtein = 0.0;
				}

			}

			i = 0;
			for (proteinIter = proteins.begin(); proteinIter != proteins.end(); proteinIter++)
			{
				i = 0;
				if (proteinIter->m_iPeptidesNumber >= UniquePeptidesTestThreshold)
				{
					proteinIter->m_bIfCalculateSPEAQInPeptides = true;
					if (proteinIter->m_iPeptidesNumber > UniquePeptidesCorrectThreshold)
					{
						for (peptidesIter = proteinIter->m_vPeptidesNativeIntensity.begin(); peptidesIter != proteinIter->m_vPeptidesNativeIntensity.end(); peptidesIter++) //use vPeptidesNativeIntensity instead of vPeptidesIntensity! by CC 20150527
						{
							//cout << index << endl;			
							dPeptideIntensityTemp = *peptidesIter;
							proteinIter->m_dProteinSPEAQInProtein += dPeptideIntensityTemp;
							if (proteinIter->m_vPeptidesCorrectionFactor.at(i) != 0.0)
								dPeptideIntensityTemp = dPeptideIntensityTemp / proteinIter->m_vPeptidesCorrectionFactor.at(i);
							else
								dPeptideIntensityTemp = 0.0;
							i++;
							proteinIter->m_vPeptidesAdjustIntensity.push_back(dPeptideIntensityTemp);
							proteinIter->m_dProteinSPEAQInPeptides += dPeptideIntensityTemp;
						}
					}
					else
					{ //not enough peptides for correction
						for (peptidesIter = proteinIter->m_vPeptidesNativeIntensity.begin(); peptidesIter != proteinIter->m_vPeptidesNativeIntensity.end(); peptidesIter++) //use vPeptidesNativeIntensity instead of vPeptidesIntensity! by CC 20150527
						{
							//cout << index << endl;			
							dPeptideIntensityTemp = *peptidesIter;
							proteinIter->m_dProteinSPEAQInProtein += dPeptideIntensityTemp;
							proteinIter->m_vPeptidesAdjustIntensity.push_back(dPeptideIntensityTemp);
							proteinIter->m_dProteinSPEAQInPeptides += dPeptideIntensityTemp;
						}

					}

				if (proteinIter->m_iPeptidesNumber != 0)
					{
						proteinIter->m_dProteinSPEAQInPeptides = proteinIter->m_dProteinSPEAQInPeptides / proteinIter->m_iPeptidesNumber;
						proteinIter->m_dLFQ = proteinIter->m_dLFQ / proteinIter->m_iPeptidesNumber;
						PeptidesCorrectionFactorIter = proteinIter->m_vPeptidesCorrectionFactor.begin();
						proteinIter->m_dProteinCorrectionFactor = 0.0;
						for (; PeptidesCorrectionFactorIter != proteinIter->m_vPeptidesCorrectionFactor.end(); PeptidesCorrectionFactorIter++)
						{
							proteinIter->m_dProteinCorrectionFactor += *PeptidesCorrectionFactorIter;
						}

						proteinIter->m_dProteinSPEAQInProtein = proteinIter->m_dProteinSPEAQInProtein / proteinIter->m_dProteinCorrectionFactor;
					}
					else
					{
						proteinIter->m_dProteinSPEAQInPeptides = 0.0;

					}
				}
			}// end for proteins

			if (trainparam.m_bIfCotainStardProtein)
			{
				vecUPS2IDIdentified.clear();
				vecUPS2logMolsIdentified.clear();
				vecUPS2logPredictIdentified.clear();
				for (proteinIter = proteins.begin(); proteinIter != proteins.end(); proteinIter++)
				{
					if (proteinIter->m_bIfCalculateSPEAQInPeptides == true)
					{
						strTemp = proteinIter->m_strProteinID;
						stPosition = strTemp.find("ups");
						if (stPosition != strTemp.npos)
						{
							strTemp = strTemp.substr(0, stPosition);
							//if (strTemp == "P07339")
							//	strTemp = "O76070";
							//if (strTemp == "P08311")
							//	strTemp = "P01579";
							Ups2MolsIter = mapUPS2Id2Mols.find(strTemp);
							if (Ups2MolsIter == mapUPS2Id2Mols.end())
							{
								cout << "Error! Can not find " << strTemp << " in UPS2 File\n";
								log << "Error! Can not find " << strTemp << " in UPS2 File\n";
								continue;
							}
							Ups2WMIter = mapUPS2Id2MW.find(strTemp);

							//cout << strTemp << endl;
							vecUPS2IDIdentified.push_back(strTemp);

							vecUPS2logMolsIdentified.push_back(log10(Ups2MolsIter->second));
							vecUPS2logPredictIdentified.push_back(log10(proteinIter->m_dProteinSPEAQInProtein));

						}

					}

				}

				dUPS2PearsonCorr = spearsonCorrelation(vecUPS2logPredictIdentified, vecUPS2logMolsIdentified);
				cout << "alpha= " << alpha << "beta = " << beta << endl;
				cout << "The Correlation between " << vecUPS2IDIdentified.size() << " Prediction UPS2 and Real UPS2 Mols is " << dUPS2PearsonCorr << endl;
				log << "alpha= " << alpha << "beta = " << beta << endl;
				log << "The Correlation between " << vecUPS2IDIdentified.size() << " Prediction UPS2 and Real UPS2 Mols is " << dUPS2PearsonCorr << endl;

			}
			if (trainparam.m_bIfCotainStardProtein)
			{
				if (dUPS2PearsonCorr > dBestUPS2SpearsonCorr)
				{
					vecSpearsonCorrelation.push_back(dSpearsonCorrelationTemp);
					dBestSpearsonCorrelation = dSpearsonCorrelationTemp;
					vecUPS2SpearsonCorr.push_back(dUPS2PearsonCorr);
					dBestUPS2SpearsonCorr = dUPS2PearsonCorr;
					dbestLamda = lambda;
					dBestNu = nu;
					dBestalpha = alpha;
					dBestBeta = beta;
					dBestm = m;
					cout << "alpha= " << alpha << "beta = " << beta << endl;
					cout << "The Correlation between " << vecUPS2IDIdentified.size() << " Prediction UPS2 and Real UPS2 Mols is " << dUPS2PearsonCorr << endl;
					log << "alpha= " << alpha << "beta = " << beta << endl;
					log << "The Correlation between " << vecUPS2IDIdentified.size() << " Prediction UPS2 and Real UPS2 Mols is " << dUPS2PearsonCorr << endl;

				}
			}
			else
			{
				if (dSpearsonCorrelationTemp > dBestSpearsonCorrelation)
				{
					vecSpearsonCorrelation.push_back(dSpearsonCorrelationTemp);
					dBestSpearsonCorrelation = dSpearsonCorrelationTemp;

					dbestLamda = lambda;
					dBestNu = nu;
					dBestalpha = alpha;
					dBestBeta = beta;
					dBestm = m;
				}
			}

			ParameterSelectByCVResult << lambda << "\t" << kfac << "\t";
			ParameterSelectByCVResult << nu << "\t" << m << "\t";
			ParameterSelectByCVResult << alpha << "\t" << beta << "\t";
			ParameterSelectByCVResult << dSpearsonCorrelationTemp << "\t";
			if (trainparam.m_bIfCotainStardProtein)
			{
				ParameterSelectByCVResult << dUPS2PearsonCorr << "\t";
			}			
			ParameterSelectByCVResult << endl;

		} //end for parameters 

		
	cout << "the best Parameters of Bayes addictiv regression trees is :\n";
	cout << "dbestLamda = " << dbestLamda << endl;
	cout << "dBestNu = " << dBestNu << endl;
	cout << "dBestalpha = " << dBestalpha << endl;
	cout << "dBestBeta = " << dBestBeta << endl;
	log << "the best Parameters of Bayes addictiv regression trees is :\n";
	log << "dbestLamda = " << dbestLamda << endl;
	log << "dBestNu = " << dBestNu << endl;
	log << "dBestalpha = " << dBestalpha << endl;
	log << "dBestBeta = " << dBestBeta << endl;

}

void CBART::mf_saveTestXY(std::ofstream& log, string path)
{
	ofstream ofile;
	ofile.open(path.c_str());
	cout << "Save XY in " << path << endl;
	log << "Save XY in " << path << endl;
	vector<CProtein>::iterator ProteinIter;
	int i = 0;
	int j = 0;
	int t = 0;

	for (int i = 0; i<m_iTestSampleNumber; i++)
	{
		for (int j = 0; j < m_iFeatureNumber; j++)
		{
			ofile << m_vTestx[j + i*m_iFeatureNumber] << "\t";
			
		}
		ofile << m_vTesty[i];
		ofile << "\n";
	}

	ofile.close();

}

void CBART::mf_saveRegressionResult(std::ofstream& log, string path) //输出回归系数与各统计量到文件并显示
{
	int i, j;

	ofstream fout(path.c_str(), ios::out);
	if (!fout)
	{
		cout << "\nCan not open  " << path << endl; 
		log << "Error:\tCan not open  " << path << endl;
		exit(0);
	}
	cout << "Save the BART regression result in " << path << endl;
	log << "Save the BART regression result in " << path << endl;

	double RMSE = 0.0;
	for (int i = 0; i < m_iTrainSampleNumber; i++)
		RMSE += (m_vTrainy[i] - m_vTestPredicty[i])*(m_vTrainy[i] - m_vTrainPredicty[i]);
	RMSE = sqrt(RMSE / m_iTrainSampleNumber);
	fout << "The RMSE of BART regression=" << RMSE << endl;

	double Correlation = spearsonCorrelation(m_vTrainy,m_vTrainPredicty);
	fout << "The Correlation between Prediction Y and Real Y=" << Correlation << endl;
	cout << "The Correlation between Prediction Y and Real Y=" << Correlation << endl;
	log << "The Correlation between Prediction Y and Real Y=" << Correlation << endl;

	fout.close();
	cout << "The result has Saved in " << path << endl;
	log << "The result has Saved in " << path << endl;

}

void CBART::mf_saveTrainXY(std::ofstream& log, string path)
{
	ofstream ofile;
	ofile.open(path.c_str());
	cout << "Save XY in " << path << endl;
	log << "Save XY in " << path << endl;
	vector<CProtein>::iterator ProteinIter;
	int i = 0;
	int j = 0;
	int t = 0;

	for (int i = 0; i<m_iTrainSampleNumber; i++)
	{
		for (int j = 0; j < m_iFeatureNumber; j++)
		{
			ofile << m_vTrainx[j + i*m_iFeatureNumber] << "\t";
		}
		ofile << m_vTrainy[i];
		ofile << "\n";
	}
	cout << "\tEnd of Save " << endl;
	log << "\tEnd of Save " << endl;

	ofile.close();

}


bool CBART::mf_BARTIteration(std::ofstream& log, vector<CProtein> &proteins, size_t burn, size_t nd, size_t m, double lambda, double nu, double kfac, double alpha, double beta)
{

// observe the changed trend of AutoCorrelation and peptides'CV with the iteration;
	ofstream pPeptidesCVFile("PeptidesCV.txt");
	if (!pPeptidesCVFile.is_open())
	{
		cout << "can not open peptidesCV.txt" << endl;
		log << "can not open peptidesCV.txt" << endl;
		exit(0);
	}
	pPeptidesCVFile << "Iterations\tFirstQuantile\tMedian\tThirdQuantile\tAutoCorrelation\tPearsonCorrelationUPS2PredictedWithMols\n";
	double dFirstQuan=0.0;
	double dMedian=0.0;
	double dThirdQuan=0.0;

	mf_ConstructXY(log, proteins);
	mf_BARTAutoRegressionRun(log, proteins, burn, nd, m, lambda, nu, kfac, alpha, beta);
	mf_AdjustPeptidesIntensity(proteins);
	mf_SaveUPS2AnalysisResult(log, proteins, 0);
	mf_CalculatePeptidesCVofProteins(proteins, dFirstQuan, dMedian, dThirdQuan);
	
	pPeptidesCVFile << 0 << "\t" << dFirstQuan << "\t" << dMedian << "\t" << dThirdQuan << "\t";
	pPeptidesCVFile << m_dAutoCorrelation << "\t" << m_dPearsonCorrelationUPS2PredictedWithMols << "\n";

	for (int i = 1; i < 20; i++)
	{
		mf_ConstructXY(log,proteins);
		mf_BARTAutoRegressionRun(log,proteins, burn, nd, m, lambda, nu, kfac, alpha, beta);
		mf_AdjustPeptidesIntensity(proteins);
		mf_SaveUPS2AnalysisResult(log,proteins, i);
		mf_CalculatePeptidesCVofProteins(proteins, dFirstQuan, dMedian, dThirdQuan);
		
		pPeptidesCVFile << i << "\t" << dFirstQuan<<"\t"<<dMedian<<"\t"<<dThirdQuan <<"\t";
		pPeptidesCVFile << m_dAutoCorrelation << "\t" << m_dPearsonCorrelationUPS2PredictedWithMols << "\n";
	}
	return true;
}

bool CBART::mf_BARTAutoRegressionRun(std::ofstream& log, vector<CProtein> &proteins, size_t burn, size_t nd, size_t m, double lambda, double nu, double kfac, double alpha, double beta)
{

	cout << "Bayes addictive regression trees\n";
	cout << "burn,nd,number of trees: " << burn << ", " << nd << ", " << m << endl;
	cout << "lambda,nu,kfac: " << lambda << ", " << nu << ", " << kfac << endl;
	log << "Bayes addictive regression trees\n";
	log << "burn,nd,number of trees: " << burn << ", " << nd << ", " << m << endl;
	log << "lambda,nu,kfac: " << lambda << ", " << nu << ", " << kfac << endl;

	double* pAllfit;//sum of fit of all trees
	double* pResidual; //y-(pAllfit-ftemp) = y-pAllfit+ftemp
	double* pFittemp; //fit of current tree
	//in sample fit
	double* pPredictTrainMean; //posterior mean of in-sample fit, sum draws,then divide

	//x cutpoints
	size_t nc = 100; //100 equally spaced cutpoints from min to max.
	makexinfo(m_iFeatureNumber, m_iTrainSampleNumber, &m_vTrainx[0], xi, nc);

	//trees
	tree cTreeTemp;
	cTreeTemp.setm(m_dMeany / m);
	m_vTrees.clear();
	for (int i = 0; i < m; i++)
	{
		m_vTrees.push_back(cTreeTemp);
	}

	pi.pb = .5; //prob of birth given  birth/death

	pi.alpha = alpha; //prior prob a bot node splits is alpha/(1+d)^beta, d is depth of node
	pi.beta = beta; //2 for bart means it is harder to build big trees.
	pi.tau = (m_dMaxy - m_dMiny) / (2 * kfac*sqrt((double)m));
	pi.sigma = m_dStdy;
	cout << "\nalpha, beta: " << pi.alpha << ", " << pi.beta << endl;
	cout << "sigma, tau: " << pi.sigma << ", " << pi.tau << endl;
	log << "\nalpha, beta: " << pi.alpha << ", " << pi.beta << endl;
	log << "sigma, tau: " << pi.sigma << ", " << pi.tau << endl;

	//--------------------------------------------------
	//dinfo
	pAllfit = new double[m_iTrainSampleNumber]; //sum of fit of all trees
	for (size_t i = 0; i<m_iTrainSampleNumber; i++) pAllfit[i] = m_dMeany;
	pResidual = new double[m_iTrainSampleNumber]; //y-(m_pAllfit-ftemp) = y-m_pAllfit+ftemp
	pFittemp = new double[m_iTrainSampleNumber]; //fit of current tree
	di.n = m_iTrainSampleNumber; di.p = m_iFeatureNumber; di.x = &m_vTrainx[0]; di.y = pResidual; //the y for each draw will be the residual 

	//--------------------------------------------------
	//storage for ouput
	//in sample fit
	pPredictTrainMean = new double[m_iTrainSampleNumber]; //posterior mean of in-sample fit, sum draws,then divide
	for (size_t i = 0; i<m_iTrainSampleNumber; i++) pPredictTrainMean[i] = 0.0;

	//for sigma draw
	double rss;  //residual sum of squares
	double restemp; //a residual
	std::ofstream bsdf("bart-sd.txt"); //note that we write all burn+nd draws to this file.

	std::ofstream btsf("bart-ave-tree-size.txt"); //note that we write all burn+nd draws to this file.
	double ats; //place for average tree size

	double anb; //place for average number of bottom nodes
	std::ofstream bnbf("bart-ave-num-bots.txt"); //note that we write all burn+nd draws to this file.

	//--------------------------------------------------
	//mcmc
	//random number generation
	uint seed = 99;
	RNG gen(seed); //this one random number generator is used in all draws

	cout << "\nMCMC:\n";
	log << "\nMCMC:\n";

	clock_t tp;
	tp = clock();
	for (size_t i = 0; i<(nd + burn); i++)
	{
		if (i % 100 == 0)
		{
			cout << "i: " << i << endl;
			log << "i: " << i << endl;

		}
		//draw trees
		for (size_t j = 0; j<m; j++)
		{
			fit(m_vTrees[j], xi, di, pFittemp);
			for (size_t k = 0; k<m_iTrainSampleNumber; k++)
			{
				pAllfit[k] = pAllfit[k] - pFittemp[k];
				pResidual[k] = m_vTrainy[k] - pAllfit[k];
			}
			bd(log,m_vTrees[j], xi, di, pi, gen);
			drmu(m_vTrees[j], xi, di, pi, gen);
			fit(m_vTrees[j], xi, di, pFittemp);
			for (size_t k = 0; k<m_iTrainSampleNumber; k++) pAllfit[k] += pFittemp[k];
		}
		//draw sigma
		rss = 0.0;
		for (size_t k = 0; k<m_iTrainSampleNumber; k++)
		{
			restemp = m_vTrainy[k] - pAllfit[k];
			rss += restemp*restemp;
		}
		pi.sigma = sqrt((nu*lambda + rss) / gen.chi_square(nu + m_iTrainSampleNumber));
		bsdf << pi.sigma << endl;
		ats = 0.0; anb = 0.0;
		for (size_t k = 0; k<m; k++)
		{
			ats += m_vTrees[k].treesize();
			anb += m_vTrees[k].nbots();
		}
		btsf << ats / m << endl;
		bnbf << anb / m << endl;
		if (i >= burn)
		{
			for (size_t k = 0; k<m_iTrainSampleNumber; k++)
				pPredictTrainMean[k] += pAllfit[k];

		}
	}
	tp = clock() - tp;
	double thetime = (double)(tp) / (double)(CLOCKS_PER_SEC);
	cout << "time for loop: " << thetime << endl;	
	log << "time for loop: " << thetime << endl;

	std::ofstream timef("time.txt");
	timef << thetime << endl;

	m_vTrainPredicty.clear();
	for (size_t k = 0; k < m_iTrainSampleNumber; k++)
	{
		pPredictTrainMean[k] /= (double)nd;
		m_vTrainPredicty.push_back(pPredictTrainMean[k]);
	}

	std::ofstream treef("trees.txt");
	//treef << xi << endl; //the cutpoints
	mf_xitofile(treef, xi);
	treef << m << endl;  //number of trees
	treef << m_iFeatureNumber << endl;  //dimension of x's
	for (size_t j = 0; j<m; j++) treef << m_vTrees[j] << endl;  //all the trees

	m_dAutoCorrelation = spearsonCorrelation(m_vTrainy, pPredictTrainMean);
	cout << "The Correlation between Prediction Y and Real Y is " << m_dAutoCorrelation << endl;
	log << "The Correlation between Prediction Y and Real Y is " << m_dAutoCorrelation << endl;

	Discretization(pPredictTrainMean, m_iTrainSampleNumber);
	//Assertion(pPredictTestMean, m_iTestSampleNumber);

	vector<CProtein>::iterator proteinIter;
	vector<double>::iterator peptidesIter;
	int index = 0;
	for (proteinIter = proteins.begin(); proteinIter != proteins.end(); proteinIter++)
	{
		if (proteinIter->m_iPeptidesNumber >= UniquePeptidesTrainThreshold)
		{
			proteinIter->m_vPeptidesCorrectionFactor.clear();
			for (peptidesIter = proteinIter->m_vPeptidesNativeIntensity.begin(); peptidesIter != proteinIter->m_vPeptidesNativeIntensity.end(); peptidesIter++) //use vPeptidesNativeIntensity instead of vPeptidesIntensity! by CC 20150527
			{
				proteinIter->m_vPeptidesCorrectionFactor.push_back(pPredictTrainMean[index]);
				index++;
			}
		}

	}


	delete pAllfit;
	pAllfit = NULL;
	delete pResidual;
	pResidual = NULL;
	delete pFittemp;
	pFittemp = NULL;
	delete pPredictTrainMean;
	pPredictTrainMean = NULL;

	return 0;


}

bool CBART::mf_CalculatePeptidesCVofProteins(vector<CProtein> proteins, double &FirstQ, double &median, double &ThirdQ)
{//calculate the peptides CV of every protein, then get the meadian and First and third quartile of all CVs;

	vector<double> vecPeptidesCVs;
	vector<double> vecPeptidesIntensityOfOneProtein;

	vector<CProtein>::iterator proteinIter;
	vector<double>::iterator peptidesIter;
	vector<double>::iterator PeptidesCorrectionFactorIter;

	double dPeptideIntensityTemp;
	double dCvTemp;

	//PepAttributeWorker pepworker;
	int i = 0;
	for (proteinIter = proteins.begin(); proteinIter != proteins.end(); proteinIter++)
	{
		i = 0;
		if (proteinIter->m_iPeptidesNumber >= UniquePeptidesTrainThreshold)
		{
			for (peptidesIter = proteinIter->m_vPeptidesAdjustIntensity.begin(); peptidesIter != proteinIter->m_vPeptidesAdjustIntensity.end(); peptidesIter++) //use vPeptidesNativeIntensity instead of vPeptidesIntensity! by CC 20150527
			{
				dPeptideIntensityTemp = *peptidesIter;
				vecPeptidesIntensityOfOneProtein.push_back(dPeptideIntensityTemp);
			}
			dCvTemp = CalculateCV(vecPeptidesIntensityOfOneProtein);
			vecPeptidesIntensityOfOneProtein.clear();
			vecPeptidesCVs.push_back(dCvTemp);
		}
		
	}

	GetQuantiles(vecPeptidesCVs, FirstQ, median, ThirdQ);
	return true;

}

void CBART::mf_SaveUPS2AnalysisResult(std::ofstream& log, vector<CProtein >proteins, int fold)
{

	char Buffer[BUFFERLENGTH];
	char *pstr;
	char *pstr1;
	string strTemp;

	cout << "\tLoading UPS2 mols\n";
	log << "\tLoading UPS2 mols\n";

	string UPS2path = "UPS2_mols.ini";
	map<string, double> mapUPS2Id2Mols;
	map<string, double> mapUPS2Id2MW;
	map<string, double> mapUPS2iBAQs;
	CProteinWorker proteinworker;
	proteinworker.mf_LoadUPS2Mols(log, UPS2path, mapUPS2Id2Mols, mapUPS2Id2MW);

	double dPearsonCorr = 0.0;
	vector<CProtein>::iterator proteinIter;
	map<string, double>::iterator Ups2MolsIter;
	map<string, double>::iterator Ups2WMIter;
	map<string, double>::iterator UPS2IBAQIter;

	double meanReal = 0.0;
	double meanTest = 0.0;
	double sd1 = 0.0, sd2 = 0.0;
	vector<string> vecUPS2IDIdentified;
	vector<double> vecUPS2logPredictIdentified;
	vector<double> vecUPS2logMolsIdentified;

	int iNumberofUPS2Identified;
	string::size_type stPosition;
	string::size_type stStart;

	ofstream ofile;
	char Suffix[10];
	itoa(fold,Suffix,10);
	string strSuffix=Suffix;
	string CorrelationPath ="CheckCorrelationBetweenPredicted_realUPS2Identified_"+strSuffix+".txt";

	ofile.open(CorrelationPath);
	ofile << "UPS2 Protein\tPeptidesNumber\tlogUPS2Mols\tLogProteinSPEAQInPeptides\n";
	for (proteinIter = proteins.begin(); proteinIter != proteins.end(); proteinIter++)
	{
		if (proteinIter->m_bIfCalculateSPEAQInPeptides == true)
		{
			strTemp = proteinIter->m_strProteinID;
			stStart = strTemp.find("|");
			stPosition = strTemp.find("ups");
			if (stPosition != strTemp.npos)
			{
				if (stStart != strTemp.npos)
				{
					strTemp = strTemp.substr(stStart+1, stPosition-stStart-1);
				}
				else
				{
					strTemp = strTemp.substr(0, stPosition);
				}
				
				//if (strTemp == "P07339")
				//	strTemp = "O76070";
				//if (strTemp == "P08311")
				//	strTemp = "P01579";
				Ups2MolsIter = mapUPS2Id2Mols.find(strTemp);
				if (Ups2MolsIter == mapUPS2Id2Mols.end())
				{
					cout << " Can not find " << strTemp << " in UPS2 File\n";
					log << " Can not find " << strTemp << " in UPS2 File\n";

					continue;
				}
				Ups2WMIter = mapUPS2Id2MW.find(strTemp);

				vecUPS2IDIdentified.push_back(strTemp);
				ofile << strTemp << "\t";
				ofile << proteinIter->m_iPeptidesNumber << "\t";
				vecUPS2logMolsIdentified.push_back(log10(Ups2MolsIter->second));
				ofile << log10(Ups2MolsIter->second) << "\t";
				vecUPS2logPredictIdentified.push_back(log10(proteinIter->m_dProteinSPEAQInPeptides));
				ofile << log10(proteinIter->m_dProteinSPEAQInPeptides) << "\n";

			}

		}

	}
	ofile.close();

	dPearsonCorr = spearsonCorrelation(vecUPS2logPredictIdentified,vecUPS2logMolsIdentified);
	m_dPearsonCorrelationUPS2PredictedWithMols = dPearsonCorr;
	cout << "The Correlation between " << vecUPS2IDIdentified.size() << " Prediction UPS2 and Real UPS2 Mols is " << dPearsonCorr << endl;
	log << "The Correlation between " << vecUPS2IDIdentified.size() << " Prediction UPS2 and Real UPS2 Mols is " << dPearsonCorr << endl;

}
