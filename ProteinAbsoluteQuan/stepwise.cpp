#include"stdafx.h"
#include<algorithm>
#include<math.h>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/students_t.hpp>

void Cstepwise::mf_StepwiseInit(std::ofstream& log, vector<CProtein> proteins)
{
	cout << "\tStepwise Initing! \n";
	log << "\tStepwise Initing! \n";
	int i;
	m_iNumberOfFeature = NumberOfFeature;
	vector<CProtein>::iterator ProteinIter;
	m_NumberOfTrainSamples = 0;
	for (ProteinIter = proteins.begin(); ProteinIter != proteins.end(); ProteinIter++)
	{
		m_NumberOfTrainSamples += ProteinIter->m_vPeptidesSequences.size();
	}
	m_px = new double*[m_NumberOfTrainSamples];         // dynamic allocation memory
	for (i = 0; i < m_NumberOfTrainSamples; i++) { m_px[i] = new double[m_iNumberOfFeature + 1]; }
	if (m_px == NULL)
	{
		cout << "m_px 内存分配失败" << endl;
		log << "m_px 内存分配失败" << endl;
	}
	m_pTestXY = new double*[m_NumberOfTrainSamples];         // dynamic allocation memory
	for (i = 0; i < m_NumberOfTrainSamples; i++) { m_pTestXY[i] = new double[m_iNumberOfFeature + 1]; }
	m_pTrainXForCV = new double*[m_NumberOfTrainSamples];         // dynamic allocation memory
	for (i = 0; i < m_NumberOfTrainSamples; i++) { m_pTrainXForCV[i] = new double[m_iNumberOfFeature + 1]; }
	m_pTrainPredictY = new double[m_NumberOfTrainSamples];
	m_pTestPredictY = new double[m_NumberOfTrainSamples];
	m_pNativeTestPredictY = new double[m_NumberOfTrainSamples]; //added by CC 20150527
	m_psortIndex = new int[m_NumberOfTrainSamples];
	m_pLeaveXIndex = new int[m_NumberOfTrainSamples];
	m_pr = new double*[m_iNumberOfFeature + 1];
	for (i = 0; i < m_iNumberOfFeature + 1; i++) { m_pr[i] = new double[m_iNumberOfFeature + 1]; }
	m_pxx = new double[m_iNumberOfFeature + 1];
	m_pCoefficient = new double[m_iNumberOfFeature + 1];
	m_pv = new double[m_iNumberOfFeature + 1];
	m_ps = new double[m_iNumberOfFeature + 1];
	m_pye = new double[m_NumberOfTrainSamples];
	m_pyr = new double[m_NumberOfTrainSamples];
}

bool Cstepwise::mf_Normalize(int i)
{
	if (i == 0)
		return mf_Normalize(m_px, m_NumberOfLimitedSamples, m_iNumberOfFeature);
	else if (i == 1)
		return mf_Normalize(m_pTestXY, m_NumberOfTrainSamples, m_iNumberOfFeature);
	else
		return false;
}

bool Cstepwise::mf_Normalize(double **p, int m, int n)            //using ZScore； //change the function class from void to bool by ChangCheng 20150527
{
	if (n < 1 || m < 3) //added by ChangCheng 20150527 
		return false;

	double *meanXY = new double[n + 1];
	double *sdeXY = new double[n + 1];
	for (int i = 0; i < n + 1; i++)
		meanXY[i] = 0.0;
	for (int i = 0; i < m; i++)   //Change m_NumberOfTrainSamples by m gzhq 20150528
		for (int j = 0; j < n + 1; j++)
			meanXY[j] += p[i][j];
	for (int j = 0; j < n + 1; j++)
		meanXY[j] = meanXY[j] / m;

	for (int i = 0; i < n + 1; i++)
		sdeXY[i] = 0.0;
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n + 1; j++)
			sdeXY[j] += (p[i][j] - meanXY[j])*(p[i][j] - meanXY[j]);
	for (int j = 0; j < n + 1; j++)
		sdeXY[j] = sqrt(sdeXY[j] / (m - 1));

	m_dMeanY = meanXY[n];        //add by gzhq 20150604
	m_dSdeY = sdeXY[n];
	//standarization
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n + 1; j++)
			p[i][j] = (p[i][j] - meanXY[j]) / sdeXY[j];
	//no centralization
	//for (int i = 0; i < m_NumberOfTrainSamples; i++)
	//	for (int j = 0; j < n + 1; j++)
	//		p[i][j] = p[i][j] / sdeXY[j];
	delete[]meanXY;
	delete[] sdeXY;
	return true;
}


void Cstepwise::mf_StepwiseRegression(std::ofstream& log, double ff1, double ff2, double es)//执行逐步回归分析
{
	cout << "\tstepwise Begin:\n";
	log << "\tstepwise Begin:\n";
	int i, j, ii, m, imi, imx, l, it;
	double z, phi, sd, vmi, vmx, q, fmi, fmx;
	m_dthresh1 = ff1;  m_dthresh2 = ff2;  m_dEps = es;
	m = m_iNumberOfFeature + 1; q = 0.0;
	double *Tempsd = new double[m_iNumberOfFeature + 1];
	for (j = 0; j <= m_iNumberOfFeature; j++)
	{
		z = 0.0;
		for (i = 0; i <= m_NumberOfLimitedSamples - 1; i++)
			z = z + m_px[i][j];                   
		m_pxx[j] = z / m_NumberOfLimitedSamples;
	}
	for (i = 0; i <= m_iNumberOfFeature; i++)
		for (j = 0; j <= i; j++)
		{
			z = 0.0;
			for (ii = 0; ii <= m_NumberOfLimitedSamples - 1; ii++)
				z = z + (m_px[ii][i] - m_pxx[i])*(m_px[ii][j] - m_pxx[j]);
			m_pr[i][j] = z;
		}
	for (i = 0; i <= m_iNumberOfFeature; i++)  Tempsd[i] = sqrt(m_pr[i][i]);// the standard deviation of ith feature
	for (i = 0; i <= m_iNumberOfFeature; i++)
		for (j = 0; j <= i; j++)
		{
			m_pr[i][j] = m_pr[i][j] / (Tempsd[i] * Tempsd[j]);
			m_pr[j][i] = m_pr[i][j];
		}
	phi = m_NumberOfLimitedSamples - 1.0;
	sd = Tempsd[m_iNumberOfFeature] / sqrt(m_NumberOfLimitedSamples - 1.0);  //the standard devivation of last feature;
	it = 1;
	while (it == 1)
	{
		it = 0;
		vmi = 1.0e+35; vmx = 0.0;
		imi = -1; imx = -1;
		for (i = 0; i <= m_iNumberOfFeature; i++)
		{
			m_pv[i] = 0.0;
			m_pCoefficient[i] = 0.0; 
			m_ps[i] = 0.0;
		}
		for (i = 0; i <= m_iNumberOfFeature - 1; i++)
			if (m_pr[i][i] >= m_dEps)
			{
				m_pv[i] = m_pr[i][m_iNumberOfFeature] * m_pr[m_iNumberOfFeature][i] / m_pr[i][i];
				if (m_pv[i] >= 0.0)
				{
					if (m_pv[i]>vmx) { vmx = m_pv[i]; imx = i; }
				}
				else
				{
					m_pCoefficient[i] = m_pr[i][m_iNumberOfFeature] * Tempsd[m_iNumberOfFeature] / Tempsd[i];// ？
					m_ps[i] = sqrt(m_pr[i][i])*sd / Tempsd[i];// ？
					if (fabs(m_pv[i])<vmi)
					{
						vmi = fabs(m_pv[i]); imi = i;
					}
				}
			}
		if (phi != m_iNumberOfFeature - 1.0)
		{
			z = 0.0;
			for (i = 0; i <= m_iNumberOfFeature - 1; i++)  
				z = z + m_pCoefficient[i] * m_pxx[i];
			m_pCoefficient[m_iNumberOfFeature] = m_pxx[m_iNumberOfFeature] - z;
			m_ps[m_iNumberOfFeature] = sd; 
			m_pv[m_iNumberOfFeature] = q;
		}
		else
		{
			m_pCoefficient[m_iNumberOfFeature] = m_pxx[m_iNumberOfFeature]; m_ps[m_iNumberOfFeature] = sd;
		}
		fmi = vmi*phi / m_pr[m_iNumberOfFeature][m_iNumberOfFeature];
		fmx = (phi - 1.0)*vmx / (m_pr[m_iNumberOfFeature][m_iNumberOfFeature] - vmx);
		if ((fmi<m_dthresh2) || (fmx >= m_dthresh1))
		{
			if (fmi<m_dthresh2)  
			{ 
				phi = phi + 1.0;
				l = imi; 
			    cout << "Delete variable " << l << "\n"; 
				log << "Delete variable " << l << "\n";
			}
			else 
			{ 
				phi = phi - 1.0; 
				l = imx; 
				cout << "Add variable " << l << "\n";
				log << "Add variable " << l << "\n";
			}
			for (i = 0; i <= m_iNumberOfFeature; i++)
				if (i != l)
					for (j = 0; j <= m_iNumberOfFeature; j++)
						if (j != l)
							m_pr[i][j] = m_pr[i][j] - (m_pr[l][j] / m_pr[l][l])*m_pr[i][l];
			for (j = 0; j <= m_iNumberOfFeature; j++)
				if (j != l) m_pr[l][j] = m_pr[l][j] / m_pr[l][l];
			for (i = 0; i <= m_iNumberOfFeature; i++)
				if (i != l)  m_pr[i][l] = -m_pr[i][l] / m_pr[l][l];
			m_pr[l][l] = 1.0 / m_pr[l][l];
			q = m_pr[m_iNumberOfFeature][m_iNumberOfFeature] * Tempsd[m_iNumberOfFeature] * Tempsd[m_iNumberOfFeature];
			sd = sqrt(m_pr[m_iNumberOfFeature][m_iNumberOfFeature] / phi)*Tempsd[m_iNumberOfFeature];
			m_dC = sqrt(1.0 - m_pr[m_iNumberOfFeature][m_iNumberOfFeature]);
			m_dF = (phi*(1.0 - m_pr[m_iNumberOfFeature][m_iNumberOfFeature])) / ((m_NumberOfLimitedSamples - phi - 1.0)*m_pr[m_iNumberOfFeature][m_iNumberOfFeature]);
			it = 1;
		}
	}
	for (i = 0; i <= m_NumberOfLimitedSamples - 1; i++)
	{
		z = 0.0;
		for (j = 0; j <= m_iNumberOfFeature - 1; j++) z = z + m_pCoefficient[j] * m_px[i][j];
		m_pye[i] = m_pCoefficient[m_iNumberOfFeature] + z; m_pyr[i] = m_px[i][m_iNumberOfFeature] - m_pye[i];
	}
	delete[]Tempsd;
}

void Cstepwise::mf_saveStepwiseXY(std::ofstream& log, string path)
{
	ofstream ofile;
	ofile.open(path.c_str());
	cout << "Save XY in " << path << endl;
	log << "Save XY in " << path << endl;
	vector<CProtein>::iterator ProteinIter;
	int i = 0;
	int j = 0;
	int t = 0;
	for (int i = 0; i<m_NumberOfLimitedSamples; i++)
	{
		for (int j = 0; j <= m_iNumberOfFeature; j++)
			ofile << m_px[i][j] << "\t";
		ofile << "\n";
	}

	cout << "\tEnd of Save " << endl;
	log << "\tEnd of Save " << endl;
	ofile.close();

}

void Cstepwise::mf_saveRegressionResult(std::ofstream& log, string path) //输出回归系数与各统计量到文件并显示
{
	int i, j;

	ofstream fout(path.c_str(), ios::out);
	if (!fout)
	{
		cout << "Can not open " << path << endl;
		log << "Error:\tCan not open " << path << endl;
		exit(0);
	}
	cout << "Save the regression result in " << path << endl;
	log << "Save the regression result in " << path << endl;
	fout << endl;
	fout << "m_dthresh1 = " << m_dthresh1 << "    " << "m_dthresh2 = " << m_dthresh2 << endl;
	fout << "m_dEps = " << m_dEps << endl;
	fout << "square of deviance and residual sum of squares:" << endl;
	for (i = 0; i<m_iNumberOfFeature; i++)
	{
		fout << "Ex(" << i << ")=" << m_pxx[i] << "    ";
	}
	fout << "\nEy=" << m_pxx[m_iNumberOfFeature] << endl;
	fout << "regression coefficients:" << endl;
	for (i = 0; i<m_iNumberOfFeature; i++)
	{
		fout << m_pCoefficient[i] << "\t";
	}
	fout << m_pCoefficient[i];
	fout << endl << endl;
	fout << "sum of squares of partial regression:" << endl;
	for (i = 0; i<m_iNumberOfFeature; i++)
	{
		fout << m_pv[i] << "\t";
	}
	fout << endl << "residual sum of squares = " << m_pv[m_iNumberOfFeature] << endl << endl;
	fout << "the standard deviation of regression coefficients:" << endl;
	for (i = 0; i<m_iNumberOfFeature; i++)
	{
		fout << m_ps[i] << "\t";
	}
	fout << endl << "the estimated SD = " << m_ps[m_iNumberOfFeature]  << endl;
	fout << "multiple correlation coefficientsC = " << m_dC << endl << endl;
	fout << "the value of F-test = " << m_dF << endl << endl;
	fout << "estimated values of the conditional expectation of independent variables and ";
	fout<<" residuals of observed variables:" << endl;
	for (i = 0; i<m_NumberOfLimitedSamples; i++)
	{
		fout << m_pye[i] << "\t" << m_pyr[i] << endl;;
	}
	fout << endl;
	fout << "m_NumberOfLimitedSamples residuals of observed variables:" << endl;
	for (i = 0; i <= m_iNumberOfFeature; i++)
	{
		for (j = 0; j <= m_iNumberOfFeature; j++)
		{
			fout << m_pr[i][j] << " ";
		}
		fout << endl;
	}

	for (int i = 0; i < m_NumberOfLimitedSamples; i++)
	{
		m_pTrainPredictY[i] = 0.0;
		for (int j = 0; j < m_iNumberOfFeature; j++)
			m_pTrainPredictY[i] += m_pCoefficient[j] * m_px[i][j];

		m_pTrainPredictY[i] += m_pCoefficient[m_iNumberOfFeature];
	}

	m_dRMSE = 0.0;
	for (int i = 0; i < m_NumberOfLimitedSamples; i++)
		m_dRMSE += (m_px[i][m_iNumberOfFeature] - m_pTrainPredictY[i])*(m_px[i][m_iNumberOfFeature] - m_pTrainPredictY[i]);
	m_dRMSE = sqrt(m_dRMSE / m_NumberOfLimitedSamples);
	fout << "The RMSE of stepwise regression is " << m_dRMSE << endl;

	ofstream ofile;
	ofile.open("checkCorrelationRealYAndY");
	cout << "Begin Save x and Y" << endl;
	log << "Begin Save x and Y" << endl;
	for (int i = 0; i<m_NumberOfLimitedSamples; i++)
	{
		ofile << m_px[i][m_iNumberOfFeature] << "\t" << m_pTrainPredictY[i] << endl;
	}
	ofile.close();

	m_dCorrelation = 0.0;
	double meanReal = 0.0;
	double meanTest = 0.0;
	double sd1 = 0.0, sd2 = 0.0;

	meanTest = Average(m_pTrainPredictY, m_NumberOfLimitedSamples);
	for (int i = 0; i < m_NumberOfLimitedSamples; i++)
		meanReal += m_px[i][m_iNumberOfFeature];
	meanReal = meanReal / m_NumberOfLimitedSamples;
	for (int i = 0; i < m_NumberOfLimitedSamples; i++)
		m_dCorrelation += (m_px[i][m_iNumberOfFeature] - meanReal)*(m_pTrainPredictY[i] - meanTest);
	for (int i = 0; i < m_NumberOfLimitedSamples; i++)
	{
		sd1 += (m_px[i][m_iNumberOfFeature] - meanReal)*(m_px[i][m_iNumberOfFeature] - meanReal);
		sd2 += (m_pTrainPredictY[i] - meanTest)*(m_pTrainPredictY[i] - meanTest);
	}
	m_dCorrelation = m_dCorrelation / sqrt(sd1*sd2);
	fout << "The Correlation between Prediction Y and Real Y is " << m_dCorrelation << endl;
	cout << "The Correlation between Prediction Y and Real Y is " << m_dCorrelation << endl;
	log << "The Correlation between Prediction Y and Real Y is " << m_dCorrelation << endl;

	fout.close();
	cout << "The result has Saved in " << path << endl;
	log << "The result has Saved in " << path << endl;

}
bool Cstepwise::mf_ConstructXY(std::ofstream& log, vector<CProtein> proteins)
{
	cout << "\tConstructing XY\n\n";	
	log << "\tConstructing XY\n\n";

	vector<CProtein>::iterator  ProteinIter;
	int j = 0; //feature num
	int n = 0; //peptide num per protein
	m_NumberOfLimitedSamples = 0;  //Change i by m_NumberOfLimitedSamples gzhq 20150528
	m_iNumberOfTrainProteins = 0;

	for (ProteinIter = proteins.begin(); ProteinIter != proteins.end(); ProteinIter++)
	{
		if (ProteinIter->m_iPeptidesNumber >= UniquePeptidesTrainThreshold)
		{
			m_iNumberOfTrainProteins++;
			for (n = 0; n < ProteinIter->m_iPeptidesNumber; n++)
			{
				for (j = 0; j < ProteinIter->m_vPeptidesAttributes.at(n).m_vecAttributes.size(); j++)
				{
					m_px[m_NumberOfLimitedSamples][j] = ProteinIter->m_vPeptidesAttributes.at(n).m_vecAttributes.at(j);
				}

				m_px[m_NumberOfLimitedSamples][j] = ProteinIter->m_vAdjustPeptidesIntensityFactor.at(n);
				m_NumberOfLimitedSamples++;
			}

		}

	}
	cout << "Constructed XY using " << m_NumberOfLimitedSamples << " Peptides in "<<m_iNumberOfTrainProteins <<" Proteins."<< endl; // added by CC 20150527
	log << "Constructed XY using " << m_NumberOfLimitedSamples << " Peptides in " << m_iNumberOfTrainProteins << " Proteins." << endl; // added by CC 20150527
	return true;
}

bool Cstepwise::mf_loadModel(std::ofstream& log, string m_strRegressionResult)
{
	cout << "\tLoading model!\n\n";
	log << "\tLoading model!\n\n";

	ifstream fin(m_strRegressionResult.c_str());
	if (!fin)
	{
		cout << "Can not open " << m_strRegressionResult << endl;
		log << "Error:\tCan not open " << m_strRegressionResult << endl;
		exit(0);
	}
	string strTemp1, strTemp2;
	int iTemp1 = 0, iTemp2 = 0, i = 0;

	while (getline(fin, strTemp1))
	{
		strTemp2 = strTemp1.substr(0, 8);
		if (strTemp2 == "回归系数")
		{
			getline(fin, strTemp1);

			iTemp1 = 0;
			iTemp2 = strTemp1.find("\t", iTemp1);
			while (iTemp2 != strTemp1.npos)
			{

				strTemp2 = strTemp1.substr(iTemp1, iTemp2 - iTemp1);
				m_pCoefficient[i] = atof(strTemp2.c_str());
				i++;
				iTemp1 = iTemp2 + 1;
				iTemp2 = strTemp1.find("\t", iTemp1);
			}
			strTemp2 = strTemp1.substr(iTemp1, strTemp1.size() - iTemp1);
			m_pCoefficient[i] = atof(strTemp2.c_str());
			cout << "The Model has " << i << " coefficients.\n";
			log << "The Model has " << i << " coefficients.\n";

			return true;
		}

	}
	return false;
}
void Cstepwise::mf_saveTestXY(std::ofstream& log, string path)
{
	ofstream ofile;
	ofile.open(path.c_str());
	cout << "Save TestXY in " << path << endl;
	log << "Save TestXY in " << path << endl;

	vector<CProtein>::iterator ProteinIter;
	int i = 0;
	int j = 0;
	int t = 0;

	for (int i = 0; i<m_iTestSampleNumber; i++)
	{
		for (int j = 0; j <= m_iNumberOfFeature; j++)
			ofile << m_pTestXY[i][j] << "\t";
		ofile << "\n";
	}

	cout << "\tEnd of Save " << endl;
	log << "\tEnd of Save " << endl;

	ofile.close();

}

void Cstepwise::mf_saveStepwiseX(std::ofstream& log, vector<CProtein> proteins, string path)
{
	ofstream ofile;
	ofile.open(path.c_str());
	cout << "Save XY in " << path << endl;
	log << "Save XY in " << path << endl;

	vector<CProtein>::iterator ProteinIter;
	int i = 0;
	int j = 0;
	int t = 0;

	for (ProteinIter = proteins.begin(); ProteinIter != proteins.end(); ProteinIter++)
	{

		for (j = 0; j < ProteinIter->m_iPeptidesNumber; j++)
		{
			ofile << ProteinIter->m_strProteinID << "\t";
			ofile << ProteinIter->m_vPeptidesSequences.at(j) << "\t";
			for (t = 0; t <= m_iNumberOfFeature; t++)
				ofile << m_pTestXY[i][t] << "\t";
			ofile << endl;
			i++;
		}

	}
	ofile.close();

}

void Cstepwise::mf_LoadTestXSPEAQInPeptides(std::ofstream& log, string SPEAQInPeptidesPathByOtherMethod)
{
	cout << "\tLoading SPEAQInPeptidess!\n\n";
	log << "\tLoading SPEAQInPeptidess!\n\n";

	ifstream fin(SPEAQInPeptidesPathByOtherMethod.c_str());
	if (!fin)
	{
		cout << "Can not open " << SPEAQInPeptidesPathByOtherMethod << endl;
		log << "Error:\tCan not open " << SPEAQInPeptidesPathByOtherMethod << endl;
		exit(0);

	}
	string strTemp1, strTemp2;
	int iTemp1 = 0, iTemp2 = 0, i = 0;
	int iSamplesTemp = 0;

	while (getline(fin, strTemp1))
	{
		m_pTestPredictY[iSamplesTemp] = atof(strTemp1.c_str());
		m_pNativeTestPredictY[iSamplesTemp] = m_pTestPredictY[iSamplesTemp];
		iSamplesTemp++;
	}
	if (iSamplesTemp != m_iTestSampleNumber)
	{
		cout << "The Number of SPEAQInPeptidess does not correspond to TestX's" << endl;
		log << "The Number of SPEAQInPeptidess does not correspond to TestX's" << endl;
		exit(0);
	}

}


/*
Calculate peptides quantification coefficient 
*/
void Cstepwise::mf_Prediction(vector<CProtein> &proteins)
{
	
	//************ To test another regression method **************
	//string SPEAQInPeptidesPathByOtherMethod = "E:\\Rcode\\regress\\MARS\\PredictY.txt";
	//string SPEAQInPeptidesPathByOtherMethod = "E:\\ccode\\ProteinAbsoluteQuan_20150816\\ProteinAbsoluteQuan\\Predict\\NativeY.txt";
	//string SPEAQInPeptidesPathByOtherMethod = "E:\\ccode\\ProteinAbsoluteQuan_20150816\\ProteinAbsoluteQuan\\Predict\\NativeYDivideByMax.txt";
	//string SPEAQInPeptidesPathByOtherMethod = "E:\\ccode\\ProteinAbsoluteQuan_20150816\\ProteinAbsoluteQuan\\Predict\\PredictYByBART.txt";
	//string SPEAQInPeptidesPathByOtherMethod = "E:\\ccode\\ProteinAbsoluteQuan_20150816\\ProteinAbsoluteQuan\\Predict\\PredictYByBARTMedian.txt";
	//string SPEAQInPeptidesPathByOtherMethod = "E:\\ccode\\ProteinAbsoluteQuan_20150816\\ProteinAbsoluteQuan\\PredictY.txt";
	//string SPEAQInPeptidesPathByOtherMethod = "E:\\ccode\\BARTCVersion\\BARTCVersion\\bart-fit.txt";
	//string SPEAQInPeptidesPathByOtherMethod = "E:\\ccode\\ProteinAbsoluteQuan20150906\\ProteinAbsoluteQuan\\ProteinAbsoluteQuan\\BARTPredictY.txt";
	//mf_LoadTestXSPEAQInPeptides(SPEAQInPeptidesPathByOtherMethod);

	
	for (int i = 0; i < m_iTestSampleNumber; i++)
	{
		m_pTestPredictY[i] = 0.0;
		for (int j = 0; j < m_iNumberOfFeature; j++)
		{
			m_pTestPredictY[i] += m_pCoefficient[j] * m_pTestXY[i][j];  

		}
		m_pTestPredictY[i] += m_pCoefficient[m_iNumberOfFeature];
	
	}
	for (int i = 0; i < m_iTestSampleNumber; i++)
	{
		m_pNativeTestPredictY[i] = m_pTestPredictY[i]; //added by CC 20150527
	}
	


	//need discretization
	Discretization(m_pTestPredictY, m_iTestSampleNumber);
	//Assertion(m_pTestPredictY, m_iTestSampleNumber);

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
				proteinIter->m_vPeptidesCorrectionFactor.push_back(m_pTestPredictY[index]);
				index++;
			}
		}
			
	}

}

bool Cstepwise::mf_ConstructTestX(std::ofstream& log, vector<CProtein> proteins)
{                 //add this function by gzhq 20150603
	cout << "\tConstructing TestXY\n";
	log << "\tConstructing TestXY\n";
	vector<CProtein>::iterator  ProteinIter;
	int j = 0; //feature num
	int n = 0; //peptide num per protein
	m_iTestSampleNumber = 0;
	m_iNumberOfTestProteins = 0;

	map<string, CPeptideAttribute> mapPeptideSequenceAndAttribute;


	for (ProteinIter = proteins.begin(); ProteinIter != proteins.end(); ProteinIter++)
	{
		
		if (ProteinIter->m_iPeptidesNumber >= UniquePeptidesTestThreshold)
		{
			m_iNumberOfTestProteins++;
			for (n = 0; n < ProteinIter->m_iPeptidesNumber; n++)
			{
				for (j = 0; j < ProteinIter->m_vPeptidesAttributes.at(n).m_vecAttributes.size(); j++)
				{
					m_pTestXY[m_iTestSampleNumber][j] = ProteinIter->m_vPeptidesAttributes.at(n).m_vecAttributes.at(j);
				}
				m_pTestXY[m_iTestSampleNumber][j] = ProteinIter->m_vAdjustPeptidesIntensityFactor.at(n);  //change m_vnativeintensity to m_vPeptidesIntensity by gzhq 20150706
				m_iTestSampleNumber++;
			}
		}
	}
	cout << "ConstructTestXY using " << m_iTestSampleNumber <<" peptides in "<<m_iNumberOfTestProteins<<" proteins. " <<endl; // added by CC 20150527
	log << "ConstructTestXY using " << m_iTestSampleNumber << " peptides in " << m_iNumberOfTestProteins << " proteins. " << endl; // added by CC 20150527
	return true;
}

void Cstepwise::mf_Predict(std::ofstream& log, vector<CProtein> &proteins, string RegressionResultPath)
{
	mf_loadModel(log,RegressionResultPath);
	mf_ConstructTestX(log,proteins);

	//string checkXpath = "CheckPredictX.txt";
	//mf_saveStepwiseX(log,proteins, checkXpath);

	int iPosTemp = RegressionResultPath.find_last_of("\\");
	iPosTemp = RegressionResultPath.find_last_of("\\", iPosTemp - 1);
	mf_Prediction(proteins);
	mf_ShowPredictionAnalysis(log);

}

void Cstepwise::mf_ShowPredictionAnalysis(std::ofstream& log)
{
	m_dCorrelation = 0.0;
	double meanReal = 0.0;
	double meanTest = 0.0;
	double sd1 = 0.0, sd2 = 0.0;

	ofstream ofile;
	ofile.open("checkCorrelationYandPredictedY.txt");
	cout << " Save Y and Predicted Y for checking in checkCorrelationYandPredictedY.txt " << endl;
	log << " Save Y and Predicted Y for checking in checkCorrelationYandPredictedY.txt " << endl;

	for (int i = 0; i<m_iTestSampleNumber; i++)
	{
		ofile << m_pTestXY[i][m_iNumberOfFeature] << "\t" << m_pNativeTestPredictY[i] << endl;
	}
	ofile.close();

	meanTest = Average(m_pNativeTestPredictY, m_iTestSampleNumber);
	for (int i = 0; i < m_iTestSampleNumber; i++)
		meanReal += m_pTestXY[i][m_iNumberOfFeature];
	meanReal = meanReal / m_iTestSampleNumber;
	for (int i = 0; i < m_iTestSampleNumber; i++)
		m_dCorrelation += (m_pTestXY[i][m_iNumberOfFeature] - meanReal)*(m_pNativeTestPredictY[i] - meanTest);
	for (int i = 0; i < m_iTestSampleNumber; i++)
	{
		sd1 += (m_pTestXY[i][m_iNumberOfFeature] - meanReal)*(m_pTestXY[i][m_iNumberOfFeature] - meanReal);
		sd2 += (m_pNativeTestPredictY[i] - meanTest)*(m_pNativeTestPredictY[i] - meanTest);
	}
	m_dCorrelation = m_dCorrelation / sqrt(sd1*sd2);
	cout << "The Correlation between Prediction Y and Real Y is " << m_dCorrelation << endl;
	log << "The Correlation between Prediction Y and Real Y is " << m_dCorrelation << endl;

}

//int Cstepwise::mf_LeaveOneProteinOutConstructXY(int l)
//{//返回进入训练集的样本数目；留出第l个蛋白的肽段做测试集，其他的做训练集；
//	vector<CProtein>::iterator  ProteinIter;
//	int i = 0;
//	int j = 0;
//	int n = 0;
//	int iIndexTemp = 0;
//	int ProteinNum = 0;
//	m_iTestSampleNumber = 0;
//
//	//iLeaveXNumber = 0;
//	for (ProteinIter = m_vecProteins.begin(); ProteinIter != m_vecProteins.end(); ProteinIter++)
//	{
//		if (ProteinIter->m_iPeptidesNumber >= UniquePeptidesThreshold)  //added by gzhq20150529
//		{
//			if (ProteinNum != l)
//			{
//				for (n = 0; n < ProteinIter->m_iPeptidesNumber; n++)
//				{
//					for (j = 0; j < ProteinIter->m_vPeptidesAttributes.at(n).m_vecAttributes.size(); j++)
//					{
//						m_pTrainXForCV[i][j] = ProteinIter->m_vPeptidesAttributes.at(n).m_vecAttributes.at(j);
//					}
//
//					m_pTrainXForCV[i][j] = ProteinIter->m_vAdjustPeptidesIntensityFactor.at(n);
//					i++;
//					iIndexTemp++;//记录留出的蛋白对应的肽段在训练集中的位置
//				}
//			}
//			else
//			{
//				//cout << "留蛋白质 " << l << endl;
//				for (n = 0; n < ProteinIter->m_iPeptidesNumber; n++)
//				{
//					for (j = 0; j < ProteinIter->m_vPeptidesAttributes.at(n).m_vecAttributes.size(); j++)
//					{
//						m_pTestXY[m_iTestSampleNumber][j] = ProteinIter->m_vPeptidesAttributes.at(n).m_vecAttributes.at(j);
//					}
//
//					m_pTestXY[m_iTestSampleNumber][j] = ProteinIter->m_vAdjustPeptidesIntensityFactor.at(n);
//					m_pLeaveXIndex[m_iTestSampleNumber] = iIndexTemp;
//
//					m_iTestSampleNumber++;
//					iIndexTemp++;
//				}
//
//			}
//		}
//
//
//		ProteinNum++;
//		//iLeaveXNumber = m_iTestSampleNumber;
//
//	}
//
//	return i - 1;
//}

//double Cstepwise::mf_LeaveOneOutCV(double ff1, double ff2, double es)
//{
//	cout << "\tLeave One Out Crossing Validation\n";
//	int i, j, ii, m, imi, imx, l, it;
//	double z, phi, sd, vmi, vmx, q, fmi, fmx;
//	m_dthresh1 = ff1;  m_dthresh2 = ff2;  m_dEps = es;
//	int n = 0;
//	double *TestYForCV = new double[m_NumberOfLimitedSamples];   
//	for (int i = 0; i < m_NumberOfLimitedSamples; i++)
//		TestYForCV[i] = 0.0;
//	double *Tempsd = new double[m_iNumberOfFeature + 1];
//	for (int iProteinNum = 0; iProteinNum < m_vecProteins.size(); iProteinNum++)
//	{
//
//		n = mf_LeaveOneProteinOutConstructXY(iProteinNum);
//
//		//cout << "\t" << m_iTestSampleNumber << endl;
//		m = m_iNumberOfFeature + 1; q = 0.0;
//
//		for (j = 0; j <= m_iNumberOfFeature; j++)
//		{
//			z = 0.0;
//			for (i = 0; i <= n - 1; i++) z = z + m_pTrainXForCV[i][j] / n;
//			m_pxx[j] = z;
//		}
//		for (i = 0; i <= m_iNumberOfFeature; i++)
//			for (j = 0; j <= i; j++)
//			{
//				z = 0.0;
//				for (ii = 0; ii <= n - 1; ii++)
//					z = z + (m_pTrainXForCV[ii][i] - m_pxx[i])*(m_pTrainXForCV[ii][j] - m_pxx[j]);
//				m_pr[i][j] = z;
//			}
//		for (i = 0; i <= m_iNumberOfFeature; i++)  Tempsd[i] = sqrt(m_pr[i][i]);
//		for (i = 0; i <= m_iNumberOfFeature; i++)
//			for (j = 0; j <= i; j++)
//			{
//				m_pr[i][j] = m_pr[i][j] / (Tempsd[i] * Tempsd[j]);
//				m_pr[j][i] = m_pr[i][j];
//			}
//		phi = n - 1.0;
//		sd = Tempsd[m_iNumberOfFeature] / sqrt(n - 1.0);  
//		it = 1;
//		while (it == 1)
//		{
//			it = 0;
//			vmi = 1.0e+35; vmx = 0.0;
//			imi = -1; imx = -1;
//			for (i = 0; i <= m_iNumberOfFeature; i++)
//			{
//				m_pv[i] = 0.0; m_pCoefficient[i] = 0.0; m_ps[i] = 0.0;
//			}
//			for (i = 0; i <= m_iNumberOfFeature - 1; i++)
//				if (m_pr[i][i] >= m_dEps)
//				{
//					m_pv[i] = m_pr[i][m_iNumberOfFeature] * m_pr[m_iNumberOfFeature][i] / m_pr[i][i];
//					if (m_pv[i] >= 0.0)
//					{
//						if (m_pv[i] > vmx) { vmx = m_pv[i]; imx = i; }
//					}
//					else
//					{
//						m_pCoefficient[i] = m_pr[i][m_iNumberOfFeature] * Tempsd[m_iNumberOfFeature] / Tempsd[i];// ？
//						m_ps[i] = sqrt(m_pr[i][i])*sd / Tempsd[i];// ？
//						if (fabs(m_pv[i]) < vmi)
//						{
//							vmi = fabs(m_pv[i]); imi = i;
//						}
//					}
//				}
//			if (phi != m_iNumberOfFeature - 1.0)
//			{
//				z = 0.0;
//				for (i = 0; i <= m_iNumberOfFeature - 1; i++)  z = z + m_pCoefficient[i] * m_pxx[i];
//				m_pCoefficient[m_iNumberOfFeature] = m_pxx[m_iNumberOfFeature] - z; m_ps[m_iNumberOfFeature] = sd; m_pv[m_iNumberOfFeature] = q;
//			}
//			else
//			{
//				m_pCoefficient[m_iNumberOfFeature] = m_pxx[m_iNumberOfFeature]; m_ps[m_iNumberOfFeature] = sd;
//			}
//			fmi = vmi*phi / m_pr[m_iNumberOfFeature][m_iNumberOfFeature];
//			fmx = (phi - 1.0)*vmx / (m_pr[m_iNumberOfFeature][m_iNumberOfFeature] - vmx);
//			if ((fmi < m_dthresh2) || (fmx >= m_dthresh1))
//			{
//				if (fmi < m_dthresh2)  { phi = phi + 1.0; l = imi; /*cout << "变量" << l << "被剔除\n";*/ }
//				else  { phi = phi - 1.0; l = imx; /*cout << "变量" << l << "加入\n";*/ }
//				for (i = 0; i <= m_iNumberOfFeature; i++)
//					if (i != l)
//						for (j = 0; j <= m_iNumberOfFeature; j++)
//							if (j != l)
//								m_pr[i][j] = m_pr[i][j] - (m_pr[l][j] / m_pr[l][l])*m_pr[i][l];
//				for (j = 0; j <= m_iNumberOfFeature; j++)
//					if (j != l) m_pr[l][j] = m_pr[l][j] / m_pr[l][l];
//				for (i = 0; i <= m_iNumberOfFeature; i++)
//					if (i != l)  m_pr[i][l] = -m_pr[i][l] / m_pr[l][l];
//				m_pr[l][l] = 1.0 / m_pr[l][l];
//				q = m_pr[m_iNumberOfFeature][m_iNumberOfFeature] * Tempsd[m_iNumberOfFeature] * Tempsd[m_iNumberOfFeature];
//				sd = sqrt(m_pr[m_iNumberOfFeature][m_iNumberOfFeature] / phi)*Tempsd[m_iNumberOfFeature];
//				m_dC = sqrt(1.0 - m_pr[m_iNumberOfFeature][m_iNumberOfFeature]);
//				m_dF = (phi*(1.0 - m_pr[m_iNumberOfFeature][m_iNumberOfFeature])) / ((n - phi - 1.0)*m_pr[m_iNumberOfFeature][m_iNumberOfFeature]);
//				it = 1;
//			}
//		}//while //added by CC 20150527
//		cout << iProteinNum << endl; //for debug by CC 20150527
//		//delete[]Tempsd;  //del by CC 20150527
//		for (i = 0; i <= n - 1; i++)
//		{
//			z = 0.0;
//			for (j = 0; j <= m_iNumberOfFeature - 1; j++) z = z + m_pCoefficient[j] * m_pTrainXForCV[i][j];
//			m_pye[i] = m_pCoefficient[m_iNumberOfFeature] + z; m_pyr[i] = m_pTrainXForCV[i][m_iNumberOfFeature] - m_pye[i];
//		}
//		for (int i = 0; i < m_iTestSampleNumber; i++)
//		{
//			for (int j = 0; j < m_iNumberOfFeature; j++)
//				TestYForCV[m_pLeaveXIndex[i]] += m_pCoefficient[j] * m_pTrainXForCV[m_pLeaveXIndex[i]][j];
//		}
//
//	}//for (int iProteinNum = 0; iProteinNum < m_vecProteins.size(); iProteinNum++)
//	delete[]Tempsd;  //added by CC 20150527
//	Tempsd = NULL; //added by CC 20150527
//
//	m_dCorrelation = 0.0;
//	double meanReal = 0.0;
//	double meanTest = 0.0;
//	double sd1 = 0.0, sd2 = 0.0;
//	meanTest = Average(TestYForCV, m_NumberOfLimitedSamples);
//	for (int i = 0; i < m_NumberOfLimitedSamples; i++)
//		meanReal += m_pTrainXForCV[i][m_iNumberOfFeature];
//	meanReal = meanReal / m_NumberOfLimitedSamples;
//	for (int i = 0; i < m_NumberOfLimitedSamples; i++)
//		m_dCorrelation += (m_pTrainXForCV[i][m_iNumberOfFeature] - meanReal)*(TestYForCV[i] - meanTest);
//	for (int i = 0; i < m_NumberOfLimitedSamples; i++)
//	{
//		sd1 += (m_pTrainXForCV[i][m_iNumberOfFeature] - meanReal)*(m_pTrainXForCV[i][m_iNumberOfFeature] - meanReal);
//		sd2 += (TestYForCV[i] - meanTest)*(TestYForCV[i] - meanTest);
//	}
//	m_dCorrelation = m_dCorrelation / sqrt(sd1*sd2);
//
//	delete[] TestYForCV;
//	if (sd1*sd2 <= 0) // added by CC 20150527
//		return -1;
//	else	return m_dCorrelation;
//
//}

//bool Cstepwise::mf_CheckAttributeValueByChangchengsData(string infilePath, string outFilePath, string strProteinFasterFilePath, string strProteinFilePath)
//{
//	mf_LoadProteinIDAndSequence(m_mapProteinIDAndSequence, strProteinFasterFilePath);
//
//	FILE * pFile;
//	pFile = fopen(infilePath.c_str(), "r");
//	if (pFile == NULL)
//	{
//		cout << "Error when opening \"" << strProteinFilePath << endl;
//		return false;
//	}
//	ofstream ofile;
//	ofile.open(outFilePath.c_str());
//	cout << "\tBegin Save" << endl;
//	ofile << "peptideSequence\tProteinID\tAttributes\n";
//
//
//	//cout<<"open "<<strProteinFilePath<<endl;
//	char Buffer[BUFFERLENGTH];
//	char *pstr;
//	char *pstr1;
//	fgets(Buffer, BUFFERLENGTH, pFile);    //
//	int k = 0;
//	string strTemp1, strTemp2;
//	map<string, string>::iterator iter;
//	vector<double > Attributes;
//	double dIntensityTemp = 0.0;
//	map<string, double>::iterator str2doubleIter;
//
//	CPeptideAttribute peptideAttributeTemp;
//	PepAttributeWorker PeptideAttributeWorkerTemp;
//	string::size_type peptideLocation = 0, begin = 0, Getwidth = 0;
//	m_vecAAindexAttributeHeaders = PeptideAttributeWorkerTemp.mf_LoadAAindex();
//
//	while (!feof(pFile))
//	{
//		//cout << "k =" << k << endl;
//		Attributes.clear();
//		k++;
//		pstr = Buffer;
//		pstr1 = strstr(pstr, "\t");
//		*pstr1 = '\0';
//		strTemp1 = pstr;
//		ofile << pstr << "\t";
//		pstr = pstr1 + 1;
//
//		//pstr1 = strstr(pstr, "\t");
//		//pstr = pstr1 + 1;
//
//		//pstr1 = strstr(pstr, "\t");
//		//pstr = pstr1 + 1;
//		//pstr1 = strstr(pstr, "\t");
//		//pstr = pstr1 + 1;
//		//pstr1 = strstr(pstr, "\t");
//		//pstr = pstr1 + 1;
//
//		pstr1 = strstr(pstr, "\n");
//		*pstr1 = '\0';
//		//cout << pstr << endl;
//		iter = m_mapProteinIDAndSequence.find(pstr);
//		if (iter != m_mapProteinIDAndSequence.end())
//		{
//			strTemp2 = iter->second;
//			//str2doubleIter = m_Peptides.PeptideSequenceAndIntensity.find(strTemp1);
//			//if (str2doubleIter != m_Peptides.PeptideSequenceAndIntensity.end())
//			//{
//			//	dIntensityTemp = str2doubleIter->second;
//			//}
//
//
//			peptideLocation = 0, begin = 0, Getwidth = 0;
//			peptideLocation = iter->second.find(strTemp1);
//			if (peptideLocation != iter->second.npos)
//			{
//
//				if ((peptideLocation >= FlankingRegionWidth) && (peptideLocation <= iter->second.size() - strTemp1.size() - FlankingRegionWidth))
//				{
//					begin = peptideLocation - FlankingRegionWidth;
//					Getwidth = strTemp1.size() + 2 * FlankingRegionWidth;
//				}
//				else if (peptideLocation < FlankingRegionWidth)
//				{
//					begin = 0;
//					Getwidth = peptideLocation + strTemp1.size() + FlankingRegionWidth;
//				}
//				else if (peptideLocation > strTemp1.size() - strTemp1.size() - FlankingRegionWidth)
//				{
//					begin = peptideLocation - FlankingRegionWidth;
//					Getwidth = iter->second.size() - peptideLocation - 1 + strTemp1.size() + FlankingRegionWidth;
//				}
//				strTemp2 = iter->second.substr(begin, Getwidth);
//
//				ofile << pstr << "\t";
//				Attributes = PeptideAttributeWorkerTemp.mf_GetAttributeFromSequence(strTemp1, strTemp2);
//				for (int j = 0; j < Attributes.size(); j++)
//					ofile << Attributes.at(j) << "\t";
//			}
//		}
//		ofile << endl;
//		fgets(Buffer, BUFFERLENGTH, pFile);
//	}
//
//	fclose(pFile);
//	ofile.close();
//}

Cstepwise:: ~Cstepwise()
{
	int i;
	for (i = 0; i < m_NumberOfTrainSamples; i++) { delete[m_iNumberOfFeature + 1] m_px[i]; m_px[i] = NULL; }
	delete[m_NumberOfTrainSamples] m_px;
	m_px = NULL;
	for (i = 0; i < m_NumberOfTrainSamples; i++) { delete[m_iNumberOfFeature + 1] m_pTestXY[i]; m_pTestXY[i] = NULL; }
	for (i = 0; i < m_NumberOfTrainSamples; i++) { delete[m_iNumberOfFeature + 1] m_pTrainXForCV[i]; m_pTrainXForCV[i] = NULL; }
	delete[] m_pTrainXForCV;
	m_pTrainXForCV = NULL;

	delete[] m_pLeaveXIndex;
	delete[]m_pTrainPredictY;
	delete[]m_psortIndex;
	for (i = 0; i < m_iNumberOfFeature; i++) { delete[m_iNumberOfFeature + 1] m_pr[i]; m_pr[i] = NULL; }
	delete[m_NumberOfTrainSamples] m_pr;
	m_pr = NULL;
	delete[] m_pxx;
	delete[m_iNumberOfFeature+1]m_pCoefficient;
	delete[]m_pv;
	delete[]m_ps;
	delete[]m_pye;
	delete[]m_pyr;

	for (i = 0; i < m_iTestSampleNumber; i++) { delete[m_iNumberOfFeature + 1] m_pTestXY[i]; m_pTestXY[i] = NULL; }
	delete[m_iTestSampleNumber] m_pTestXY;
	m_pTestXY = NULL;
	delete[]m_pTestPredictY;
	delete[]m_pNativeTestPredictY; //added by CC 20150527

	cout << "*********************Program end！***************************\n";

}


