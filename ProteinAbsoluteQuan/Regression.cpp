#include "stdafx.h"
#include "Regression.h"
#include"BARTRegression.h"
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/students_t.hpp>

void CRegression::mf_StepwiseRegression(std::ofstream& log, vector<CProtein> &proteins, CTrainParam trainparam)
{
	Cstepwise stepwise;
	stepwise.mf_StepwiseInit(log,proteins);
	stepwise.mf_ConstructXY(log, proteins);

	//int iPosTemp = trainparam.m_strRegressionResult.find_last_of("\\");
	//iPosTemp = trainparam.m_strRegressionResult.find_last_of("\\", iPosTemp-1);
	//string TrainXYpath = trainparam.m_strRegressionResult.substr(0,iPosTemp+1)+"TrainXY.txt";
	//stepwise.mf_saveStepwiseXY(log, TrainXYpath);

	//string CheckXYpath = "CheckXY.txt";
	//stepwise.mf_saveStepwiseXY(CheckXYpath);

	double f1 = 0.0, f2 = 0.0;
	/*boost::math::students_t T1_distribution(stepwise.m_NumberOfLimitedSamples - 1);
	boost::math::students_t T2_distribution(stepwise.m_NumberOfLimitedSamples);
	f1 = boost::math::quantile(T1_distribution, trainParam.m_dAlpha1);
	f2 = boost::math::quantile(T2_distribution, trainParam.m_dAlpha2);*/
	boost::math::fisher_f_distribution<double> F1_distribution(1, stepwise.m_NumberOfLimitedSamples - 1);
	boost::math::fisher_f_distribution<double> F2_distribution(1, stepwise.m_NumberOfLimitedSamples);
	f1 = boost::math::quantile(F1_distribution, trainparam.m_dAlpha1);
	f2 = boost::math::quantile(F2_distribution, trainparam.m_dAlpha2);
	cout << "\tf1= " << f1 << " \tf2= " << f2 << endl;
	log << "\tf1= " << f1 << " \tf2= " << f2 << endl;

	//cout << "\tbest multiple correlation coefficient is " << bestParam.m_dBestCorrelation << endl;
	//string CheckXYpathnormalized = "CheckXYnormalized.txt";
	//stepwise.mf_saveStepwiseXY(CheckXYpathnormalized);
	stepwise.mf_StepwiseRegression(log, f1, f2, 1.0e-32); //执行逐步回归分析
	stepwise.mf_saveRegressionResult(log, trainparam.m_strRegressionResult);       //输出回归系数与各统计量到文件并显示

	stepwise.mf_Predict(log, proteins, trainparam.m_strRegressionResult);

}

void CRegression::mf_BARTRegression(std::ofstream& log, vector<CProtein> &proteins, CTrainParam trainparam, string strExperimentName)
{
	
	CBART bartRegression;
	bartRegression.m_strQuantiResultPath = trainparam.m_strQuantiResultPath;
	bartRegression.m_strSelectedFeaturePath = bartRegression.m_strQuantiResultPath + "\\SelectedFeature"+strExperimentName+".txt";

	if (trainparam.m_bIfOptimizeParameters)
	{
		// CrossValidation for parameters
		bartRegression.mf_ConstructXY(log,proteins);
		bartRegression.mf_ChooseBestParameter(log,proteins,trainparam);
		bartRegression.mf_AdjustPeptidesIntensity(proteins);
	}
	else
	{
		bartRegression.m_dBestAlpha = trainparam.m_dAlpha;
		bartRegression.m_dBestBeta = trainparam.m_dBeta;
		bartRegression.m_iBestNumberOfTrees = trainparam.m_iNumberOfTrees;
	}

	//bartRegression.mf_BARTIteration(proteins,100, 100, 200, 1.0, 3.0, 2.0, bartRegression.m_dBestAlpha, bartRegression.m_dBestBeta);
	bartRegression.mf_ConstructXY(log, proteins);
	bartRegression.mf_ConstructTestX(log, proteins);
	bartRegression.mf_BartRegressionRun(log, proteins, 100, 100, bartRegression.m_iBestNumberOfTrees, 1.0, 3.0, 2.0, bartRegression.m_dBestAlpha, bartRegression.m_dBestBeta);
	string strRegressionResultPath = trainparam.m_strRegressionResult + strExperimentName + ".txt";
	bartRegression.mf_saveRegressionResult(log, strRegressionResultPath);

		
}

void CRegression::mf_RegressionRun(std::ofstream& log, vector<CProtein> &proteins, CTrainParam trainparam, CPredictParam predicParam, string strExperimentName)
{
	if (trainparam.m_strRegressionMethod == "stepwise")
	{
		mf_StepwiseRegression(log,proteins, trainparam);
	}
	else
	{
		mf_BARTRegression(log, proteins, trainparam, strExperimentName);
	}
}

