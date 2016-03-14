#include"QuantificationParam.h"

// this class is in charge of regression to adjust peptides intensity;
class CRegression
{
public:
	void mf_RegressionRun(std::ofstream& log, vector<CProtein> &proteins, CTrainParam trainparam, CPredictParam predicParam, string strExperimentName);

	void mf_StepwiseRegression(std::ofstream& log, vector<CProtein> &proteins, CTrainParam trainparam);
	void mf_BARTRegression(std::ofstream& log, vector<CProtein> &proteins, CTrainParam trainparam, string strExperimentName);

};
