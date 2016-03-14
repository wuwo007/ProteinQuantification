#include"BasicClass.h"
#include<string>
using namespace std;

class Cstepwise
{
private:

	double  **m_px, m_dthresh1, m_dthresh2, m_dEps;
	double  *m_pxx;// The average value of indenpendent variables and dependent variables
	double **m_pTrainXForCV;  //Storing X in train dataset for Cross Validation£»added by gzhq 20150529
	double * m_pTrainPredictY;
	int *m_psortIndex;
	int *m_pLeaveXIndex;
	double *m_pCoefficient; //regression coefficients
	double *m_pv;//sum of squares of partial regression and residual sum of squares
	double *m_ps;  // the standard deviation of regression coefficients and the estimated SD
	double m_dC;  //multiple correlation coefficients
	double m_dF;  //the value of F-test
	double *m_pye;// m_NumberOfLimitedSamples estimated values of the conditional expectation of independent variables
	double *m_pyr;//m_NumberOfLimitedSamples residuals of observed variables
	double **m_pr;  //normalized coefficient correlation matrix

	int m_iTestSampleNumber;//the number of samples in test dataset for cross validation
	double m_dCorrelation;
	double m_dBestCorrelation;

	double m_dRMSE;
	double m_dMeanY;
	double m_dSdeY;

	int m_iNumberOfFeature;

	double **m_pTestXY;      //Storing X in test dataset for Cross Validation and prediction£»  change TestX to TestXY by gzhq 20150604
	double *m_pTestPredictY;
	double *m_pNativeTestPredictY; //to save the predicted y before discretization added by CC 20150527 ;change name m_pNativeTestY to m_pNativeTestPredictY by gzhq 20150625

	map<string, string> m_mapProteinIDAndSequence;

public:
	int m_iNumberOfTrainProteins;
	int m_iNumberOfTestProteins;
	int m_NumberOfTrainSamples;
	int m_NumberOfLimitedSamples;  //the number of samples limited by unique peptides threshold; Added by gzhq 20150528
	Cstepwise()
	{
		m_NumberOfTrainSamples = 0;
		m_NumberOfLimitedSamples = 0;
	}

	void mf_StepwiseInit(std::ofstream& log, vector<CProtein> proteins);


	bool mf_Normalize(double **p, int m, int n);// Normalization For X and Y in int **x; change the param by gzhq20150603
	bool mf_Normalize(int n);

	void mf_StepwiseRegression(std::ofstream& log, double, double, double); 
	void mf_showXY(std::ofstream& log );  
	void mf_saveStepwiseXY(std::ofstream& log, string path);

	bool mf_ConstructXY(std::ofstream& log, vector<CProtein> proteins);
	bool mf_ConstructTestX(std::ofstream& log);
	void mf_saveRegressionResult(std::ofstream& log, string path);  


	void mf_Predict(std::ofstream& log, vector<CProtein> &proteins, string regressionResultPath);
	bool mf_loadModel(std::ofstream& log, string m_strRegressionResult);

	bool mf_ConstructTestX(std::ofstream& log, vector<CProtein> proteins);
	void mf_saveStepwiseX(std::ofstream& log, vector<CProtein> proteins, string path);
	void mf_LoadTestXSPEAQInPeptides(std::ofstream& log, string SPEAQInPeptidesPathByOtherMethod);

	void mf_Prediction(vector<CProtein> &proteins);  
	void mf_saveTestXY(std::ofstream& log, string path);
	void mf_ShowPredictionAnalysis(std::ofstream& log); 
	void mf_ShowUPS2AnalysisResult(std::ofstream& log, vector<CProtein> proteins);



	
	//bool mf_CheckAttributeValueByChangchengsData(string infilePath, string outFilePath);

	~Cstepwise();
};
