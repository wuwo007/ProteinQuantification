#include <vector>
#include <ctime>

#include "rng.h"
#include "tree.h"
#include "info.h"
#include "funs.h"
#include"BasicClass.h"
#include"QuantificationParam.h"

class CBART
{

public:
	// Iteration for better result
	bool mf_BARTIteration(std::ofstream& log, vector<CProtein> &proteins, size_t burn = 100, size_t nd = 100, size_t m = 200, double lambda = 1.0, double nu = 3.0, double kfac = 2.0, double alpha = 0.95, double beta = 1.0);
	bool mf_BARTAutoRegressionRun(std::ofstream& log, vector<CProtein> &proteins, size_t burn = 100, size_t nd = 100, size_t m = 200, double lambda = 1.0, double nu = 3.0, double kfac = 2.0, double alpha = 0.95, double beta = 1.0);
	bool mf_CalculatePeptidesCVofProteins(vector<CProtein> proteins, double &FirstQ, double &median, double &ThirdQ);
	void mf_SaveUPS2AnalysisResult(std::ofstream& log, vector<CProtein >proteins, int i);
	void mf_xitofile(std::ofstream& fnm, xinfo xi); 
	bool mf_ConstructXY(std::ofstream& log, vector<CProtein> proteins );
	bool mf_ConstructTestX(std::ofstream& log, vector<CProtein> proteins);
	bool mf_BartRegressionRun(std::ofstream& log, vector<CProtein>&proteins, size_t burn = 100, size_t nd = 100, size_t m = 200, double lambda = 1.0, double nu = 1.0, double kfac = 2.0, double alpha = 0.95, double beta = 1.0);
	void mf_AdjustPeptidesIntensity(vector<CProtein>& proteins);
	void mf_saveRegressionResult(std::ofstream& log, string path);  

	void mf_SelectFeature(std::ofstream& log, vector<tree> trees, map<size_t, size_t> &FeatureCount);
	void mf_saveTestXY(std::ofstream& log, string path);
	void mf_saveTrainXY(std::ofstream& log, string path);

	// for ten-fold cross-validation
	void mf_ChooseBestParameter(std::ofstream& log, vector<CProtein> &proteins, CTrainParam trainparam);
	void mf_MakeIndexPartition(size_t SampleNumber,vector<size_t> &IndexPartation); //make index partation for dividing dataset
	void mf_ConstructCVData( std::ofstream& log, vector<CProtein> proteins, vector<size_t>IndexPartation, size_t fold, vector<double> &TrainxForCV, vector<double> &TrainYForCV, vector<double>&TestxForCV, vector<double>&TestyForCV);
	void mf_TrainAndPredictForCV(std::ofstream& log, vector<CProtein> &proteins, vector<double> &TrainxForCV, vector<double> &TrainYForCV, vector<double>&TestxForCV, vector<double> &PredictedY, size_t burn = 100, size_t nd = 100, size_t m = 200, double lambda = 1.0, double nu = 3.0, double kfac = 2.0, double alpha = 0.95, double beta = 1.0);

	double m_dBestAlpha;
	double m_dBestBeta;
	size_t m_iBestNumberOfTrees;
	string m_strQuantiResultPath;
	string m_strSelectedFeaturePath;

private:
	std::vector<double> m_vTrainx;
	std::vector<double> m_vTrainy;
	size_t m_iFeatureNumber;
	size_t m_iTrainSampleNumber;
	double m_dAutoCorrelation;

	double m_dMiny;
	double m_dMaxy;
	sinfo m_cAllys;

	double m_dMeany;
	double	m_dStdy;

	std::vector<double> m_vTestx;     //prediction observations, stored like x
	std::vector<double> m_vTesty;
	std::vector<double> m_vTestPredicty;
	vector<double> m_vTrainPredicty;
	dinfo dip; //data information for prediction
	size_t m_iTestSampleNumber;
	size_t m_iNumberOfTrainProteins;
	size_t m_iNumberOfTestProteins;



	xinfo xi;
	size_t nc = 100; //100 equally spaced cutpoints from min to max.

	//trees
	std::vector<tree> m_vTrees;

	//prior and mcmc
	pinfo pi;
	dinfo di;

	//for sigma draw
	double m_dRss;  //residual sum of squares
	double m_dRestemp; //a residual

	double m_dPearsonCorrelationUPS2PredictedWithMols;

	
};