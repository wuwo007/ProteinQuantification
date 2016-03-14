#include"BasicClass.h"
#include"QuantificationParam.h"
// this class implements the proteis quantification 
class CProteinWorker
{
public:
	CProteinWorker()
	{
		m_bIfIBAQIntensityExist = false; //change IfIBAQIntensityExit to IfIBAQIntensityExist by ChangCheng 20150527
		m_iPeptidesNumberOfProteins = 0;
	}



	bool mf_LoadProteinFasta(std::ofstream& log, string strProteinFasterFilePath, FastaType fastatype);
	bool mf_LoadProteins(std::ofstream& log, CTrainParam trainParam, int iExperimentIndex);
	bool mf_LoadProteinSPEAQInPeptides(std::ofstream& log, string path);
	bool mf_LoadUPS2Mols(std::ofstream& log, string path, map<string, double>&mapUPS2Id2Mols, map<string, double>&mapUPS2Id2MW);
	bool mf_AnalysisPeptidesCVForSameUPS2Mols(std::ofstream& log, string path, vector<CProtein>& vecProteins, map<string, double>&mapUPS2Id2Mols);

	int mf_GetExperimentNames(std::ofstream& log, string experimentalDesignPath);
	bool mf_CalculateProteinIntensity(std::ofstream& log, int iExperimentIndex);
	void mf_CalculateProteinIBAQ(std::ofstream& log);
	void mf_CalculateProteinSPEAQInPeptides(std::ofstream& log, int iExperimentIndex);


	void mf_SaveProteinsNativePeptides(std::ofstream& log, string path);
	void mf_saveProteinSPEAQInPeptides(std::ofstream& log, string ProteinSPEAQInPeptidesResultPath);
	void mf_MergeProteinSPEAQInPeptides(std::ofstream& log, string ProteinQsoreresultPath);
	void mf_MergeRegressionResult(std::ofstream& log, string RegressionResultPath);
	void mf_ShowUPS2AnalysisResult(std::ofstream& log, int iExperimentIndex);
	void mf_Clear();

	vector<string> m_vecStrExperiments;  // experiment names in order;
	CTrainParam m_trainParam;
	CPredictParam m_predictParam;
	vector<CProtein> m_vecProteins;
	
private:
	int m_iPeptidesNumberOfProteins;  //The number of peptides of proteins loaded
	map<string, string> m_mapProteinsIdAndSequence;
	vector<string> m_vecAAindexAttributeHeaders;

	int m_iNumberOfExperiments;
	bool m_bIfIBAQIntensityExist;
	bool m_bIFLFQIntensityExist;
	 

};
