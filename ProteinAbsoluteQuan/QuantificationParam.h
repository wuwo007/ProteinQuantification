
#pragma once

//File version: 2015-12-04
/*this class contain all parameters which used in the stage of training and predicting
*/
#include<string>
using namespace std;
class CQuantificationParam
{
public:
	string m_strProteinsPath;
	string m_strRegressionResult;
	string strIndentifySoftwareType;
	string m_strProteinFastaPath;
	FastaType m_eFastaType;
	//string m_strResultPath;
	string m_strMaxquantResultPath;
	string m_strQuantiResultPath;
	string m_strExperimentDesignPath;
	bool m_bIfCotainStardProtein;
	string m_strIdentifierOfStandPro;
	string m_strAbundanceOfStandProteinsPath;
};
class CTrainParam:public CQuantificationParam
{
public:
	int mf_setparameters(std::ofstream& log, string parafilepath);

	// for stepwise
	double m_dAlpha1;
	double m_dAlpha2;
	// for BART 
	double m_dAlpha;
	double m_dBeta;
	double m_iNumberOfTrees;
	int m_ik;
	string m_strRegressionMethod;

	bool m_bIfOptimizeParameters;
};

class CPredictParam:public CQuantificationParam
{
public:
	void mf_setparameters(std::ofstream& log, string parafilepath);


	string m_strProteinSPEAQInPeptidesPath;
	string m_strCutPosition;//the position of the protein where enzyme works
	bool m_blr;//cut to the left position or to the right, T-->left,  False-->right
	int m_imaxMissedClevage;//the maxium missed clevage allowed
	int m_iMinPepLength;//the required minimum length of the peptides, which is used to select peptides, known as L1 in [L1, L2]
	int m_iMaxPepLength;//the required maxium length of the peptides, which is used to select peptides, known as L2 in [L1, L2]

};