#pragma once

#include<map>
#include<string>
#include"PepAttributeWorker.h"
#include <vector>
using namespace std;


class CProtein
{
public:
	CProtein();
	void Clear();
	string m_strProteinID;
	string m_strProteinFullName;
	vector<int> m_vPeptidesIDs;
	vector<string> m_vPeptidesSequencesWithFlankingRegion;
	vector<string> m_vPeptidesSequences;
	vector<double> m_vdPeptidesMW;
	map<string, vector<string>> m_mapPeptidesSequenceAndCorrespondRawFiles;
	vector<CPeptideAttribute> m_vPeptidesAttributes;
	vector<double> m_vAdjustPeptidesIntensityFactor; //modified a spelling error by CC 20150527
	vector<double> m_vPeptidesNativeIntensity; //modified a spelling error by CC 20150527
	vector<double> m_vPeptidesAdjustIntensity;
	vector<double> m_vPeptidesCorrectionFactor;
	double m_dPeptidesIntensityMedian;
	double m_dPeptidesIntensityMax;
	double m_dProteinSPEAQInPeptides;
	double m_dProteinSPEAQInProtein;
	double m_dProteinCorrectionFactor;
	double m_dMaxquantIBAQ;
	double m_dReCalculateIBAQ;
	double m_dLFQ;
	double m_dCVNative;
	double m_dCVAfterAdjustAfter;
	double m_dPeptidesIntensitySum;
	int m_iPeptidesNumber;
	int m_iNumberOfTheoreticEnzyme;
	bool m_bIfCalculateSPEAQInPeptides;

	string m_strProteinSequence;
	string mf_GetPeptidesAdjacentSequence(std::ofstream& log, string PeptideSequence);
};
class Cpeptides
{
public:
	void mf_Clear();
	map<int, string> m_mapPeptideIDAndSequence;
	map<int, double> m_mapPeptideIDAndIntensity;
	map<string, double> m_mapPeptideSequenceAndIntensity;
	map<string, bool> m_mapPeptideSequencuAndIfUse;
	map<string, CPeptideAttribute> m_mapPeptideSequenceAndAttribute;
	string m_strPeptideIntensityName;

};


class CGene
{
public:
	map<string, vector<string>> m_mapProteinIDAndPeptidesSequenceWithMod;
	map<string, vector<double>> m_mapPRoteinIDAndPeptidesIntensitys;
	map<string, double> m_mapProteinIDAndIntensity;
	map<string, int> m_mapProteinIdAndDigestNumber;
	double m_dExpressionIntensity;
	string m_strgeneID;
	void mf_CaluculateIntensityFromProteinsIntensity1();
	void mf_CaluculateIntensityFromProteinsIntensity2();
	void mf_Clear()
	{
		m_dExpressionIntensity = 0.0;
		m_mapProteinIDAndPeptidesSequenceWithMod.clear();
		m_mapPRoteinIDAndPeptidesIntensitys.clear();
		m_mapProteinIdAndDigestNumber.clear();
		m_mapProteinIDAndIntensity.clear();
	}

};

class CMergedProtein
{
public:
	void Clear();
	string m_strPrtoteinName;
	vector<string> m_vecExperiments;
	vector<double> m_vecMaxquantIBAQOfExperiments;
	vector<double> m_vecSPEAQInPeptidesOfExperiments;
	vector<double> m_vecMaxquantLFQOfExperiments;
	vector<double> m_vecSPEAQInProteinOfExperiments;
};

class CMergedProteins
{
public:
	void mf_ExperimentsAnalysis();
	vector<CMergedProtein> m_vecMergedProteins;
	vector<string> m_strExpreriments;
	vector<double> m_vecLFQCV;
	vector<double> m_vecMaxquantIBAQCV;
	vector<double> m_vecSPEAQInPeptidesCV;
	vector<double> m_vecSPEAQInProteinCV;
};

double CalculateCV(vector<double>s);
double Average(double *Array, int Number);
double Average(vector<double> Array);
void Discretization(double *TestY, int m_iTestSampleNumber);
void Discretization(vector<double> &TestY);
double GetMaxIntensity(const vector<double> PeptidesIntensity);
double GetMaxIntensity(vector<double> PeptidesIntensity);
double GetMedianIntensity(vector<double> PeptidesIntensity);

void GetQuantiles(vector<double> PeptidesIntensity,double &FirstQ,double &median,double &ThirdQ);

//Ensure that there is not negative and zeros;
void Assertion(vector<double> &TestY);
void Assertion(double *TestY, int m_iTestSampleNumber);
void quickSort(double s[], int index[], int l, int r);  // A problem: Stack overflow may happen in the sort function. 
void quickSort(vector<size_t> &s, vector<size_t> &index, int l, int r);

double spearsonCorrelation(double s1[], double s2[],int numberOfSample);
double spearsonCorrelation(vector<double> s1, double s2[]);
double spearsonCorrelation(vector<double> s1, vector<double> s2);

void GetAttributesFromFirstRow(char * pstr, map<string, int>& mapAttributesAndColumns);

bool fStringToBool(string str);