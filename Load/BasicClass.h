#pragma once
#include<map>
#include<string>
#include <vector>
using namespace std;
#define BUF_LENGTH 10000
const int FlankingRegionWidth = 15;
const int NumberOfFeature = 587;
class CProtein
{
public:
	CProtein();
	void Clear();
	string m_strProteinID;
	string m_strProteinFullName;
	vector<int> m_vPeptidesIDs;
	vector<string> m_vPeptidesSequences;
	map<string,vector<string>> m_mapPeptidesSequenceAndCorrespondRawFiles;
	map < string, vector<double>> m_mapExperimentAndPeptidesIntensity;
	double m_dProteinSPEAQInPeptides;
	vector<double> m_vecdMaxquantIBAQ;
	vector<double> m_vecdMaxquantLFQ;
	double m_dCVNative;
	double m_dCVAfterAdjust;
	double m_dPeptidesIntensitySum;
	int m_iPeptidesNumber;
	int m_iNumberOfTheoreticEnzyme;

	string m_strProteinSequence;
	string mf_GetPeptidesAdjacentSequence(std::ofstream& log, string PeptideSequence);
};
class Cpeptides
{
public:
	void mf_Clear();
	map<int, string> m_mapPeptideIDAndSequence;
	map<int, vector<double>> m_mapPeptideIDAndIntensitys;
	map<string, vector<double>> m_mapPeptideSequenceAndIntensitys;
	map<string, bool> m_mapPeptideSequencuAndIfUse;
};


double Average(double *Array, int Number);
double Average(vector<double> Array);

double spearsonCorrelation(double s1[], double s2[]);
double CalculateCV(vector<double>s);
double IntegralFromVector(std::ofstream& log, vector<double> Values, vector<double> times);

bool fStringToBool(string str);

string strim(string &str);