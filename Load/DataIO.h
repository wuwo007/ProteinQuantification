#pragma once
#include"BasicClass.h"
#include"loadParam.h"
#include"stdafx.h"
#include<map>
#include<string>
using namespace std;

class CLoadIO
{
public:
	CLoadIO();
	vector<CProtein> m_vecProteins;
	Cpeptides m_Peptides;
	//vector<string> m_vecAAindexAttributeHeaders;
	bool m_bIfIBAQIntensityExist;
	bool m_bIfLFQIntensityExist;
	bool m_bIfExistReverseProtein;
	bool m_bIfExistContanProtein;
	string m_strReversePrefix;
	string m_strContaminantfix;
	int m_NumberOfPeptides;
	int m_iNumberOfExperiments;
	map<string, string> m_mapProteinIDAndSequence;
	map<string, int> m_mapPeptideSequenceAndAppearNumber;

	string mf_GetPeptidesSequenceFromID(std::ofstream& log, int ID);
	double mf_GetPeptidesIntensityFromID(std::ofstream& log, int ID, int iExperimentIndex);
	double mf_GetMedianIntensity(std::ofstream& log, vector<double> PeptidesIntensity);
	void mf_saveProteins(std::ofstream& log, string proteinspath);
	//void mf_savePeptides(string peptidesFeaturePath);
	bool mf_LoadProteinIDAndSequence(std::ofstream& log, map<string, string> &IdAndSequence, string strProteinFasterFilePath, FastaType fastatype);
	
	
	bool mf_LoadPeptides(std::ofstream& log, string strPeptideFilePath);
	bool mf_LoadProteins(std::ofstream& log, string strProteinFilePath, string strProteinFasterFilePath);//从proteinGroups.txt载入蛋白质，并且在fasta文件中取得相应的序列

	void mf_GetAttributesFromFirstRow(char * pstr, map<string, int>& mapAttributesAndColumns);

};

class CLoadMaxquantIO :public CLoadIO
{
public:
	CLoadMaxquantIO(std::ofstream& log, const CLoadParam loadParam);

	int mf_GetPeptideAttributecolumn(std::ofstream& log, string AttributeName);
	
	int mf_GetProteinAttributecolumn(std::ofstream& log, map<string, int>mapAttrtibuteAndcolumns,
									int iReversecolumnNum, int iPoteinalContaminantcolumnNum, 
									int iPeptideIDscolumnNum, int  iPeptideRazorcolumnNum, 
									vector<int> veciBAQcolumnNum ,vector<int> veciLFQintensitycolumnNum);  //analysis the first row to get the attribute columns;
	void mf_GetAttributesName(std::ofstream& log, string ExperimentDesignPath);

	//string m_strIBAQ_IntensityName;//IBAQ_intrensity、LFQ_intensity、Peptides_intensity对应的列名不固定，因此需根据输入文件而定；
	//string m_strLFQ_IntensityName;
	map<string, string> m_mapExperimentNameAndIBAQIntensityName;
	map<string, string> m_mapExperimentNameAndLFQIntensityName;
	map<string, string> m_mapExperimentNameAndPeptideIntensityName; 
	//string m_strPeptideIntensityName;


	bool mf_LoadProteinFasta(std::ofstream& log, string strProteinFasterFilePath, FastaType fastatype);
	bool mf_LoadPeptides(std::ofstream& log, string strPeptideFilePath);
	bool mf_LoadProteins(std::ofstream& log, string strProteinFilePath);//从proteinGroups.txt载入蛋白质，并且在fasta文件中取得相应的序列
	//string mf_GetPeptidesFlankingRegionFromProteinSequence(string PeptideSequence);
	bool mf_LoadPeptidesEithoutRedunPeptides(std::ofstream& log, string strPeptideFilePath);

	void mf_saveProteins(std::ofstream& log, string proteinspath);
	//void mf_savePeptides(string peptidesFeaturePath);
};




