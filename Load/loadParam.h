#pragma once
#include"stdafx.h"
class CLoadParam
{
public:
	CLoadParam();
	void mf_setparameters(std::ofstream& log, string parafilepath);

	void mf_showparam();
	string strLoadFileType;
	DataType m_eDataType;
	FastaType m_eFastaType;

	string m_strMaxquantResultPath;
	string m_strPeptideFilePath;
	string m_strProteinFilePath;
	string m_strProteinFasterFilePath;
	string m_strExprimentDesignPath;


	string m_strQuantiResultPath;
	string m_strProteinsPath;
	string m_strpeptidesFeaturePath;

	bool m_bIfIBAQIntensityExist;
	bool m_bIfLFQIntensityExist;
	bool m_bIfExistReverseProtein;
	bool m_bIfExistContanProtein;
	string m_strReverseProteinIDPrefix;  //The prefix of Reverse protein ID;
	string m_strContaminantProteinIDPrefix; //The prefix of Contaminant protein ID;
	string m_strGeneProteinsPath;
	string m_strGeneExpressionQuantificationResultPath;
	//string m_strCutPosition;//the position of the protein where enzyme works
	//bool m_blr;//cut to the left position or to the right, T-->left,  False-->right
	//int m_imaxMissedClevage;//the maxium missed clevage allowed
	//int m_iMinPepLength;//the required minimum length of the peptides, which is used to select peptides, known as L1 in [L1, L2]
	//int m_iMaxPepLength;//the required maxium length of the peptides, which is used to select peptides, known as L2 in [L1, L2]


};
