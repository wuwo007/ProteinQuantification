#pragma once
#include<iostream>
#include<string>
#include<vector>
#include<algorithm>
using namespace std;
class CPredictParam;

class CProteinDigestion
{
public:
	CProteinDigestion(CPredictParam param);
	~CProteinDigestion(void);

	/*-----Get Method----------------*/
	vector <string>  getAllPossiblePeptides(string ProteinSequence); //so get all possible peptides without consider missing cut
	vector <string>  getAllPossiblePeptidesOpt(string ProteinSequence);//optimized method
	bool isInCutPosition(char aChar, string cutPosition);

	int mf_CalculateOptDigestionNumber(string ProteinSequence);
	int mf_CalculateDigestionNumberWithoutMissingCut(string ProteinSeqence);
private:
	string m_strCutPosition;//the position of the protein where enzyme works
	bool m_blr;//cut to the left position or to the right, T-->left,  False-->right
	int m_imaxMissedClevage;//the maxium missed clevage allowed
	int m_iPeptideCount;//total peptide that might be extracted from proteinSequence
	int m_iMinPepLength;//the required minimum length of the peptides, which is used to select peptides, known as L1 in [L1, L2]
	int m_iMaxPepLength;//the required maxium length of the peptides, which is used to select peptides, known as L2 in [L1, L2]

	/*-----------Below are added as the Optimized Protein Digestion added by GZHQ 20160221------*/


};

