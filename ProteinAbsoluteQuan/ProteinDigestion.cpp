#include "StdAfx.h"
#include "ProteinDigestion.h"
#include"QuantificationParam.h"

CProteinDigestion::CProteinDigestion(CPredictParam param)
{
	transform(param.m_strCutPosition.begin(), param.m_strCutPosition.end(), param.m_strCutPosition.begin(), toupper);
	this->m_strCutPosition = param.m_strCutPosition;
	this->m_blr = param.m_blr;
	this->m_imaxMissedClevage = param.m_imaxMissedClevage;
	this->m_iPeptideCount = 0;

	//these below 2 lines are added on 19:29 Dec 20, 2013
	this->m_iMinPepLength = param.m_iMinPepLength;
	this->m_iMaxPepLength = param.m_iMaxPepLength;
}
int CProteinDigestion::mf_CalculateOptDigestionNumber(string ProteinSequence)
{
	vector<string> vecPeptidesTemp;
	vecPeptidesTemp = getAllPossiblePeptidesOpt(ProteinSequence);
	return (int)vecPeptidesTemp.size();
}
int CProteinDigestion::mf_CalculateDigestionNumberWithoutMissingCut(string ProteinSequence)
{
	vector<string> vecPeptidesTemp;
	vecPeptidesTemp = getAllPossiblePeptides(ProteinSequence);
	return (int)vecPeptidesTemp.size();
}


CProteinDigestion::~CProteinDigestion(void)
{

}



vector<string> CProteinDigestion::getAllPossiblePeptidesOpt(string proteinSequence)
{
	vector<int> fk;//the positions where enzyme works.
	fk.clear();
	int index;

	vector<string> vecPeptidesTemp;
	if (this->m_blr)//if cut to the left, then we should reverse the protein sequence and cut to the right
	{
		reverse(proteinSequence.begin(), proteinSequence.end());
	}

	/*----Below is cut to the right of the position----*/
	for (index = 0; index < (int)proteinSequence.length(); index++)//all the positions where enzyme works.
	{
		if (this->m_strCutPosition.find(proteinSequence.at(index)) != string::npos)//if found
			fk.push_back(index);
	}

	if (fk.empty())//if the vector is empty(i.e. so far, no cut positions found)
	{
		fk.push_back(proteinSequence.length() - 1);
		vecPeptidesTemp.push_back(proteinSequence);
		return vecPeptidesTemp;
	}
	else if (fk.back() != proteinSequence.length() - 1)//to see if the last AA of protein sequence is in the fk(i.e. the last element of the vector fk)
		fk.push_back(proteinSequence.length() - 1);

	int missedCleavage;
	int finish;//the end index of a new peptide in the protein sequence
	int start;//the start index of a new peptide in the protein sequence
	string strpeptide;
	int pepCount = 0;
	int peptideLength;
	if (((int)fk.size() - 1) > this->m_imaxMissedClevage)
	{
		for (missedCleavage = 0; missedCleavage <= this->m_imaxMissedClevage; missedCleavage++)
		{
			start = 0;
			for (index = 0; index < (int)fk.size() - missedCleavage; index++)//Attention!! size() returns an 'unsigned int', so need to cast it into 'int'
			{
				finish = fk.at(index + missedCleavage);
				if (finish + 1 != proteinSequence.size())
				{
					if ((proteinSequence.at(start + 1) != 'P') && (proteinSequence.at(finish + 1) != 'P'))
					{
						strpeptide = proteinSequence.substr(start, finish - start + 1);

						peptideLength = strpeptide.length();
						if (peptideLength >= this->m_iMinPepLength && peptideLength <= this->m_iMaxPepLength)//if the peptide's length is appropriate
						{
							if (this->m_blr)//if cut to the left, need to reverse the peptide, because right now, it is cut to the right of the position
								reverse(strpeptide.begin(), strpeptide.end());

							vector <string>::iterator iter = std::find(vecPeptidesTemp.begin(), vecPeptidesTemp.end(), strpeptide);
							if (iter == vecPeptidesTemp.end())//if not found before
							{
								vecPeptidesTemp.push_back(strpeptide);
							}
						}//if
					}
					start = fk.at(index) + 1;
				}
				else if (start + 1 < proteinSequence.size())
				{
					if (proteinSequence.at(start + 1) != 'P')
					{
						strpeptide = proteinSequence.substr(start, finish - start + 1);

						peptideLength = strpeptide.length();
						if (peptideLength >= this->m_iMinPepLength && peptideLength <= this->m_iMaxPepLength)//if the peptide's length is appropriate
						{
							if (this->m_blr)//if cut to the left, need to reverse the peptide, because right now, it is cut to the right of the position
								reverse(strpeptide.begin(), strpeptide.end());

							vector <string>::iterator iter = std::find(vecPeptidesTemp.begin(), vecPeptidesTemp.end(), strpeptide);
							if (iter == vecPeptidesTemp.end())//if not found before
							{
								vecPeptidesTemp.push_back(strpeptide);
							}
						}//if
					}
					start = fk.at(index) + 1;
				}

			}//for
		}//for

		if (this->m_blr)//if cut to the left, then we should reverse the protein sequence and cut to the right
		{
			reverse(proteinSequence.begin(), proteinSequence.end());
		}

		//below this block is added @ 14:46 Jan 6, 2014, to handle the protein sequence starting with the AA 'M, meanwhile cut to the right position'
		if (!proteinSequence.empty() && proteinSequence.at(0) == 'M' && (!this->m_blr))
		{
			start = 1;//the AA 'M' at this->proteinSequence.at(0) is easily to get loose and fall apart.
			for (index = 0; index <= this->m_imaxMissedClevage; index++)
			{
				finish = fk.at(index);
				peptideLength = finish - start + 1;

				if (peptideLength >= this->m_iMinPepLength && peptideLength <= this->m_iMaxPepLength)//restraint on peptide
				{
					strpeptide = proteinSequence.substr(start, peptideLength);
					vector <string>::iterator iter = std::find(vecPeptidesTemp.begin(), vecPeptidesTemp.end(), strpeptide);
					if (iter == vecPeptidesTemp.end())//if not found before
					{
						vecPeptidesTemp.push_back(strpeptide);
					}
				}
			}//for
		}//if branch
	}
	else
	{
		int iMaxMissedCleavageTemp = fk.size() - 1;
		for (missedCleavage = 0; missedCleavage <= iMaxMissedCleavageTemp; missedCleavage++)
		{
			start = 0;
			for (index = 0; index < (int)fk.size() - missedCleavage; index++)//Attention!! size() returns an 'unsigned int', so need to cast it into 'int'
			{
				finish = fk.at(index + missedCleavage);
				if (finish + 1 != proteinSequence.size())
				{
					if ((proteinSequence.at(start + 1) != 'P') && (proteinSequence.at(finish + 1) != 'P'))
					{
						strpeptide = proteinSequence.substr(start, finish - start + 1);

						peptideLength = strpeptide.length();
						if (peptideLength >= this->m_iMinPepLength && peptideLength <= this->m_iMaxPepLength)//if the peptide's length is appropriate
						{
							if (this->m_blr)//if cut to the left, need to reverse the peptide, because right now, it is cut to the right of the position
								reverse(strpeptide.begin(), strpeptide.end());

							vector <string>::iterator iter = std::find(vecPeptidesTemp.begin(), vecPeptidesTemp.end(), strpeptide);
							if (iter == vecPeptidesTemp.end())//if not found before
							{
								vecPeptidesTemp.push_back(strpeptide);
							}
						}//if
					}
					start = fk.at(index) + 1;
				}
				else if (start + 1<proteinSequence.size())
				{
					if (proteinSequence.at(start + 1) != 'P')
					{
						strpeptide = proteinSequence.substr(start, finish - start + 1);

						peptideLength = strpeptide.length();
						if (peptideLength >= this->m_iMinPepLength && peptideLength <= this->m_iMaxPepLength)//if the peptide's length is appropriate
						{
							if (this->m_blr)//if cut to the left, need to reverse the peptide, because right now, it is cut to the right of the position
								reverse(strpeptide.begin(), strpeptide.end());

							vector <string>::iterator iter = std::find(vecPeptidesTemp.begin(), vecPeptidesTemp.end(), strpeptide);
							if (iter == vecPeptidesTemp.end())//if not found before
							{
								vecPeptidesTemp.push_back(strpeptide);
							}
						}//if
					}
					start = fk.at(index) + 1;
				}

			}//for
		}//for

		if (this->m_blr)//if cut to the left, then we should reverse the protein sequence and cut to the right
		{
			reverse(proteinSequence.begin(), proteinSequence.end());
		}

		//below this block is added @ 14:46 Jan 6, 2014, to handle the protein sequence starting with the AA 'M, meanwhile cut to the right position'
		if (!proteinSequence.empty() && proteinSequence.at(0) == 'M' && (!this->m_blr))
		{
			start = 1;//the AA 'M' at this->proteinSequence.at(0) is easily to get loose and fall apart.
			for (index = 0; index <= iMaxMissedCleavageTemp; index++)
			{
				finish = fk.at(index);
				peptideLength = finish - start + 1;

				if (peptideLength >= this->m_iMinPepLength && peptideLength <= this->m_iMaxPepLength)//restraint on peptide
				{
					strpeptide = proteinSequence.substr(start, peptideLength);
					vector <string>::iterator iter = std::find(vecPeptidesTemp.begin(), vecPeptidesTemp.end(), strpeptide);
					if (iter == vecPeptidesTemp.end())//if not found before
					{
						vecPeptidesTemp.push_back(strpeptide);
					}
				}
			}//for
		}//if branch
	}


	return vecPeptidesTemp;
	//Commented @ 3-3-2014
}


/*
do not consider missing cut
constrain peptide lenth 
*/
vector <string>  CProteinDigestion::getAllPossiblePeptides(string proteinSequence)
{
	vector<int> fk;//the positions where enzyme works.
	fk.clear();
	int index;

	vector<string> vecPeptidesTemp;
	if (this->m_blr)//if cut to the left, then we should reverse the protein sequence and cut to the right
	{
		reverse(proteinSequence.begin(), proteinSequence.end());
	}
	
	/*----Below is cut to the right of the position----*/
	for (index = 0; index < (int)proteinSequence.length(); index++)//all the positions where enzyme works.
	{
		if (this->m_strCutPosition.find(proteinSequence.at(index)) != string::npos)//if found
			fk.push_back(index);
	}

	if (fk.empty())//if the vector is empty(i.e. so far, no cut positions found)
	{
		
		vecPeptidesTemp.push_back(proteinSequence);
		return vecPeptidesTemp;
	}
	else if (fk.back() != proteinSequence.length() - 1)//to see if the last AA of protein sequence is in the fk(i.e. the last element of the vector fk)
		fk.push_back(proteinSequence.length() - 1);

	//Below this block is added @ 15:14 Jan 6, 2013. To prevent malicious attack.
	//if (this->m_imaxMissedClevage > ((int)fk.size() - 1))
	//{
	//	cout << "Max-Missed-Cleavage Is Not Appropriate For protein sequence "<<proteinSequence << endl;

	//	system("pause");
	//}
	string strpeptide;
	int finish;//the end index of a new peptide in the protein sequence
	int start;//the start index of a new peptide in the protein sequence
	start = 0;
	int peptideLength;
	for (index = 0; index < (int)fk.size(); index++)//Attention!! size() returns an 'unsigned int', so need to cast it into 'int'
	{
		finish = fk.at(index );
		if (finish + 1 != proteinSequence.size())
		{
			if ((proteinSequence.at(start + 1) != 'P') && (proteinSequence.at(finish + 1) != 'P'))
			{
				strpeptide = proteinSequence.substr(start, finish - start + 1);

				peptideLength = strpeptide.length();
				if (peptideLength >= this->m_iMinPepLength && peptideLength <= this->m_iMaxPepLength)//if the peptide's length is appropriate
				{
					if (this->m_blr)//if cut to the left, need to reverse the peptide, because right now, it is cut to the right of the position
						reverse(strpeptide.begin(), strpeptide.end());

					vector <string>::iterator iter = std::find(vecPeptidesTemp.begin(), vecPeptidesTemp.end(), strpeptide);
					if (iter == vecPeptidesTemp.end())//if not found before
					{
						vecPeptidesTemp.push_back(strpeptide);
					}
				}//if
			}
			start = fk.at(index) + 1;
		}
		else if (start + 1 < proteinSequence.size())
		{
			if (proteinSequence.at(start + 1) != 'P')
			{
				strpeptide = proteinSequence.substr(start, finish - start + 1);

				peptideLength = strpeptide.length();
				if (peptideLength >= this->m_iMinPepLength && peptideLength <= this->m_iMaxPepLength)//if the peptide's length is appropriate
				{
					if (this->m_blr)//if cut to the left, need to reverse the peptide, because right now, it is cut to the right of the position
						reverse(strpeptide.begin(), strpeptide.end());

					vector <string>::iterator iter = std::find(vecPeptidesTemp.begin(), vecPeptidesTemp.end(), strpeptide);
					if (iter == vecPeptidesTemp.end())//if not found before
					{
						vecPeptidesTemp.push_back(strpeptide);
					}
				}//if
			}
			start = fk.at(index) + 1;
		}
	}

	return vecPeptidesTemp;
}
bool CProteinDigestion::isInCutPosition(char aChar, string cutPosition)
{
	if (cutPosition.find(aChar) != string::npos)//if found
		return true;
	return false;
}




