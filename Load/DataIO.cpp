#include"stdafx.h"
#include"DataIO.h"
#include <windows.h>
#include<algorithm>
#include<math.h>
#include <sys/stat.h>

#ifndef S_ISDIR
#define S_ISDIR(mode)  (((mode) & S_IFMT) == S_IFDIR)
#endif

#ifndef S_ISREG
#define S_ISREG(mode)  (((mode) & S_IFMT) == S_IFREG)
#endif
CLoadIO::CLoadIO()
{
	m_bIfIBAQIntensityExist = false;
	m_bIfLFQIntensityExist = false;
}

string CLoadIO::mf_GetPeptidesSequenceFromID(std::ofstream& log, int ID)
{
	map<int, string>::iterator Iter;
	Iter = m_Peptides.m_mapPeptideIDAndSequence.find(ID);
	if (Iter != m_Peptides.m_mapPeptideIDAndSequence.end())
		return Iter->second;
	else
	{
		cout << "Can not find peptide " << ID << " in  filtered peptides!" << endl;
		log << "Can not find peptide " << ID << " in  filtered peptides!" << endl;
		return "null";
	}


}
double CLoadIO::mf_GetPeptidesIntensityFromID(std::ofstream& log, int ID, int iExperimentIndex)
{
	map<int, vector<double>>::iterator Iter;
	Iter = m_Peptides.m_mapPeptideIDAndIntensitys.find(ID);
	if (Iter != m_Peptides.m_mapPeptideIDAndIntensitys.end())
		return Iter->second.at(iExperimentIndex);
	else
	{
		cout << "Can not find peptide's intensity " << ID << endl;
		log << "Can not find peptide's intensity " << ID << endl;
		return 0;
	}

}

bool CLoadIO::mf_LoadProteinIDAndSequence(std::ofstream& log, map<string, string> &IdAndSequence, string strProteinFasterFilePath, FastaType fastatype)
{
	char Buffer[BUFFERLENGTH];
	FILE * pFile;
	pFile = fopen(strProteinFasterFilePath.c_str(), "r");
	if (pFile == NULL)
	{
		cout << "Can not open " << strProteinFasterFilePath << endl;
		log << "Error:\tCan not open " << strProteinFasterFilePath << endl;
		exit(0);
	}
	char *pstr;
	char *pstr1;
	string strIDTemp;
	string strSequenceTemp;
	string strTemp;
	fgets(Buffer, BUFFERLENGTH, pFile);
	while ((Buffer[0] == '\n') && (!feof(pFile)))  //allow empty row;
	{
		fgets(Buffer, BUFFERLENGTH, pFile);
	}
	
	std::match_results<std::string::const_iterator> Matchresult;
	bool valid;
	string strFirstRow;
	pstr = Buffer;
	if (Buffer[0] == '>')  //get the protein ID according to the regular expression fastatype;
	{
		strFirstRow = Buffer;
		valid = std::regex_search(strFirstRow, Matchresult, fastatype);
		if (valid == true)
		{
			strIDTemp = Matchresult[1];
		}
		else
		{
			cout << "Can not sparse the fasta file by the regular expression"<< endl;
			log << "Error:\tCan not sparse the fasta file by the regular expression" << endl;
			exit(0);
		}
	}

	fgets(Buffer, BUFFERLENGTH, pFile);

	while (Buffer[0] != '>')  // get the sequence of the protein id;
	{

		strTemp = Buffer;
		strSequenceTemp = strSequenceTemp + strTemp.substr(0, strTemp.size() - 1);
		fgets(Buffer, BUFFERLENGTH, pFile);
	}
	while (!feof(pFile))
	{
		if (Buffer[0] == '\0')   //allow the blank lines
			continue;
		pstr = Buffer;

		if (Buffer[0] == '>')  //get the protein ID according to the regular expression fastatype;
		{
			IdAndSequence.insert(pair<string, string>(strIDTemp, strSequenceTemp));
			strSequenceTemp.clear();
			strFirstRow = Buffer;
			valid = std::regex_search(strFirstRow, Matchresult, fastatype);
			if (valid == true)
			{
				strIDTemp = Matchresult[1];
			}
			else
			{
				cout << "Can not sparse the fasta file by the regular expression" << endl;
				log << "Error:\tCan not sparse the fasta file by the regular expression" << endl;
				exit(0);
			}

			fgets(Buffer, BUFFERLENGTH, pFile);
		}
		else // get the sequence of the protein id;
		{
			strTemp.clear();
			strTemp = pstr;
			strSequenceTemp = strSequenceTemp + strTemp.substr(0, strTemp.size() - 1);
			fgets(Buffer, BUFFERLENGTH, pFile);
		}
	}
	IdAndSequence.insert(pair<string, string>(strIDTemp, strSequenceTemp));

	fclose(pFile);
	return 1;
}

double CLoadIO::mf_GetMedianIntensity(std::ofstream& log, vector<double> PeptidesIntensity)
{
	sort(PeptidesIntensity.begin(), PeptidesIntensity.end());
	int len = PeptidesIntensity.size();
	if (len == 0)             //added by gzhq 20150529
		return 0.0;
	if (len % 2 == 0)
		return (PeptidesIntensity.at(len / 2) + PeptidesIntensity.at((len - 2) / 2)) / 2;
	else
		return PeptidesIntensity.at(floor(len / 2));
}
void CLoadIO::mf_saveProteins(std::ofstream& log, string proteinspath)
{
	ofstream ofile(proteinspath.c_str());
	if (!ofile)
	{
		cout << "can not open " << proteinspath << endl;
		log << "can not open " << proteinspath << endl;
	}
	else
	{
		cout << "\tsave Proteinpeptides in " << proteinspath << endl;
		log << "\tsave Proteinpeptides in " << proteinspath << endl;
	}
	ofile << "ProteinID\tMajority protein IDs\tPeptidesNumber\tvPeptidesSequence\tvPeptidesIntensitys\tbIfIBAQExist\tdMaxquantIBAQ\tbIfMaxquantLFQExist\tdMaxquantLFQ\n";

	vector<CProtein>::iterator  ProteinIter;
	vector<string>::iterator vPeptideIter;
	vector<double>::iterator vPeptideIntensityIter;
	for (ProteinIter = m_vecProteins.begin(); ProteinIter != m_vecProteins.end(); ProteinIter++)
	{
		ofile << ProteinIter->m_strProteinID << "\t" << ProteinIter->m_strProteinFullName << "\t";
		ofile << ProteinIter->m_iPeptidesNumber << "\t";
		for (vPeptideIter = ProteinIter->m_vPeptidesSequences.begin(); vPeptideIter != ProteinIter->m_vPeptidesSequences.end() - 1; vPeptideIter++)
			ofile << *vPeptideIter << ";";
		ofile << *vPeptideIter << "\t";
		ofile << *vPeptideIntensityIter << "\t";

		if (m_bIfIBAQIntensityExist)
		{
			ofile << 1 << "\t";
		}
		else
		{
			ofile << 0 << "\t";
			ofile << 0 << "\t";
		}
		if (m_bIfLFQIntensityExist)
		{
			ofile << 1 << "\t";
		}
		else
		{
			ofile << 0 << "\t";
			ofile << 0 << "\n";
		}

	}
	ofile.close();
}

void CLoadIO::mf_GetAttributesFromFirstRow(char * pstr, map<string, int>& mapAttributesAndColumns)
{
	char *pstr1;
	int icolumns = 0;

	pstr1 = strstr(pstr, "\t");
	while (pstr1 != NULL)
	{
		if (pstr1 != NULL)
		{
			*pstr1 = '\0';
		}
		else
		{
			cout << "Error:\tThe format of proteinGroups.txt  is wrong.\n";
			exit(0);
		}
		mapAttributesAndColumns.insert(pair<string, int>(pstr, icolumns));
		icolumns++;
		pstr = pstr1 + 1;
		pstr1 = strstr(pstr, "\t");
	}
	pstr1 = strstr(pstr, "\n");
	if (pstr1 != NULL)
	{
		*pstr1 = '\0';
	}
	else
	{
		cout << "Error:\tThe format of proteinGroups.txt  is wrong.\n";
		exit(0);
	}
	mapAttributesAndColumns.insert(pair<string, int>(pstr, icolumns));
}
CLoadMaxquantIO::CLoadMaxquantIO(std::ofstream& log, const CLoadParam loadParam)
{
	m_NumberOfPeptides = 0;
	m_bIfExistReverseProtein = loadParam.m_bIfExistReverseProtein;
	m_bIfExistContanProtein = loadParam.m_bIfExistContanProtein;
	m_strReversePrefix = loadParam.m_strReverseProteinIDPrefix;
	m_strContaminantfix = loadParam.m_strContaminantProteinIDPrefix;
	mf_GetAttributesName(log,loadParam.m_strExprimentDesignPath);
}


/*mf_GetAttributesName
Because some attribute name change with the experiment design,
so we need determine the attribute name according to the experiment design file.
*/
void CLoadMaxquantIO::mf_GetAttributesName(std::ofstream& log, string ExperimentDesignPath)
{
	// get the string between the second tab and the line break in the rows except the first row；
	map<string, int> mapExperiments;
	map<string, int>::iterator mapExperimentsIter;
	ifstream fin(ExperimentDesignPath.c_str());
	if (!fin)
	{
		cout << "Can not open file " << ExperimentDesignPath << endl;
		log << "Error:\tCan not open file " << ExperimentDesignPath << endl;
		exit(2);
	}
	string strTemp1, strTemp2;
	int iTemp1 = 0, iTemp2 = 0;
	m_iNumberOfExperiments = 0;
	getline(fin, strTemp1);//jump the first row；	
	strTemp2 = "";
	mapExperiments.insert(pair<string, int>(strTemp2, 0)); //add by gzhq 20160225
	while (getline(fin, strTemp1))
	{
		iTemp1 = strTemp1.find("\t");
		iTemp1 = strTemp1.find("\t", iTemp1 + 1);
		iTemp2 = strTemp1.find("\n", iTemp1 + 1);
		strTemp2 = strTemp1.substr(iTemp1 + 1, iTemp2 - iTemp1 - 1);
		mapExperiments.insert(pair<string, int>(strTemp2, 0));
		
	}
	m_iNumberOfExperiments = mapExperiments.size();

	string strPeptideIntensityName;
	string strIBAQ_IntensityName;
	string strLFQ_IntensityName;

	//if (mapExperiments.size() > 0)
	//{
		mapExperimentsIter = mapExperiments.begin();
		for (; mapExperimentsIter != mapExperiments.end(); mapExperimentsIter++)
		{
			strPeptideIntensityName = "Intensity " + mapExperimentsIter->first;
			m_mapExperimentNameAndPeptideIntensityName.insert(pair<string, string>(mapExperimentsIter->first, strPeptideIntensityName));
			strIBAQ_IntensityName = "iBAQ " + mapExperimentsIter->first;
			m_mapExperimentNameAndIBAQIntensityName.insert(pair<string, string>(mapExperimentsIter->first, strIBAQ_IntensityName));
			if (mapExperimentsIter->first != "")
			{
				strLFQ_IntensityName = "LFQ intensity " + mapExperimentsIter->first;
				m_mapExperimentNameAndLFQIntensityName.insert(pair<string, string>(mapExperimentsIter->first, strLFQ_IntensityName));
			}
		}						   

	//}
/*	else
	{
		strPeptideIntensityName = "Intensity";
		strIBAQ_IntensityName = "iBAQ";
		strLFQ_IntensityName = "LFQ intensity";

	}*/	
	fin.close();
}
/*mf_GetPeptideAttributecolumn
    determine the column number of peptide attributes which we need
	 "Sequence"，"id ", "Proteins",strPeptideIntensityName
*/
int CLoadMaxquantIO::mf_GetPeptideAttributecolumn(std::ofstream& log, std::string AttributeName)
{
	if (AttributeName == "Sequence")
		return 0;
	else if (AttributeName == "id")
		return 1;
	else if (AttributeName == "Proteins")  //change by gzhq 20150602
		return 2;
	int icolumn = 2;
	map<string, string>::iterator mapExperimentAndPeptideIntensityIter;
	mapExperimentAndPeptideIntensityIter = m_mapExperimentNameAndPeptideIntensityName.begin();
	for (; mapExperimentAndPeptideIntensityIter != m_mapExperimentNameAndPeptideIntensityName.end(); mapExperimentAndPeptideIntensityIter++)
	{
		if (AttributeName == mapExperimentAndPeptideIntensityIter->second)
		{
			icolumn++;
			return icolumn;
		}
	}
	 
	cout<<"Can not find the column: "<<AttributeName<<endl;
	log << "Can not find the column: " << AttributeName << endl;
	return -1;

}


/*mf_LoadPeptides
extract peptides information from peptides.txt ,and save it in m_Peptides.
informations include:
"Sequence", "id", "Proteins", "Leading razor protein", peptide intensity;
deleting constraint:
1, column " Proteins"  has more than one protein;
2,
*/

bool CLoadMaxquantIO::mf_LoadPeptides(std::ofstream& log, string strPeptideFilePath)
{
	//firstly, link up the proteins'sequences, save the start location of every protein
	string strProteinsSequenceTemp;

	map<string, string>::iterator ProteinIDAndSequenceIter=m_mapProteinIDAndSequence.begin();
	int iPositionTemp = 0, nextPosition = 0;
	
	for (; ProteinIDAndSequenceIter != m_mapProteinIDAndSequence.end(); ProteinIDAndSequenceIter++)
	{
		strProteinsSequenceTemp = strProteinsSequenceTemp +"?"+ ProteinIDAndSequenceIter->second;
	}

	int ibeginTemp = 0,icutPosition = 0;
	string ProteinId1;

	cout << "\tLoading Peptides\n";
	log << "\tLoading Peptides\n";
	FILE * pFile;
	pFile = fopen(strPeptideFilePath.c_str(), "r");
	if (pFile == NULL)
	{
		cout << "Can not open " << strPeptideFilePath << endl;
		log << "Error:\tCan not open " << strPeptideFilePath << endl;
		exit(0);

	}
	char Buffer[BUFFERLENGTH];
	char *pstr;
	char *pstr1;
	fgets(Buffer, BUFFERLENGTH, pFile);    

	int iSequencecolumnNum = 0, iIDcolumnNum = 0, iProteinsNamecolumnNum = 0, iIntensitycolumnNum = 0;
	bool bIfFindSequence = false, bIfFindID = false, bIfFindProteinNames = false, bIffindIntensity = false;
	pstr = Buffer;
	int icolumns = 0, count = 0;
	pstr1 = strstr(pstr, "\t");
	while (count<4&&pstr1!=NULL)
	{
		if (pstr1 != NULL)
		{
			*pstr1 = '\0';
		}
		else
		{
			log << "Error:\tThe format of peptides.txt  is wrong.\n";
			exit(0);
		}
		switch (mf_GetPeptideAttributecolumn(log,pstr))
		{
		case 0:
			iSequencecolumnNum = icolumns;
			bIfFindSequence = true;
			count++;
			break;
		case 1:
			iIDcolumnNum = icolumns;
			bIfFindID = true;
			count++;
			break;
		case 2:
			iProteinsNamecolumnNum = icolumns;
			bIfFindProteinNames = true;
			count++;
			break;
		case 3:
			iIntensitycolumnNum = icolumns;
			bIffindIntensity = true;
			count++;
			break;
		}
		icolumns++;
		pstr = pstr1 + 1;
		pstr1 = strstr(pstr, "\t");
	}
	if (count < 4)  
	{
		if (!bIfFindID)
		{
			cout << "Can not find the column \"id\" in the peptides.txt" << endl;
			log << "Error:\tCan not find the column \"id\" in the peptides.txt" << endl;
			exit(0);
		}
		if (!bIffindIntensity)
		{
			cout << "Can not find the column  m_strPeptideIntensityName in the peptides.txt" << endl;
			log << "Error:\tCan not find the column  m_strPeptideIntensityName in the peptides.txt" << endl;
			exit(0);
		}
		if (!bIfFindProteinNames)
		{
			cout << "Can not find the column \"Proteins\" in the peptides.txt" << endl;
			log << "Error:\tCan not find the column \"Proteins\" in the peptides.txt" << endl;
			exit(0);
		}
		if (!bIfFindSequence)
		{
			cout << "Can not find the column \"Sequence\" in the peptides.txt" << endl;
			log << "Error:\tCan not find the column \"Sequence\" in the peptides.txt" << endl;
			exit(0);
		}
	}
	icolumns = 0;
	count = 0;

	bool bIfDelete;
	string strSequenceTemp;
	double dIntensityTemp;
	int iIDTemp;
	map<string, string>::iterator mapProteinIDAndSequenceTempIter; 
	int iNormalProteinsNumber;  
	bool b_ifExistReverseProtein = true;
	bool b_ifExistcontaminantProtein = true;
	//string strReversePrefix="REV_";
	//string strcontaminantfix = "CON_";
	string strProteinNameTemp;
	char * pstr2; 
	fgets(Buffer, BUFFERLENGTH, pFile);
	pstr = Buffer;
	while (!feof(pFile))
	{
		bIfDelete = false;
		count = 0;
		icolumns = 0;
		
		while (count < 4)
		{
			pstr1 = strstr(pstr, "\t");
			if (pstr1 != NULL)
			{
				*pstr1 = '\0';
			}
			else
			{
				log << "Error:\tThe format of peptides.txt  is wrong.\n";
				exit(0);
			}
			if (icolumns == iSequencecolumnNum)
			{
				strSequenceTemp = pstr;
				count++;
			}
			if (icolumns == iIDcolumnNum)
			{
				iIDTemp = atoi(pstr);
				count++;
			}
			if (icolumns == iProteinsNamecolumnNum)
			{
				if (b_ifExistReverseProtein)
					if (strstr(pstr, m_strReversePrefix.c_str()) != NULL || (pstr == ""))
						bIfDelete = true;                      // delete the peptide, if its corresponding protein group contain reverse protein
				
				if (!b_ifExistcontaminantProtein)
				{  //do not contain contamination proteins
					if (strstr(pstr, ";") != NULL) //remove the shared peptides added by ChangCheng 20150527
						bIfDelete = true;
				}
				else
				{//may contain contamination proteins
					iNormalProteinsNumber = 0;
					pstr2 = strstr(pstr, ";");
					while (pstr2 != NULL)
					{

						*pstr2 = '\0';
						strProteinNameTemp = pstr;
						strProteinNameTemp = strProteinNameTemp.substr(0, m_strContaminantfix.length()); //change the constant 4 to the length of strContaminatefix by CC 20150727
						if (strProteinNameTemp != m_strContaminantfix)
							iNormalProteinsNumber++;
						pstr = pstr2 + 1;
						pstr2 = strstr(pstr, ";");
					}
					strProteinNameTemp = pstr;
					strProteinNameTemp = strProteinNameTemp.substr(0, m_strContaminantfix.length());
					if ((strProteinNameTemp != m_strContaminantfix) && (strProteinNameTemp != ""))
						iNormalProteinsNumber++;

					if (iNormalProteinsNumber != 1)
						bIfDelete = true;
				}
				
				count++;
			}
			if (icolumns == iIntensitycolumnNum)
			{
				dIntensityTemp = atof(pstr);
				count++;
			}
			icolumns++;
			pstr = pstr1 + 1;
		}

		if (!bIfDelete)
		{
			ibeginTemp = 0;
			iPositionTemp = strProteinsSequenceTemp.find(strSequenceTemp, ibeginTemp);
			if (iPositionTemp == strProteinsSequenceTemp.npos)
			{
				cerr << "Error:\tThere is not peptide " << strSequenceTemp << " in Protein fasta" << endl;
				exit(0);
			}
			else
			{
				icutPosition = strProteinsSequenceTemp.find("?", iPositionTemp + 1);
				if (icutPosition != strProteinsSequenceTemp.npos)
				{
					iPositionTemp = strProteinsSequenceTemp.find(strSequenceTemp, icutPosition + 1);
					if (iPositionTemp != strProteinsSequenceTemp.npos)
						bIfDelete = true;
				}
			}
		}
	 
		if (!bIfDelete)
		{
			m_Peptides.m_mapPeptideIDAndSequence.insert(pair<int, string>(iIDTemp, strSequenceTemp));	
		}
		fgets(Buffer, BUFFERLENGTH, pFile);
		pstr = Buffer;
	}
	fclose(pFile);
	return 1;
}

/*mf_LoadPeptidesEithoutRedunPeptides
extract peptides information from peptides.txt ,and save it in m_Peptides. but without redundancy removal 
informations include:
"Sequence", "id", "Proteins", "Leading razor protein", peptide intensity;
deleting constraint:
1, column " Proteins"  has more than one protein;
2,
*/
bool CLoadMaxquantIO::mf_LoadPeptidesEithoutRedunPeptides(std::ofstream& log, string strPeptideFilePath)
{
	cout << "\tLoading Peptides\n";
	log << "\tLoading Peptides\n";

	// open the file
	FILE * pFile;
	pFile = fopen(strPeptideFilePath.c_str(), "r");
	if (pFile == NULL)
	{
		cout << "Can not open " << strPeptideFilePath << endl;
		log << "Error:\tCan not open " << strPeptideFilePath << endl;
		exit(0);
	}


	char Buffer[BUFFERLENGTH];
	char *pstr;
	char *pstr1;
	//get the columns of sequence、ID、Proteins、Intensity UPS2_yeast by analysising the first row;
	fgets(Buffer, BUFFERLENGTH, pFile);    

	int iSequencecolumnNum = 0, iIDcolumnNum = 0, iProteinsNamecolumnNum = 0;
	vector<int> veciIntensitycolumnNum;
	bool bIfFindSequence = false, bIfFindID = false, bIfFindProteinNames = false;
	vector<bool> vecbIffindIntensity;
	int icolumns = 0, count = 0;
	int iNumberOfcolumns = 3 + m_iNumberOfExperiments;
	map<string, int> mapAttrtibuteAndcolumns;
	map<string, int>::iterator mapAttrtibuteAndcolumnsIter;

	pstr = Buffer;
	mf_GetAttributesFromFirstRow(pstr, mapAttrtibuteAndcolumns);
	
	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Sequence");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iSequencecolumnNum = mapAttrtibuteAndcolumnsIter->second;
		bIfFindSequence = true;
		count++;
	}
	else
	{
		bIfFindSequence = false;
	}		

	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("id");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iIDcolumnNum = mapAttrtibuteAndcolumnsIter->second;
		bIfFindID = true;
		count++;
	}
	else
	{
		bIfFindID = false;
	}

	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Proteins");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iProteinsNamecolumnNum= mapAttrtibuteAndcolumnsIter->second;
		bIfFindProteinNames = true;
		count++;
	}
	else
	{
		bIfFindProteinNames = false;
	}

	map<string, string>::iterator mapExperimentAndPeptideIntensityIter;
	mapExperimentAndPeptideIntensityIter = m_mapExperimentNameAndPeptideIntensityName.begin();
	for (; mapExperimentAndPeptideIntensityIter != m_mapExperimentNameAndPeptideIntensityName.end(); mapExperimentAndPeptideIntensityIter++)
	{
		mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find(strim(mapExperimentAndPeptideIntensityIter->second));
		if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
		{
			veciIntensitycolumnNum.push_back(mapAttrtibuteAndcolumnsIter->second);
			vecbIffindIntensity.push_back(true);
			count++;
		}
		else
		{
			vecbIffindIntensity.push_back(false);
		}
	}

	if (count < iNumberOfcolumns)
	{
		if (!bIfFindID)
		{
			cout << "Can not find the column \"id\" in the peptides.txt" << endl;
			log << "Error:\tCan not find the column \"id\" in the peptides.txt" << endl;
			exit(0);
		}
		if (!bIfFindProteinNames)
		{
			cout << "Can not find the column \"Proteins\" in the peptides.txt" << endl;
			log << "Error:\tCan not find the column \"Proteins\" in the peptides.txt" << endl;
			exit(0);
		}
		if (!bIfFindSequence)
		{
			cout << "Can not find the column \"Sequence\" in the peptides.txt" << endl;
			log << "Error:\tCan not find the column \"Sequence\" in the peptides.txt" << endl;
			exit(0);
		}
		for (int i = 0; i < vecbIffindIntensity.size(); i++)
		{
			if (!vecbIffindIntensity.at(i))
			{
				cout << "Can not find the column PeptideIntensityName in the peptides.txt" << endl;
				log << "Error:\tCan not find the column PeptideIntensityName in the peptides.txt" << endl;
				exit(0);

			}
		}

	}

	icolumns = 0;
	count = 0;
	bool bIfDelete;
	string strSequenceTemp;
	vector<double> vecdIntensityTemp;
	int iIDTemp;
	map<string, string>::iterator mapProteinIDAndSequenceTempIter;
	int iCountTemp;
	int iNormalProteinsNumber;  
	bool b_ifExistReverseProtein = true;
	bool b_ifExistcontaminantProtein = true;
	//string strReversePrefix = "REV_";
	//string strcontaminantfix = "CON_";
	string strProteinNameTemp;
	char * pstr2; 
	fgets(Buffer, BUFFERLENGTH, pFile);
	pstr = Buffer;
	while (!feof(pFile))
	{	
		bIfDelete = false;
		count = 0;
		icolumns = 0;
		vecdIntensityTemp.clear();
		while (count < iNumberOfcolumns)
		{
			pstr1 = strstr(pstr, "\t");
			if (pstr1 != NULL)
			{
				*pstr1 = '\0';
			}
			else
			{
				log << "Error:\tThe format of peptides.txt  is wrong.\n";
				exit(0);
			}
			if (icolumns == iSequencecolumnNum)
			{
				strSequenceTemp = pstr;
				count++;
			}
			if (icolumns == iIDcolumnNum)
			{
				iIDTemp = atoi(pstr);
				count++;
			}
			if (icolumns == iProteinsNamecolumnNum)
			{
				if (b_ifExistReverseProtein)
					if (strstr(pstr, m_strReversePrefix.c_str()) != NULL || (pstr == ""))
						bIfDelete = true;

				//  change by gzhq 20150909   // delete the peptide, if its corresponding protein group contain reverse protein
				if (!b_ifExistcontaminantProtein)
				{  
					if (strstr(pstr, ";") != NULL) //remove the shared peptides added by ChangCheng 20150527
						bIfDelete = true;
				}
				else
				{
					iNormalProteinsNumber = 0;
					pstr2 = strstr(pstr, ";");
					while (pstr2 != NULL)
					{
						*pstr2 = '\0';
						strProteinNameTemp = pstr;
						strProteinNameTemp = strProteinNameTemp.substr(0, m_strContaminantfix.length()); //前缀 //change the constant 4 to the length of strContaminatefix by CC 20150727
						if (strProteinNameTemp != m_strContaminantfix)
							iNormalProteinsNumber++;
						pstr = pstr2 + 1;
						pstr2 = strstr(pstr, ";");
					}
					strProteinNameTemp = pstr;
					strProteinNameTemp = strProteinNameTemp.substr(0, m_strContaminantfix.length());
					if (strProteinNameTemp != m_strContaminantfix)
						iNormalProteinsNumber++;

					if (iNormalProteinsNumber != 1)
						bIfDelete = true;
				}
				count++;
			}
			for (int i = 0; i < veciIntensitycolumnNum.size(); i++)
			{
				if (icolumns == veciIntensitycolumnNum.at(i))
				{
					vecdIntensityTemp.push_back(atof(pstr));
					count++;
				}
			}
			icolumns++;
			pstr = pstr1 + 1;
		}

		if (!bIfDelete)
		{
			m_Peptides.m_mapPeptideIDAndIntensitys.insert(pair<int, vector<double>>(iIDTemp, vecdIntensityTemp));
			m_Peptides.m_mapPeptideIDAndSequence.insert(pair<int, string>(iIDTemp, strSequenceTemp));
			m_Peptides.m_mapPeptideSequenceAndIntensitys.insert(pair<string, vector<double>>(strSequenceTemp, vecdIntensityTemp));
		}
		fgets(Buffer, BUFFERLENGTH, pFile);
		pstr = Buffer;
	}
	fclose(pFile);
	cout << "\tLoaded " << m_Peptides.m_mapPeptideIDAndSequence.size() << " Peptides!\n\n";
	log << "\tLoaded " << m_Peptides.m_mapPeptideIDAndSequence.size() << " Peptides!\n\n";
	return 1;
}


bool CLoadMaxquantIO::mf_LoadProteinFasta(std::ofstream& log, string strProteinFasterFilePath, FastaType fastatype)
{
	if (mf_LoadProteinIDAndSequence(log,m_mapProteinIDAndSequence, strProteinFasterFilePath, fastatype))
		return true;
	else
		return false;
}

/*mf_LoadProteins
extract proteins information from proteinGroups.txt ,and save it in m_vecProteins.
informations include:
"Majority protein IDs", "Unique peptides", "iBAQ UPS2_...", "LFQ intensity UPS2_...",
"Reverse", "Potential contaminant", "Peptide IDs", "Peptide is razor",
deleting constraint:
1, column " Reverse" contain "+";
3, proteinID is P07339 and P08311;
4, "Peptide is razor" do not contain true;
5, peptide intensity equal or less than 0, delete it; 
   if the number of protein's peptides don't  greater than 0, delete it;
*/
bool CLoadMaxquantIO::mf_LoadProteins(std::ofstream& log, string strProteinFilePath)
{
	cout << "\tLoading Protein!\n";
	log << "\tLoading Protein!\n";

	//open the file
	FILE * pFile;
	pFile = fopen(strProteinFilePath.c_str(), "r");
	if (pFile == NULL)
	{
		cout << "Can not open " << strProteinFilePath << endl;
		log << "Error:\tCan not open " << strProteinFilePath << endl;
		exit(0);
	}

	char Buffer[BUFFERLENGTH];
	char *pstr;
	char *pstr1;
	map<string, int> mapAttrtibuteAndcolumns;
	map<string, int>::iterator mapAttrtibuteAndcolumnsIter;
	int icolumns = 0;
	int count = 0;

	int iMajorProteinIDcolumnNum = 0;
	vector<int> veciBAQcolumnNum;
	vector<int> veciLFQintensitycolumnNum;
	int iReversecolumnNum = 0;
	int iPoteinalContaminantcolumnNum = 0;
	int iPeptideIDscolumnNum = 0;
	int iPeptideRazorcolumnNum = 0;

	bool bIfFindProteinID = false;
	vector<bool> vecbIfFindIBAQ;
	vector<bool> vecbIfFindLFQIntensity;
	bool bIFFindReverse = false;
	bool bIfFindPoteinalContaminant = false;
	bool bIfFindPeptideIDs = false;
	bool bIfFindPeptideRazor = false;

	int iAttributeNumbers = 5+2*m_iNumberOfExperiments-1;// The number of Attributes want to get from proteingroups.txt Add by gzhq20150528; There is not "LFQ intensity";


	fgets(Buffer, BUFFERLENGTH, pFile);
	pstr = Buffer;
	mf_GetAttributesFromFirstRow(pstr, mapAttrtibuteAndcolumns);

	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Majority protein IDs");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iMajorProteinIDcolumnNum = mapAttrtibuteAndcolumnsIter->second;
		bIfFindProteinID = true;
		count++;
	}
	else
	{
		bIfFindProteinID = false;
	}

	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Reverse");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iReversecolumnNum = mapAttrtibuteAndcolumnsIter->second;
		bIFFindReverse = true;
		count++;
	}
	else
	{
		bIFFindReverse = false;
	}

	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Potential contaminant");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iPoteinalContaminantcolumnNum = mapAttrtibuteAndcolumnsIter->second;
		bIfFindPoteinalContaminant = true;
		count++;
	}
	else
	{
		bIfFindPoteinalContaminant = false;
	}

	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Peptide IDs");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iPeptideIDscolumnNum = mapAttrtibuteAndcolumnsIter->second;
		bIfFindPeptideIDs = true;
		count++;
	}
	else
	{
		bIfFindPeptideIDs = false;
	}

	mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find("Peptide is razor");
	if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
	{
		iPeptideRazorcolumnNum = mapAttrtibuteAndcolumnsIter->second;
		bIfFindPeptideRazor = true;
		count++;
	}
	else
	{
		bIfFindPeptideRazor = false;
	}
	
	
	map<string, string>::iterator mapExperimentNameAndIBAQIntensityNameIter;
	mapExperimentNameAndIBAQIntensityNameIter = m_mapExperimentNameAndIBAQIntensityName.begin();
	for (; mapExperimentNameAndIBAQIntensityNameIter != m_mapExperimentNameAndIBAQIntensityName.end(); mapExperimentNameAndIBAQIntensityNameIter++)
	{
		mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find(strim(mapExperimentNameAndIBAQIntensityNameIter->second));
		if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
		{
			m_bIfIBAQIntensityExist = true;
			veciBAQcolumnNum.push_back(mapAttrtibuteAndcolumnsIter->second);
			vecbIfFindIBAQ.push_back(true);
			count++;
		}
		else
		{
			vecbIfFindIBAQ.push_back(false);
		}
	}

	map<string, string>::iterator mapExperimentNameAndLFQIntensityNameIter;
	mapExperimentNameAndLFQIntensityNameIter = m_mapExperimentNameAndLFQIntensityName.begin();
	for (; mapExperimentNameAndLFQIntensityNameIter != m_mapExperimentNameAndLFQIntensityName.end(); mapExperimentNameAndLFQIntensityNameIter++)
	{
		mapAttrtibuteAndcolumnsIter = mapAttrtibuteAndcolumns.find(strim(mapExperimentNameAndLFQIntensityNameIter->second));
		if (mapAttrtibuteAndcolumnsIter != mapAttrtibuteAndcolumns.end())
		{
			m_bIfLFQIntensityExist = true;
			veciLFQintensitycolumnNum.push_back(mapAttrtibuteAndcolumnsIter->second);
			vecbIfFindLFQIntensity.push_back(true);
			count++;
		}
		else
		{
			vecbIfFindLFQIntensity.push_back(false);
		}
	}

	if (!m_bIfIBAQIntensityExist)  iAttributeNumbers -= m_iNumberOfExperiments;
	if (!m_bIfLFQIntensityExist)   iAttributeNumbers -= (m_iNumberOfExperiments-1);

	if (count < iAttributeNumbers)
	{
		if (!bIfFindPeptideIDs)
		{
			cout << "Can not find the column \"Peptide IDs\" in the proteinGroups.txt" << endl;
			log << "Error:\tCan not find the column \"Peptide IDs\" in the proteinGroups.txt" << endl;
			exit(0);
		}
		if (!bIfFindPeptideRazor)
		{
			cout << "Can not find the column \"Peptide is razor\" in the proteinGroups.txt" << endl;
			log << "Error:\tCan not find the column \"Peptide is razor\" in the proteinGroups.txt" << endl;
			exit(0);
		}
		if (!bIfFindPoteinalContaminant)
		{
			cout << "Can not find the column \"Potential contaminant\" in the proteinGroups.txt" << endl;
			log << "Error:\tCan not find the column \"Potential contaminant\" in the proteinGroups.txt" << endl;
			exit(0);
		}
		if (!bIfFindProteinID)
		{
			cout << "Can not find the column \"Majority protein IDs\" in the proteinGroups.txt" << endl;
			log << "Error:\tCan not find the column \"Majority protein IDs\" in the proteinGroups.txt" << endl;
			exit(0);
		}
		if (!bIFFindReverse)
		{
			cout << "Can not find the column \"Reverse\" in the proteinGroups.txt" << endl;
			log << "Error:\tCan not find the column \"Reverse\" in the proteinGroups.txt" << endl;
			exit(0);
		}
		if (m_bIfIBAQIntensityExist)
		{
			for (int i = 0; i < vecbIfFindIBAQ.size(); i++)
			{
				if (!vecbIfFindIBAQ.at(i))
				{
					cout << "Can not find the column IBAQPeptideIntensityName in the proteinGroups.txt" << endl;
					log << "Error:\tCan not find the column IBAQPeptideIntensityName in the proteinGroups.txt" << endl;
					exit(0);

				}
			}
		}
		if (m_bIfLFQIntensityExist)
		{
			for (int i = 0; i < vecbIfFindLFQIntensity.size(); i++)
			{
				if (!vecbIfFindLFQIntensity.at(i))
				{
					cout << "Can not find the column LFQPeptideIntensityName in the proteinGroups.txt" << endl;
					log << "Can not find the column LFQPeptideIntensityName in the proteinGroups.txt" << endl;
				}
			}
		}
	}


	//
	icolumns = 0;
	count = 0;
	bool bIfDelete;
	CProtein cprotein;
	string strProteinIDTemp;
	string strSequenceTemp;
	vector<double> vecdIBAQ_IntensityTemp;
	vector<double> vecdLFQ_IntensityTemp;
	string strPeptidesIDTemp;
	string strPeptidesIDRazor;
	string strtemp;
	int itemp;
	double dIntensityTemp = 0.0;
	vector<double> vecdIntensityTemp;
	int begin1 = 0, end1 = 0;
	int begin2 = 0, end2 = 0;
	int NumberOfSharePeptides = 0;
	map<string, string>::iterator ProteinIDSequenceIter;

	int iNormalProteinsNumber;  
	bool b_ifExistReverseProtein = true;
	bool b_ifExistcontaminantProtein = true;
	//string strReversePrefix = "REV_";
	//string strcontaminantfix = "CON_";
	string strProteinNameTemp;
	char * pstr2; 
	bool bIfContainUPS2;
	int i = 0;
	bool bIfIBAQAllZeros;
	double dIBAQTemp;

	map<string, string>::iterator mapExperimentNamesIter;
	map<string, vector<double>>::iterator mapExperimentAndIntensitysIter;
	int iExperimentIndex = 0;


	fgets(Buffer, BUFFERLENGTH, pFile);
	pstr = Buffer;
	map<string, bool>::iterator mapSequenceAndIfUseIter;
	while (!feof(pFile))
	{
		bIfDelete = false;
		bIfIBAQAllZeros = true;
		count = 0;
		icolumns = 0;
		vecdIBAQ_IntensityTemp.clear();
		vecdLFQ_IntensityTemp.clear();
		while (count < iAttributeNumbers)
		{
			pstr1 = strstr(pstr, "\t");
			if (pstr1 != NULL)
			{
				*pstr1 = '\0';
			}
			else
			{
				log << "Error:\tThe format of proteinGroups.txt  is wrong.\n";
				exit(0);
			}

			if (icolumns == iMajorProteinIDcolumnNum)
			{				//delete the protein group if it contain reverse protein or only contain contaminate proteins。If left, take the first protein as its  proxy.
				cprotein.m_strProteinFullName = pstr;
				if (b_ifExistReverseProtein)
					if (strstr(pstr, m_strReversePrefix.c_str()) != NULL)
						bIfDelete = true;                      // delete the peptide, if its corresponding protein group contain reverse protein

				bIfContainUPS2 = false;
				//pstr2 = strstr(pstr, "ups");  // NO UPS2 first if take this code block into annotation 
				//if (pstr2 != NULL)
				//{
				//	bIfContainUPS2 = true;
				//}
				//pstr2 = strstr(pstr, "UPS");
				//if (pstr2 != NULL)
				//{
				//	bIfContainUPS2 = true;
				//}

				if (!b_ifExistcontaminantProtein)
				{
					if (bIfContainUPS2)
					{
						pstr2 = strstr(pstr, ";");
						if (pstr2 != NULL)
						{
							while (pstr2 != NULL)
							{
								*pstr2 = '\0';
								strProteinIDTemp = pstr;
								if ((strProteinIDTemp.find("ups") < strProteinIDTemp.size()) || (strProteinIDTemp.find("UPS") < strProteinIDTemp.size()))
								{
									break;
								}
								pstr = pstr2 + 1;
								pstr2 = strstr(pstr, ";");
							}
						}
						else
						{
							strProteinIDTemp = pstr;
						}
					}
					else
					{
						pstr2 = strstr(pstr, ";");
						if (pstr2 != NULL)
						{
							*pstr2 = '\0';
							strProteinIDTemp = pstr;
							pstr = pstr2 + 1;
						}
						else
						{
							strProteinIDTemp = pstr;
						}
					}

				}
				else
				{
					iNormalProteinsNumber = 0;
					if (bIfContainUPS2)
					{
						pstr2 = strstr(pstr, ";");
						if (pstr2 != NULL)
						{
							while (pstr2 != NULL)
							{

								*pstr2 = '\0';
								strProteinNameTemp = pstr;
								if ((strProteinNameTemp.find("ups") < strProteinNameTemp.size()) || (strProteinNameTemp.find("UPS") < strProteinNameTemp.size()))
								{
									strProteinIDTemp = pstr;
								}
								strProteinNameTemp = strProteinNameTemp.substr(0, m_strContaminantfix.length()); //前缀
								if (strProteinNameTemp != m_strContaminantfix)
									iNormalProteinsNumber++;

								pstr = pstr2 + 1;
								pstr2 = strstr(pstr, ";");
							}
							strProteinNameTemp = pstr;
							if ((strProteinNameTemp.find("ups") < strProteinNameTemp.size()) || (strProteinNameTemp.find("UPS") < strProteinNameTemp.size()))
							{
								strProteinIDTemp = pstr;
							}
							strProteinNameTemp = strProteinNameTemp.substr(0, m_strContaminantfix.length());
							if (strProteinNameTemp != m_strContaminantfix)
								iNormalProteinsNumber++;
						}
						else
						{   //only contain one protein
							strProteinNameTemp = pstr;
							strProteinIDTemp = pstr;
							strProteinNameTemp = strProteinNameTemp.substr(0, m_strContaminantfix.length()); //前缀
							if (strProteinNameTemp != m_strContaminantfix)
								iNormalProteinsNumber++;
						}

						if (iNormalProteinsNumber == 0)
							bIfDelete = true;
					}
					else
					{ // no UPS2, get the first ID 
						pstr2 = strstr(pstr, ";");
						if (pstr2 != NULL)
						{
							while (pstr2 != NULL)
							{
								*pstr2 = '\0';
								strProteinNameTemp = pstr;
								strProteinNameTemp = strProteinNameTemp.substr(0, m_strContaminantfix.length()); //前缀
								if (strProteinNameTemp != m_strContaminantfix)
								{
									strProteinIDTemp = pstr;
									iNormalProteinsNumber++;
									break;
								}
								pstr = pstr2 + 1;
								pstr2 = strstr(pstr, ";");
							}
							strProteinNameTemp = pstr;
							strProteinNameTemp = strProteinNameTemp.substr(0, m_strContaminantfix.length());
							if ((strProteinNameTemp != m_strContaminantfix) && (iNormalProteinsNumber == 0))
							{
								strProteinIDTemp = pstr;
								iNormalProteinsNumber++;
							}
						}
						else
						{   //only contain one protein
							strProteinNameTemp = pstr;
							strProteinIDTemp = pstr;
							strProteinNameTemp = strProteinNameTemp.substr(0, m_strContaminantfix.length()); //前缀
							if (strProteinNameTemp != m_strContaminantfix)
								iNormalProteinsNumber++;
						}
						if (iNormalProteinsNumber == 0)
							bIfDelete = true;
					}
				}
				count++;
			}
			
			for (i = 0; i < veciBAQcolumnNum.size(); i++)
			{
				if (icolumns == veciBAQcolumnNum.at(i))
				{
					dIBAQTemp = atof(pstr); 
					if (dIBAQTemp != 0.0)
					{
						bIfIBAQAllZeros = false;
					}
					vecdIBAQ_IntensityTemp.push_back(atof(pstr));
					count++;
				}
			}

			for (i = 0; i < veciLFQintensitycolumnNum.size(); i++)
			{
				if (icolumns == veciLFQintensitycolumnNum.at(i))
				{
					vecdLFQ_IntensityTemp.push_back(atof(pstr));
					count++;
				}
			}
			if (icolumns == iReversecolumnNum)
			{
				count++;
			}
			if (icolumns == iPoteinalContaminantcolumnNum)
			{
				count++;
			}
			if (icolumns == iPeptideIDscolumnNum)
			{
				strPeptidesIDTemp = pstr;
				count++;
			}
			if (icolumns == iPeptideRazorcolumnNum)
			{
				strPeptidesIDRazor = pstr;
				count++;
			}
			icolumns++;
			pstr = pstr1 + 1;
		}
		if (bIfIBAQAllZeros == true)
			bIfDelete = true;

		if (!bIfDelete)
		{		
			cprotein.m_strProteinID = strProteinIDTemp;
			int itemp1 = strProteinIDTemp.find("|");
			int itemp2 = strProteinIDTemp.find("|", itemp1 + 1);
			strProteinIDTemp = strProteinIDTemp.substr(itemp1 + 1, itemp2 - itemp1 - 1);
			if (strProteinIDTemp == "P07339ups" || strProteinIDTemp == "P08311ups") //change by gzhq 20151028 
			{  //deleteP07339 and P08311 in UPS2;add by gzhq 20150709
				fgets(Buffer, BUFFERLENGTH, pFile);
				pstr = Buffer;
				continue;
			}
			ProteinIDSequenceIter = m_mapProteinIDAndSequence.find(cprotein.m_strProteinID);
			if (ProteinIDSequenceIter != m_mapProteinIDAndSequence.end())
			{
				cprotein.m_strProteinSequence = ProteinIDSequenceIter->second;
			}
			else
			{
				cout << "Can not find protein: " << cprotein.m_strProteinID <<" in fasta file."<< endl;
				log << "Can not find protein: " << cprotein.m_strProteinID << " in fasta file." << endl;
				fgets(Buffer, BUFFERLENGTH, pFile);
				pstr = Buffer;
				continue;
			}
			for (i = 0; i < vecdIBAQ_IntensityTemp.size(); i++)
			{
				cprotein.m_vecdMaxquantIBAQ.push_back(vecdIBAQ_IntensityTemp.at(i));
			}
			for (i = 0; i < vecdLFQ_IntensityTemp.size(); i++)
			{
				cprotein.m_vecdMaxquantLFQ.push_back(vecdLFQ_IntensityTemp.at(i));
			}

			begin1 = 0;
			begin2 = 0;

			end1 = strPeptidesIDRazor.find(";", begin1 + 1);
			strtemp.clear();
			while (end1 != strPeptidesIDRazor.npos)
			{
				if (begin1 == 0)
					strtemp = strPeptidesIDRazor.substr(begin1, end1 - begin1);
				else
					strtemp = strPeptidesIDRazor.substr(begin1 + 1, end1 - begin1 - 1);
				end2 = strPeptidesIDTemp.find(";", begin2 + 1);
				if (begin2 == 0)
					itemp = atoi(strPeptidesIDTemp.substr(begin2, end2 - begin2).c_str());
				else
					itemp = atoi(strPeptidesIDTemp.substr(begin2 + 1, end2 - begin2 - end2 - 1).c_str());
				begin1 = end1;
				begin2 = end2;
				transform(strtemp.begin(), strtemp.end(), strtemp.begin(), ::tolower);//change to lowcase
				if (strtemp == "true")
				{
					cprotein.m_vPeptidesIDs.push_back(itemp);
					strtemp = mf_GetPeptidesSequenceFromID(log,itemp);
					if (strtemp == "null")
					{
						cprotein.m_vPeptidesIDs.pop_back();
						end1 = strPeptidesIDRazor.find(";", begin1 + 1);
						continue;
					}
					mapExperimentNamesIter = m_mapExperimentNameAndPeptideIntensityName.begin();
					iExperimentIndex = 0;
					for (; mapExperimentNamesIter != m_mapExperimentNameAndPeptideIntensityName.end(); mapExperimentNamesIter++)
					{
						dIntensityTemp = mf_GetPeptidesIntensityFromID(log,itemp, iExperimentIndex);
						iExperimentIndex++;
						mapExperimentAndIntensitysIter = cprotein.m_mapExperimentAndPeptidesIntensity.find(mapExperimentNamesIter->first);
						if (mapExperimentAndIntensitysIter == cprotein.m_mapExperimentAndPeptidesIntensity.end())
						{
							vecdIntensityTemp.clear();
							vecdIntensityTemp.push_back(dIntensityTemp);
							cprotein.m_mapExperimentAndPeptidesIntensity.insert(pair<string, vector<double>>(mapExperimentNamesIter->first, vecdIntensityTemp));
						}
						else
						{
							mapExperimentAndIntensitysIter->second.push_back(dIntensityTemp);
						}
					}//end for experiments

					cprotein.m_vPeptidesSequences.push_back(strtemp);
					m_NumberOfPeptides++;
	
				}

				end1 = strPeptidesIDRazor.find(";", begin1 + 1);

			}// end while

			if (begin2 == 0)
			{  //only one peptide
				strtemp = strPeptidesIDRazor.substr(begin1, strPeptidesIDRazor.size() - begin1);
				itemp = atoi(strPeptidesIDTemp.substr(begin2, strPeptidesIDTemp.size() - begin2).c_str());
			}
			else
			{
				strtemp = strPeptidesIDRazor.substr(begin1 + 1, strPeptidesIDRazor.size() - begin1 - 1);
				itemp = atoi(strPeptidesIDTemp.substr(begin2 + 1, strPeptidesIDTemp.size() - begin2 - 1).c_str());
			}

			transform(strtemp.begin(), strtemp.end(), strtemp.begin(), ::tolower);
			if (strtemp == "true")
			{
				strtemp = mf_GetPeptidesSequenceFromID(log,itemp);
				if (strtemp != "null")
				{
					cprotein.m_vPeptidesIDs.push_back(itemp);
					mapExperimentNamesIter = m_mapExperimentNameAndPeptideIntensityName.begin();
					iExperimentIndex = 0;
					for (; mapExperimentNamesIter != m_mapExperimentNameAndPeptideIntensityName.end(); mapExperimentNamesIter++)
					{
						dIntensityTemp = mf_GetPeptidesIntensityFromID(log,itemp, iExperimentIndex);
						iExperimentIndex++;
						mapExperimentAndIntensitysIter = cprotein.m_mapExperimentAndPeptidesIntensity.find(mapExperimentNamesIter->first);
						if (mapExperimentAndIntensitysIter == cprotein.m_mapExperimentAndPeptidesIntensity.end())
						{
							vecdIntensityTemp.clear();
							vecdIntensityTemp.push_back(dIntensityTemp);
							cprotein.m_mapExperimentAndPeptidesIntensity.insert(pair<string, vector<double>>(mapExperimentNamesIter->first, vecdIntensityTemp));
						}
						else
						{
							mapExperimentAndIntensitysIter->second.push_back(dIntensityTemp);
						}
					} //end for experiments
					cprotein.m_vPeptidesSequences.push_back(strtemp);
					m_NumberOfPeptides++;

				}

			}
			cprotein.m_iPeptidesNumber = cprotein.m_vPeptidesSequences.size();
			if (cprotein.m_vPeptidesSequences.size() > 0)
			{
				m_vecProteins.push_back(cprotein);
			}
			cprotein.Clear();

		} // end for if (strtemp == "true")

		fgets(Buffer, BUFFERLENGTH, pFile);
		pstr = Buffer;
	}
	fclose(pFile);
	cout << "\tLoaded " << m_vecProteins.size() << " Proteins  and " << m_NumberOfPeptides << " peptides! " << endl;
	log << "\tLoaded " << m_vecProteins.size() << " Proteins  and " << m_NumberOfPeptides << " peptides! " << endl;
	return 1;
}

void CLoadMaxquantIO::mf_saveProteins(std::ofstream& log, string proteinspath)
{
	ofstream ofile(proteinspath.c_str());
	if (!ofile)
	{
		cout << "Can not open " << proteinspath << endl;
		log << "Error:\tCan not open " << proteinspath << endl;
		exit(0);
	}
	else
	{
		cout << "\tSave Proteinpeptides in " << proteinspath << endl;
		log << "\tSave Proteinpeptides in " << proteinspath << endl;
	}
	ofile << "ProteinID\tMajority protein IDs\tPeptidesNumber\tvPeptidesSequence\t";
	map<string, string>::iterator mapExperimentIter;
	mapExperimentIter = m_mapExperimentNameAndPeptideIntensityName.begin();
	for (; mapExperimentIter!= m_mapExperimentNameAndPeptideIntensityName.end(); mapExperimentIter++)
	{
		ofile << "vPeptidesIntensitys " + mapExperimentIter->first+"\t";
	}
	ofile << "bIfIBAQExist\t";
	if (m_bIfIBAQIntensityExist)
	{
		mapExperimentIter = m_mapExperimentNameAndPeptideIntensityName.begin();
		for (; mapExperimentIter != m_mapExperimentNameAndPeptideIntensityName.end(); mapExperimentIter++)
		{
			ofile << "dMaxquantIBAQ " + mapExperimentIter->first + "\t";
		}

	}
	ofile << "bIfMaxquantLFQExist\t";
	if (m_bIfLFQIntensityExist)
	{
		mapExperimentIter = m_mapExperimentNameAndPeptideIntensityName.begin();
		for (mapExperimentIter++; mapExperimentIter != m_mapExperimentNameAndPeptideIntensityName.end(); mapExperimentIter++)
		{
			ofile << "dMaxquantLFQ " + mapExperimentIter->first + "\t";
		}
	}
	ofile << "\n";

	vector<CProtein>::iterator  ProteinIter;
	vector<string>::iterator vPeptideIter;
	vector<double>::iterator vPeptideIntensityIter;
	int i = 0;
	map<string, vector<double>>::iterator mapExperimentAndPeptideIntensitysIter;
	for (ProteinIter = m_vecProteins.begin(); ProteinIter != m_vecProteins.end(); ProteinIter++)
	{
		ofile << ProteinIter->m_strProteinID << "\t" << ProteinIter->m_strProteinFullName << "\t";
		ofile << ProteinIter->m_iPeptidesNumber << "\t";
		for (vPeptideIter = ProteinIter->m_vPeptidesSequences.begin(); vPeptideIter != ProteinIter->m_vPeptidesSequences.end(); vPeptideIter++)
			ofile << *vPeptideIter << ";";
		ofile << "\t";

		mapExperimentAndPeptideIntensitysIter = ProteinIter->m_mapExperimentAndPeptidesIntensity.begin();
		for (; mapExperimentAndPeptideIntensitysIter != ProteinIter->m_mapExperimentAndPeptidesIntensity.end(); mapExperimentAndPeptideIntensitysIter++)
		{
			for (i = 0; i < mapExperimentAndPeptideIntensitysIter->second.size(); i++)
			{
				ofile << mapExperimentAndPeptideIntensitysIter->second.at(i) << ";";
			}
			ofile << "\t";
		}
		
		if (m_bIfIBAQIntensityExist)
		{
			ofile << 1 << "\t";
			for (i = 0; i < ProteinIter->m_vecdMaxquantIBAQ.size(); i++)
			{
				ofile << ProteinIter->m_vecdMaxquantIBAQ.at(i) << "\t";
			}
		}
		else
		{
			ofile << 0 << "\t";
		}
		if (m_bIfLFQIntensityExist)
		{
			ofile << 1 << "\t";
			for (i = 0; i < ProteinIter->m_vecdMaxquantLFQ.size(); i++)
			{
				ofile << ProteinIter->m_vecdMaxquantLFQ.at(i) << "\t";
			}
		}
		else
		{
			ofile << 0 << "\t";
		}
		ofile << "\n";

	}
	ofile.close();
}

