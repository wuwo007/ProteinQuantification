#include"stdafx.h"
#include"loadParam.h"
#include"BasicClass.h"

CLoadParam::CLoadParam()
{
}
/*
load the parameters
*/
void  CLoadParam::mf_setparameters(std::ofstream& log, string parafilepath)
{
	cout << "\tSetting Parameters from " << parafilepath << "\n";
	log << "\tSetting Parameters from " << parafilepath << "\n";
	ifstream fin(parafilepath.c_str());
	if (!fin)
	{
		cout << "can not open the file " << parafilepath << endl;
		log << "Error:\tcan not open the file " << parafilepath << endl;
		exit(0);
	}
	string strTemp1, strTemp2;
	int iTemp1 = 0, iTemp2 = 0;
	string strFastaTypetemp;

	while (getline(fin, strTemp1))
	{
		iTemp1 = strTemp1.find("=");
		if (iTemp1 != strTemp1.npos)
		{
			strTemp2 = strTemp1.substr(0, iTemp1);
			transform(strTemp2.begin(), strTemp2.end(), strTemp2.begin(), ::tolower);//change to lowcase
			if (strTemp2 == "identificationfiletype")
			{
				iTemp1 = strTemp1.find("\"", iTemp1);
				if (iTemp1 == strTemp1.npos)
				{
					cout << "The format of parameter identificationfiletype is not correct!" << endl;
					log << "Error:\tThe format of parameter identificationfiletype is not correct!" << endl;
					exit(0);
				}
				iTemp2 = strTemp1.find("\"", iTemp1 + 1);
				if (iTemp2 == strTemp1.npos)
				{
					cout << "The format of parameter identificationfiletype is not correct!" << endl;
					log << "Error:\tThe format of parameter identificationfiletype is not correct!" << endl;
					exit(0);
				}
				strLoadFileType = strTemp1.substr(iTemp1 + 1, iTemp2 - iTemp1 - 1);
				if (strLoadFileType == "maxquant")
				{
					m_eDataType = MaxquantTpye;
				}
				else if (strLoadFileType == "pfind")
				{
					m_eDataType = pfindType;
				}
				else
				{
					cout << "The format of parameter identificationfiletype is not correct!" << endl;
					log << "Error:\tThe format of parameter identificationfiletype is not correct!" << endl;
					exit(0);
				}
			}
		}

		iTemp1 = strTemp1.find("=");
		if (iTemp1 != strTemp1.npos)
		{
			strTemp2 = strTemp1.substr(0, iTemp1);
			transform(strTemp2.begin(), strTemp2.end(), strTemp2.begin(), ::tolower);//change to lowcase
			if (strTemp2 == "identification  result file path")
			{
				iTemp1 = strTemp1.find("\"", iTemp1);
				if (iTemp1 == strTemp1.npos)
				{
					cout << "The format of parameter identificationfiletype is not correct!" << endl;
					log << "Error:\tThe format of parameter identificationfiletype is not correct!" << endl;
					exit(0);
				}
				iTemp2 = strTemp1.find("\"", iTemp1 + 1);
				if (iTemp2 == strTemp1.npos)
				{
					cout << "The format of parameter identificationfiletype is not correct!" << endl;
					log << "Error:\tThe format of parameter identificationfiletype is not correct!" << endl;
					exit(0);
				}
				m_strMaxquantResultPath = strTemp1.substr(iTemp1 + 1, iTemp2 - iTemp1 - 1);
			}
		}
		iTemp1 = strTemp1.find("=");
		if (iTemp1 != strTemp1.npos)
		{
			strTemp2 = strTemp1.substr(0, iTemp1);
			transform(strTemp2.begin(), strTemp2.end(), strTemp2.begin(), ::tolower);//change to lowcase
			if (strTemp2 == "fastapath")
			{
				iTemp1 = strTemp1.find("\"", iTemp1);
				if (iTemp1 == strTemp1.npos)
				{
					cout << "The format of parameter Fastapath is not correct!" << endl;
					log << "Error:\tThe format of parameter Fastapath is not correct!" << endl;
					exit(0);
				}
				iTemp2 = strTemp1.find("\"", iTemp1 + 1);
				if (iTemp2 == strTemp1.npos)
				{
					cout << "The format of parameter Fastapath is not correct!" << endl;
					log << "Error:\tThe format of parameter Fastapath is not correct!" << endl;
					exit(0);
				}
				m_strProteinFasterFilePath = strTemp1.substr(iTemp1 + 1, iTemp2 - iTemp1 - 1);
			}
		}

		iTemp1 = strTemp1.find("=");
		if (iTemp1 != strTemp1.npos)
		{
			strTemp2 = strTemp1.substr(0, iTemp1);
			transform(strTemp2.begin(), strTemp2.end(), strTemp2.begin(), ::tolower);//change to lowcase
			if (strTemp2 == "fastatype")
			{
				iTemp1 = strTemp1.find("\"", iTemp1);
				if (iTemp1 == strTemp1.npos)
				{
					cout << "The format of parameter fastatype is not correct!" << endl;
					log << "Error:\tThe format of parameter fastatype is not correct!" << endl;
					exit(0);
				}
				iTemp2 = strTemp1.find("\"", iTemp1 + 1);
				if (iTemp2 == strTemp1.npos)
				{
					cout << "The format of parameter fastatype is not correct!" << endl;
					log << "Error:\tThe format of parameter fastatype is not correct!" << endl;
					exit(0);
				}
				strFastaTypetemp = strTemp1.substr(iTemp1 + 1, iTemp2 - iTemp1 - 1);
				m_eFastaType = strFastaTypetemp;

			}
		}

		iTemp1 = strTemp1.find("=");
		if (iTemp1 != strTemp1.npos)
		{
			strTemp2 = strTemp1.substr(0, iTemp1);
			transform(strTemp2.begin(), strTemp2.end(), strTemp2.begin(), ::tolower);//change to lowcase
			if (strTemp2 == "ifexistreverseproteins")
			{
				iTemp1 = strTemp1.find("\"", iTemp1);
				if (iTemp1 == strTemp1.npos)
				{
					cout << "The format of parameter IfExistReverseProteins is not correct!" << endl;
					log << "Error:\tThe format of parameter IfExistReverseProteins is not correct!" << endl;
					exit(0);
				}
				iTemp2 = strTemp1.find("\"", iTemp1 + 1);
				if (iTemp2 == strTemp1.npos)
				{
					cout << "The format of parameter IfExistReverseProteins is not correct!" << endl;
					log << "Error:\tThe format of parameter IfExistReverseProteins is not correct!" << endl;
					exit(0);
				}
				m_bIfExistReverseProtein = fStringToBool(strTemp1.substr(iTemp1 + 1, iTemp2 - iTemp1 - 1));
			}
		}

		iTemp1 = strTemp1.find("=");
		if (iTemp1 != strTemp1.npos)
		{
			strTemp2 = strTemp1.substr(0, iTemp1);
			transform(strTemp2.begin(), strTemp2.end(), strTemp2.begin(), ::tolower);//change to lowcase
			if (strTemp2 == "prefixofreverseprotein")
			{
				iTemp1 = strTemp1.find("\"", iTemp1);
				if (iTemp1 == strTemp1.npos)
				{
					cout << "The format of parameter PrefixOfReverseProtein is not correct!" << endl;
					log << "Error:\tThe format of parameter PrefixOfReverseProtein is not correct!" << endl;
					exit(0);
				}
				iTemp2 = strTemp1.find("\"", iTemp1 + 1);
				if (iTemp2 == strTemp1.npos)
				{
					cout << "The format of parameter PrefixOfReverseProtein is not correct!" << endl;
					log << "Error:\tThe format of parameter PrefixOfReverseProtein is not correct!" << endl;
					exit(0);
				}
				m_strReverseProteinIDPrefix = strTemp1.substr(iTemp1 + 1, iTemp2 - iTemp1 - 1);
			}
		}

		iTemp1 = strTemp1.find("=");
		if (iTemp1 != strTemp1.npos)
		{
			strTemp2 = strTemp1.substr(0, iTemp1);
			transform(strTemp2.begin(), strTemp2.end(), strTemp2.begin(), ::tolower);//change to lowcase
			if (strTemp2 == "ifexistcontaminantproteins")
			{
				iTemp1 = strTemp1.find("\"", iTemp1);
				if (iTemp1 == strTemp1.npos)
				{
					cout << "The format of parameter IfExistContaminantProteins is not correct!" << endl;
					log << "Error:\tThe format of parameter IfExistContaminantProteins is not correct!" << endl;
					exit(0);
				}
				iTemp2 = strTemp1.find("\"", iTemp1 + 1);
				if (iTemp2 == strTemp1.npos)
				{
					cout << "The format of parameter IfExistContaminantProteins is not correct!" << endl;
					log << "Error:\tThe format of parameter IfExistContaminantProteins is not correct!" << endl;
					exit(0);
				}
				m_bIfExistContanProtein = fStringToBool(strTemp1.substr(iTemp1 + 1, iTemp2 - iTemp1 - 1));
			}
		}

		iTemp1 = strTemp1.find("=");
		if (iTemp1 != strTemp1.npos)
		{
			strTemp2 = strTemp1.substr(0, iTemp1);
			transform(strTemp2.begin(), strTemp2.end(), strTemp2.begin(), ::tolower);//change to lowcase
			if (strTemp2 == "prefixofcontaminantprotein")
			{
				iTemp1 = strTemp1.find("\"", iTemp1);
				if (iTemp1 == strTemp1.npos)
				{
					cout << "The format of parameter PrefixOfContaminantProtein is not correct!" << endl;
					log << "Error:\tThe format of parameter PrefixOfContaminantProtein is not correct!" << endl;
					exit(0);
				}
				iTemp2 = strTemp1.find("\"", iTemp1 + 1);
				if (iTemp2 == strTemp1.npos)
				{
					cout << "The format of parameter PrefixOfContaminantProtein is not correct!" << endl;
					log << "Error:\tThe format of parameter PrefixOfContaminantProtein is not correct!" << endl;
					exit(0);
				}
				m_strContaminantProteinIDPrefix = strTemp1.substr(iTemp1 + 1, iTemp2 - iTemp1 - 1);
			}
		}

	
		iTemp1 = strTemp1.find("=");
		if (iTemp1 != strTemp1.npos)
		{
			strTemp2 = strTemp1.substr(0, iTemp1);
			transform(strTemp2.begin(), strTemp2.end(), strTemp2.begin(), ::tolower);//change to lowcase
			if (strTemp2 == "resultpath")
			{
				iTemp1 = strTemp1.find("\"", iTemp1);
				if (iTemp1 == strTemp1.npos)
				{
					cout << "The format of parameter ResultPath is not correct!" << endl;
					log << "Error:\tThe format of parameter ResultPath is not correct!" << endl;
					exit(0);
				}
				iTemp2 = strTemp1.find("\"", iTemp1 + 1);
				if (iTemp2 == strTemp1.npos)
				{
					cout << "The format of parameter ResultPath is not correct!" << endl;
					log << "Error:\tThe format of parameter ResultPath is not correct!" << endl;
					exit(0);
				}
				m_strQuantiResultPath = strTemp1.substr(iTemp1 + 1, iTemp2 - iTemp1 - 1);

			}
		}
		m_strPeptideFilePath = m_strMaxquantResultPath + "\\peptides.txt";
		m_strProteinFilePath = m_strMaxquantResultPath + "\\proteinGroups.txt";
		m_strExprimentDesignPath = m_strMaxquantResultPath + "\\experimentalDesignTemplate.txt";
		m_strProteinsPath = m_strQuantiResultPath + "\\Proteins.txt";
	}
}

void CLoadParam::mf_showparam()
{
	cout << "Peptides come from " << m_strPeptideFilePath << endl;
	cout << "Proteins come from " << m_strProteinFilePath << endl;
	cout << "And the corresponding fasta file is " << m_strProteinFasterFilePath << endl;
	cout << "We need Expriment Design File " << m_strExprimentDesignPath << endl;
}