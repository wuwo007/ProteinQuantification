#include"stdafx.h"
#include"QuantificationParam.h"
#include <sys/stat.h>
#include "BasicClass.h"
/*parameters
string m_strProteinsPath;
string m_strPeptidesInfoPath;
double m_dBestAlpha1;
double m_dBestAlpha2;
*/
int  CTrainParam::mf_setparameters(std::ofstream& log, string parafilepath)
{
	cout << "\tSetting Parameters\n";
	log << "\tSetting Parameters\n";

	ifstream fin(parafilepath.c_str());
	if (!fin)
	{
		cout << "Can not open " << parafilepath << endl; 
		log << "Error\t:Can not open " << parafilepath << endl;
		return 2;
		//exit(0);
	}
	map <string, string>  mapArttributesInFile;
	string strAttributeName, strAtributeValue;
	string strTemp;
	int iTemp1 = 0, iTemp2 = 0;

	while (getline(fin, strTemp))
	{
		iTemp2 = strTemp.find("=", 0);
		strAttributeName = strTemp.substr(0, iTemp2);
		iTemp1 = strTemp.find("\"");
		if (iTemp1 == strTemp.npos)
		{
			cout << "Parameter file format is not correct" << endl;
			log << "Error:\tParameter file format is not correct" << endl;
			exit(0);
		}
		iTemp2 = strTemp.find("\"", iTemp1 + 1);
		if (iTemp2 == strTemp.npos)
		{
			cout << "Parameter file format is not correct" << endl;
			log << "Error:\tParameter file format is not correct" << endl;
			exit(0);
		}
		strAtributeValue = strTemp.substr(iTemp1 + 1, iTemp2 - iTemp1 - 1);
		mapArttributesInFile.insert(pair<string, string>(strAttributeName, strAtributeValue));
	}

	map<string, string>::iterator mapAttributesInfileIter;
	
	mapAttributesInfileIter = mapArttributesInFile.find("IdentificationFileType");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		strIndentifySoftwareType = mapAttributesInfileIter->second;
	}
	else
	{
		cout << "Parameter file format is not correct" << endl;
		log << "Error:\tParameter file format is not correct" << endl;
		exit(0);
	}

	mapAttributesInfileIter = mapArttributesInFile.find("Fastapath");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_strProteinFastaPath = mapAttributesInfileIter->second;
	}
	else
	{
		cout << "Parameter file format is not correct" << endl;
		log << "Error:\tParameter file format is not correct" << endl;
		exit(0);
	}

	string strFastaTypetemp;
	mapAttributesInfileIter = mapArttributesInFile.find("fastaType");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		strFastaTypetemp = mapAttributesInfileIter->second; 
		m_eFastaType = strFastaTypetemp;
	}
	else
	{
		cout << "Parameter file format is not correct" << endl;
		log << "Error:\tParameter file format is not correct" << endl;
		exit(0);
	}

	mapAttributesInfileIter = mapArttributesInFile.find("identification  result File path");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_strMaxquantResultPath = mapAttributesInfileIter->second;
		m_strExperimentDesignPath = m_strMaxquantResultPath + "\\experimentalDesignTemplate.txt";
	}
	else
	{
		cout << "Parameter file format is not correct" << endl;
		log << "Error:\tParameter file format is not correct" << endl;
		exit(0);
	}

	mapAttributesInfileIter = mapArttributesInFile.find("ResultPath");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_strQuantiResultPath = mapAttributesInfileIter->second;
	}
	else
	{
		cout << "Parameter file format is not correct" << endl;
		log << "Error:\tParameter file format is not correct" << endl;
		exit(0);
	}
	
	mapAttributesInfileIter = mapArttributesInFile.find("RegressionMethod");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		 m_strRegressionMethod= mapAttributesInfileIter->second;
	}
	else
	{
		cout << "Parameter file format is not correct" << endl;
		log << "Error:\tParameter file format is not correct" << endl;
		exit(0);
	}

	mapAttributesInfileIter = mapArttributesInFile.find("alpha1");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_dAlpha1 =  atof(mapAttributesInfileIter->second.c_str());
	}
	else
	{
		cout << "Parameter file format is not correct" << endl;
		log << "Error:\tParameter file format is not correct" << endl;
		exit(0);
	}

	mapAttributesInfileIter = mapArttributesInFile.find("Alpha2");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_dAlpha2 = atof(mapAttributesInfileIter->second.c_str());
	}
	else
	{
		cout << "Parameter file format is not correct" << endl;
		log << "Error:\tParameter file format is not correct" << endl;
		exit(0);
	}
	
	mapAttributesInfileIter = mapArttributesInFile.find("alpha");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_dAlpha = atof(mapAttributesInfileIter->second.c_str());
	}
	else
	{
		cout << "Parameter file format is not correct" << endl;
		log << "Error:\tParameter file format is not correct" << endl;
		exit(0);
	}

	mapAttributesInfileIter = mapArttributesInFile.find("beta");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_dBeta = atof(mapAttributesInfileIter->second.c_str());
	}
	else
	{
		cout << "Parameter file format is not correct" << endl;
		log << "Error:\tParameter file format is not correct" << endl;
		exit(0);
	}
		
	mapAttributesInfileIter = mapArttributesInFile.find("k");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		 m_ik= atoi(mapAttributesInfileIter->second.c_str());
	}
	else
	{
		cout << "Parameter file format is not correct" << endl;
		log << "Error:\tParameter file format is not correct" << endl;
		exit(0);
	}
	
	mapAttributesInfileIter = mapArttributesInFile.find("Number of trees");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_iNumberOfTrees = atoi(mapAttributesInfileIter->second.c_str());
	}
	else
	{
		cout << "Parameter file format is not correct" << endl;
		log << "Error:\tParameter file format is not correct" << endl;
		exit(0);
	}

	mapAttributesInfileIter = mapArttributesInFile.find("bIfOptimizeParameters");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_bIfOptimizeParameters = fStringToBool(mapAttributesInfileIter->second);
	}
	else
	{
		cout << "Parameter file format is not correct" << endl;
		log << "Error:\tParameter file format is not correct" << endl;
		exit(0);
	}

	mapAttributesInfileIter = mapArttributesInFile.find("bIfCotainStardProtein");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_bIfCotainStardProtein = fStringToBool(mapAttributesInfileIter->second);
	}
	else
	{
		cout << "Parameter file format is not correct" << endl;
		log << "Error:\tParameter file format is not correct" << endl;
		exit(0);
	}

	mapAttributesInfileIter = mapArttributesInFile.find("strIdentifierOfStandPro");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_strIdentifierOfStandPro= mapAttributesInfileIter->second;
	}
	else
	{
		cout << "Parameter file format is not correct" << endl;
		log << "Error:\tParameter file format is not correct" << endl;
		exit(0);
	}

	mapAttributesInfileIter = mapArttributesInFile.find("strAbundanceOfStandProteinsPath");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_strAbundanceOfStandProteinsPath = mapAttributesInfileIter->second;
	}
	else
	{
		cout << "Parameter file format is not correct" << endl;
		log << "Error:\tParameter file format is not correct" << endl;
		exit(0);
	}

	m_strProteinsPath = m_strQuantiResultPath + "\\Proteins.txt";
	m_strRegressionResult = m_strQuantiResultPath + "\\regressionResult";
	return 0;

}
void  CPredictParam::mf_setparameters(std::ofstream& log, string parafilepath)
{
	cout << "\tSetting Parameters\n";
	log << "\tSetting Parameters\n";

	ifstream fin(parafilepath.c_str());
	if (!fin)
	{
		cout << "Can not open " << parafilepath << endl;
		log << "Error:\tCan not open " << parafilepath << endl;
		exit(0);
	}

	map <string, string>  mapArttributesInFile;
	string strAttributeName, strAtributeValue;
	string strTemp;
	int iTemp1 = 0, iTemp2 = 0;

	while (getline(fin, strTemp))
	{
		iTemp2 = strTemp.find("=", 0);
		strAttributeName = strTemp.substr(0, iTemp2);
		iTemp1 = strTemp.find("\"");
		if (iTemp1 == strTemp.npos)
		{
			cout << "Parameter file format is not correct" << endl;
			log << "Error:\tParameter file format is not correct" << endl;
			exit(0);
		}
		iTemp2 = strTemp.find("\"", iTemp1 + 1);
		if (iTemp2 == strTemp.npos)
		{
			cout << "Parameter file format is not correct" << endl;
			log << "Error:\tParameter file format is not correct" << endl;
			exit(0);
		}
		strAtributeValue = strTemp.substr(iTemp1 + 1, iTemp2 - iTemp1 - 1);
		mapArttributesInFile.insert(pair<string, string>(strAttributeName, strAtributeValue));
	}

	map<string, string>::iterator mapAttributesInfileIter;

	mapAttributesInfileIter = mapArttributesInFile.find("IdentificationFileType");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		strIndentifySoftwareType = mapAttributesInfileIter->second;
	}
	else
	{
		cout << "Parameter file format is not correct" << endl;
		log << "Error:\tParameter file format is not correct" << endl;
		exit(0);
	}

	mapAttributesInfileIter = mapArttributesInFile.find("Fastapath");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_strProteinFastaPath = mapAttributesInfileIter->second;
	}
	else
	{
		cout << "Parameter file format is not correct" << endl;
		log << "Error:\tParameter file format is not correct" << endl;
		exit(0);
	}

	string strFastaTypetemp;
	mapAttributesInfileIter = mapArttributesInFile.find("fastaType");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		strFastaTypetemp = mapAttributesInfileIter->second;
		m_eFastaType = strFastaTypetemp;
	}
	else
	{
		cout << "Parameter file format is not correct" << endl;
		log << "Error:\tParameter file format is not correct" << endl;
		exit(0);
	}

	mapAttributesInfileIter = mapArttributesInFile.find("identification  result File path");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_strMaxquantResultPath = mapAttributesInfileIter->second;
		m_strExperimentDesignPath = m_strMaxquantResultPath + "\\experimentalDesignTemplate.txt";
	}
	else
	{
		cout << "Parameter file format is not correct" << endl;
		log << "Error:\tParameter file format is not correct" << endl;
		exit(0);
	}

	mapAttributesInfileIter = mapArttributesInFile.find("ResultPath");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_strQuantiResultPath = mapAttributesInfileIter->second;
	}
	else
	{
		cout << "Parameter file format is not correct" << endl;
		log << "Error:\tParameter file format is not correct" << endl;
		exit(0);
	}

	mapAttributesInfileIter = mapArttributesInFile.find("MaxMissedCleavage");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_imaxMissedClevage = atoi(mapAttributesInfileIter->second.c_str());
	}
	else
	{
		cout << "Parameter file format is not correct" << endl;
		log << "Error:\tParameter file format is not correct" << endl;
		exit(0);
	}

	mapAttributesInfileIter = mapArttributesInFile.find("PepShortestLen");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_iMinPepLength = atoi(mapAttributesInfileIter->second.c_str());
	}
	else
	{
		cout << "Parameter file format is not correct" << endl;
		log << "Error:\tParameter file format is not correct" << endl;
		exit(0);
	}

	mapAttributesInfileIter = mapArttributesInFile.find("PepLongestLen");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_iMaxPepLength = atoi(mapAttributesInfileIter->second.c_str());
	}
	else
	{
		cout << "Parameter file format is not correct" << endl;
		log << "Error:\tParameter file format is not correct" << endl;
		exit(0);
	}

	mapAttributesInfileIter = mapArttributesInFile.find("bIfCotainStardProtein");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_bIfCotainStardProtein = fStringToBool(mapAttributesInfileIter->second);
	}
	else
	{
		cout << "Parameter file format is not correct" << endl;
		log << "Error:\tParameter file format is not correct" << endl;
		exit(0);
	}

	mapAttributesInfileIter = mapArttributesInFile.find("strIdentifierOfStandPro");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_strIdentifierOfStandPro = mapAttributesInfileIter->second;
	}
	else
	{
		cout << "Parameter file format is not correct" << endl;
		log << "Error:\tParameter file format is not correct" << endl;
		exit(0);
	}

	mapAttributesInfileIter = mapArttributesInFile.find("strAbundanceOfStandProteinsPath");
	if (mapAttributesInfileIter != mapArttributesInFile.end())
	{
		m_strAbundanceOfStandProteinsPath = mapAttributesInfileIter->second;
	}
	else
	{
		cout << "Parameter file format is not correct" << endl;
		log << "Error:\tParameter file format is not correct" << endl;
		exit(0);
	}
	m_strCutPosition = "KR";
	m_blr = false;
	m_strProteinsPath = m_strQuantiResultPath + "\\Proteins.txt";
	m_strRegressionResult = m_strQuantiResultPath + "\\regressionResult";
	m_strProteinSPEAQInPeptidesPath = m_strQuantiResultPath + "\\ProteinSPEAQs";
}

