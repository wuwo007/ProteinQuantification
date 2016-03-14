//version 2015-12-04

#include"stdafx.h"
#include"ProteinWorker.h"
#include"ProteinDigestion.h"
#include"Regression.h"

// load information of proteins from proteins.txt
bool CProteinWorker::mf_LoadProteins(std::ofstream& log, CTrainParam trainParam, int iExperimentIndex)
{
	cout << "\tLoading Protein!\n";
	log << "\tLoading Protein!\n";
	//open the file
	FILE * pFile;
	pFile = fopen(trainParam.m_strProteinsPath.c_str(), "r");
	if (pFile == NULL)
	{
		cout << "Can not open " << trainParam.m_strProteinsPath << endl;
		log << "Error:\tCan not open "  << trainParam.m_strProteinsPath << endl;
		exit(0);
	}

	//temporary variable
	char Buffer[BUFFERLENGTH];
	char *pstr;
	char *pstr1;

	map<string, string> mapAttributeAndColumns;
	map<string, string>::iterator mapAttributeAndColumnsIter;

	CProtein ProteinTemp;
	CPeptideAttribute peptideAttributeTemp;
	PepAttributeWorker PeptideAttributeWorkerTemp;
	string strWithFlankingRegionTemp;
	string strtemp, substrtemp;
	int begin, end;
	map<string, string>::iterator mapProteinIDAndSequenceIter;
	map<int, string>::iterator PeptidesIDAndSequenceIter;
	int i = 0;
	map<string, double> mapPeptideAndIntensitys;
	map<string, double>::iterator mapPeptideAndIntensityIter;
	int iPeptidesNumberOfOneProtein = 0;

	fgets(Buffer, BUFFERLENGTH, pFile);    //jump the first row
	fgets(Buffer, BUFFERLENGTH, pFile);
	pstr = Buffer;

	//load the information of amino acid 
	m_vecAAindexAttributeHeaders = PeptideAttributeWorkerTemp.mf_LoadAAindex(log);
	m_vecProteins.clear();
 	m_iPeptidesNumberOfProteins = 0;

	while (!feof(pFile))
	{
		ProteinTemp.Clear();
		mapPeptideAndIntensitys.clear();
		pstr1 = strstr(pstr, "\t");
		if (pstr1 != NULL)
		{
			*pstr1 = '\0';
		}
		else
		{
			log << "Error:\tThe format of proteins.txt is wrong.\n";
			exit(0);
		}
		ProteinTemp.m_strProteinID = pstr;
		pstr = pstr1 + 1;

		mapProteinIDAndSequenceIter = m_mapProteinsIdAndSequence.find(ProteinTemp.m_strProteinID);
		if (mapProteinIDAndSequenceIter != m_mapProteinsIdAndSequence.end())
		{
			ProteinTemp.m_strProteinSequence = mapProteinIDAndSequenceIter->second;
		}
		else
		{
			cout << "Can not find Protein " << ProteinTemp.m_strProteinID << endl;
			log << "Error:\tCan not find Protein " << ProteinTemp.m_strProteinID << endl;
			exit(0);
		}

		pstr1 = strstr(pstr, "\t");
		if (pstr1 != NULL)
		{
			*pstr1 = '\0';
		}
		else
		{
			log << "Error:\tThe format of proteins.txt is wrong.\n";
			exit(0);
		}
		ProteinTemp.m_strProteinFullName = pstr;
		pstr = pstr1 + 1;

		pstr1 = strstr(pstr, "\t");
		pstr = pstr1 + 1;

		pstr1 = strstr(pstr, "\t");
		if (pstr1 != NULL)
		{
			*pstr1 = '\0';
		}
		else
		{
			log << "Error:\tThe format of proteins.txt is wrong.\n";
			exit(0);
		}
		strtemp = pstr;
		begin = 0;
		end = strtemp.find(";", begin);


		while (end != strtemp.npos)
		{   //many peptides
			substrtemp = strtemp.substr(begin, end - begin);
			mapPeptideAndIntensitys.insert(pair<string, double>(substrtemp, 0.0));
			begin = end + 1;
			end = strtemp.find(";", begin);
		}


		for (i = 0; i < iExperimentIndex + 1; i++)
		{
			pstr = pstr1 + 1;
			pstr1 = strstr(pstr, "\t");
		}
		if (pstr1 != NULL)
		{
			*pstr1 = '\0';
		}
		else
		{
			log << "Error:\tThe format of proteins.txt is wrong.\n";
			exit(0);
		}
		strtemp = pstr;
		begin = 0;
		end = strtemp.find(";", begin);

		mapPeptideAndIntensityIter = mapPeptideAndIntensitys.begin();
		while (end != strtemp.npos)
		{
			substrtemp = strtemp.substr(begin, end - begin);
			if (atof(substrtemp.c_str())==0.0)
			{
				mapPeptideAndIntensityIter=mapPeptideAndIntensitys.erase(mapPeptideAndIntensityIter);
			}
			else
			{
				mapPeptideAndIntensityIter->second = atof(substrtemp.c_str());
				mapPeptideAndIntensityIter++;
			}
			begin = end + 1;
			end = strtemp.find(";", begin);
		}
		 iPeptidesNumberOfOneProtein= mapPeptideAndIntensitys.size();
		 m_iPeptidesNumberOfProteins += iPeptidesNumberOfOneProtein;
		mapPeptideAndIntensityIter = mapPeptideAndIntensitys.begin();
		for (; mapPeptideAndIntensityIter != mapPeptideAndIntensitys.end(); mapPeptideAndIntensityIter++)
		{
			strWithFlankingRegionTemp = ProteinTemp.mf_GetPeptidesAdjacentSequence(log,mapPeptideAndIntensityIter->first);
			if (strWithFlankingRegionTemp == "NULL")
			{
				cout << "can not find peptide " << mapPeptideAndIntensityIter->first << "'s adjacent Sequence!\n";
				log << "can not find peptide " << mapPeptideAndIntensityIter->first << "'s adjacent Sequence!\n";
			}
			else
			{
				ProteinTemp.m_vPeptidesSequences.push_back(mapPeptideAndIntensityIter->first);
				ProteinTemp.m_vPeptidesSequencesWithFlankingRegion.push_back(strWithFlankingRegionTemp);                                         //根据肽段ID取得相应肽段序列和左右15位氨基酸
				peptideAttributeTemp.Clear();
				peptideAttributeTemp.m_vecAttributes = PeptideAttributeWorkerTemp.mf_GetAttributeFromSequence(mapPeptideAndIntensityIter->first, strWithFlankingRegionTemp);
				ProteinTemp.m_vPeptidesAttributes.push_back(peptideAttributeTemp);
				ProteinTemp.m_vdPeptidesMW.push_back(peptideAttributeTemp.m_vecAttributes.at(2));  //第三个属性是MW

				ProteinTemp.m_vPeptidesNativeIntensity.push_back(mapPeptideAndIntensityIter->second);
				ProteinTemp.m_vPeptidesAdjustIntensity.push_back(mapPeptideAndIntensityIter->second); //assign the same value with m_vPeptidesNativeIntensity for iterration
				ProteinTemp.m_vAdjustPeptidesIntensityFactor.push_back(0.0);
			}			

		}

		for (i = 0; i < m_iNumberOfExperiments - iExperimentIndex - 1; i++)
		{
			pstr = pstr1 + 1;
			pstr1 = strstr(pstr, "\t");
		}

		pstr = pstr1 + 1;
		pstr1 = strstr(pstr, "\t");
		if (pstr1 != NULL)
		{
			*pstr1 = '\0';
		}
		else
		{
			log << "Error:\tThe format of proteins.txt is wrong.\n";
			exit(0);
		}
		m_bIfIBAQIntensityExist = atoi(pstr);
		pstr = pstr1 + 1;

		for (i = 0; i < iExperimentIndex + 1; i++)
		{
			pstr = pstr1 + 1;
			pstr1 = strstr(pstr, "\t");
		}
		pstr1 = strstr(pstr, "\t");
		if (pstr1 != NULL)
		{
			*pstr1 = '\0';
		}
		else
		{
			log << "Error:\tThe format of proteins.txt is wrong.\n";
			exit(0);
		}
		ProteinTemp.m_dMaxquantIBAQ = atof(pstr);


		for (i = 0; i < m_iNumberOfExperiments - iExperimentIndex - 1; i++)
		{
			pstr = pstr1 + 1;
			pstr1 = strstr(pstr, "\t");
		}

		pstr = pstr1 + 1;
		pstr1 = strstr(pstr, "\t");
		if (pstr1 != NULL)
		{
			*pstr1 = '\0';
		}
		else
		{
			log << "Error:\tThe format of proteins.txt is wrong.\n";
			exit(0);
		}
		m_bIFLFQIntensityExist = atof(pstr);
		pstr = pstr1 + 1;

		for (i = 0; i < iExperimentIndex; i++)
		{
			pstr = pstr1 + 1;
			pstr1 = strstr(pstr, "\t");
		}
		if (m_bIFLFQIntensityExist&&iExperimentIndex>0)
		{
			pstr1 = strstr(pstr, "\t");
			if (pstr1 != NULL)
			{
				*pstr1 = '\0';
			}
			else
			{
				log << "Error:\tThe format of proteins.txt is wrong.\n";
				exit(0);
			}
			ProteinTemp.m_dLFQ = atof(pstr);
			pstr = pstr1 + 1;
		}
		else
		{
			ProteinTemp.m_dLFQ = 0.0;
		}

		ProteinTemp.m_iPeptidesNumber = ProteinTemp.m_vPeptidesSequences.size();

		//divide the median of peptides' intensities of one protein

		//ProteinTemp.m_dPeptidesIntensityMedian = GetMedianIntensity(ProteinTemp.m_vPeptidesNativeIntensity);
		//if (ProteinTemp.m_dPeptidesIntensityMedian != 0.0)
		//{
		//	for (size_t i = 0; i < ProteinTemp.m_vAdjustPeptidesIntensityFactor.size(); i++)
		//	ProteinTemp.m_vAdjustPeptidesIntensityFactor.at(i) = ProteinTemp.m_vPeptidesNativeIntensity.at(i) / ProteinTemp.m_dPeptidesIntensityMedian;
		//}	
		//divide the maximum of peptides' intensities of one protein	
		ProteinTemp.m_dPeptidesIntensityMax = GetMaxIntensity(ProteinTemp.m_vPeptidesNativeIntensity);
		if (ProteinTemp.m_dPeptidesIntensityMax != 0.0)
		{
			for (size_t i = 0; i < ProteinTemp.m_vAdjustPeptidesIntensityFactor.size(); i++)
				ProteinTemp.m_vAdjustPeptidesIntensityFactor.at(i) = ProteinTemp.m_vPeptidesNativeIntensity.at(i) / ProteinTemp.m_dPeptidesIntensityMax;
		}

		if (iPeptidesNumberOfOneProtein > 0)
		{
			m_vecProteins.push_back(ProteinTemp);
		}

		fgets(Buffer, BUFFERLENGTH, pFile);
		pstr = Buffer;
	}
	cout << "Loaded " << m_vecProteins.size() << " Proteins with " << m_iPeptidesNumberOfProteins << " peptides in experiment "<<m_vecStrExperiments.at(iExperimentIndex) << endl;
	log << "Loaded " << m_vecProteins.size() << " Proteins with " << m_iPeptidesNumberOfProteins << " peptides in experiment " << m_vecStrExperiments.at(iExperimentIndex)<< endl;
	return 1;

}

//bool CProteinWorker::mf_LoadProteins(std::ofstream& log, CTrainParam trainParam, int iExperimentIndex)
//{
//	cout << "\tLoading Protein!\n";
//	log << "\tLoading Protein!\n";
//	//open the file
//	FILE * pFile;
//	pFile = fopen(trainParam.m_strProteinsPath.c_str(), "r");
//	if (pFile == NULL)
//	{
//		cout << "Can not open " << trainParam.m_strProteinsPath << endl;
//		log << "Error:\tCan not open " << trainParam.m_strProteinsPath << endl;
//		exit(0);
//	}
//
//	//temporary variable
//	char Buffer[BUFFERLENGTH];
//	char *pstr;
//	char *pstr1;
//
//	map<string, int> mapAttributeAndColumns;
//	map<string, string>::iterator mapAttributeAndColumnsIter;
//	string strAttributeName;
//	int iAttributeColumn;
//	int iColumn;
//
//	fgets(Buffer, BUFFERLENGTH, pFile);
//	pstr = Buffer;
//	pstr1 = strstr(pstr, "\t");
//	iColumn = 0;
//	while (pstr1 != NULL)
//	{
//		*pstr1 = '\0';
//		strAttributeName = pstr;
//		mapAttributeAndColumns.insert(pair<string, int>(strAttributeName, iColumn));
//		iColumn++;
//		pstr = pstr1 + 1;
//		pstr1 = strstr(pstr, "\t");
//	}
//	strAttributeName = pstr;
//	mapAttributeAndColumns.insert(pair<string, int>(strAttributeName, iColumn));
//
//
//
//	CProtein ProteinTemp;
//	CPeptideAttribute peptideAttributeTemp;
//	PepAttributeWorker PeptideAttributeWorkerTemp;
//	string strWithFlankingRegionTemp;
//	string strtemp, substrtemp;
//	int begin, end;
//	map<string, string>::iterator mapProteinIDAndSequenceIter;
//	map<int, string>::iterator PeptidesIDAndSequenceIter;
//	int i = 0;
//	map<string, double> mapPeptideAndIntensitys;
//	map<string, double>::iterator mapPeptideAndIntensityIter;
//	int iPeptidesNumberOfOneProtein = 0;
//
//	fgets(Buffer, BUFFERLENGTH, pFile);    //jump the first row
//	fgets(Buffer, BUFFERLENGTH, pFile);
//	pstr = Buffer;
//
//	//load the information of amino acid 
//	m_vecAAindexAttributeHeaders = PeptideAttributeWorkerTemp.mf_LoadAAindex(log);
//	m_vecProteins.clear();
//	m_iPeptidesNumberOfProteins = 0;
//
//	int k = 0;
//	while (!feof(pFile))
//	{
//		//cout << "k= " << k << endl;
//		k++;
//		ProteinTemp.Clear();
//		mapPeptideAndIntensitys.clear();
//		pstr1 = strstr(pstr, "\t");
//		if (pstr1 != NULL)
//		{
//			*pstr1 = '\0';
//		}
//		else
//		{
//			log << "Error:\tThe format of proteins.txt is wrong.\n";
//			exit(0);
//		}
//		ProteinTemp.m_strProteinID = pstr;
//		pstr = pstr1 + 1;
//
//		mapProteinIDAndSequenceIter = m_mapProteinsIdAndSequence.find(ProteinTemp.m_strProteinID);
//		if (mapProteinIDAndSequenceIter != m_mapProteinsIdAndSequence.end())
//		{
//			ProteinTemp.m_strProteinSequence = mapProteinIDAndSequenceIter->second;
//		}
//		else
//		{
//			cout << "Can not find Protein " << ProteinTemp.m_strProteinID << endl;
//			log << "Error:\tCan not find Protein " << ProteinTemp.m_strProteinID << endl;
//			exit(0);
//		}
//
//		pstr1 = strstr(pstr, "\t");
//		if (pstr1 != NULL)
//		{
//			*pstr1 = '\0';
//		}
//		else
//		{
//			log << "Error:\tThe format of proteins.txt is wrong.\n";
//			exit(0);
//		}
//		ProteinTemp.m_strProteinFullName = pstr;
//		pstr = pstr1 + 1;
//
//		pstr1 = strstr(pstr, "\t");
//		pstr = pstr1 + 1;
//
//		pstr1 = strstr(pstr, "\t");
//		if (pstr1 != NULL)
//		{
//			*pstr1 = '\0';
//		}
//		else
//		{
//			log << "Error:\tThe format of proteins.txt is wrong.\n";
//			exit(0);
//		}
//		strtemp = pstr;
//		begin = 0;
//		end = strtemp.find(";", begin);
//
//
//		while (end != strtemp.npos)
//		{   //many peptides
//			substrtemp = strtemp.substr(begin, end - begin);
//			mapPeptideAndIntensitys.insert(pair<string, double>(substrtemp, 0.0));
//			begin = end + 1;
//			end = strtemp.find(";", begin);
//		}
//
//
//		for (i = 0; i < iExperimentIndex + 1; i++)
//		{
//			pstr = pstr1 + 1;
//			pstr1 = strstr(pstr, "\t");
//		}
//		if (pstr1 != NULL)
//		{
//			*pstr1 = '\0';
//		}
//		else
//		{
//			log << "Error:\tThe format of proteins.txt is wrong.\n";
//			exit(0);
//		}
//		strtemp = pstr;
//		begin = 0;
//		end = strtemp.find(";", begin);
//
//		mapPeptideAndIntensityIter = mapPeptideAndIntensitys.begin();
//		while (end != strtemp.npos)
//		{
//			substrtemp = strtemp.substr(begin, end - begin);
//			if (atof(substrtemp.c_str()) == 0.0)
//			{
//				mapPeptideAndIntensityIter = mapPeptideAndIntensitys.erase(mapPeptideAndIntensityIter);
//			}
//			else
//			{
//				mapPeptideAndIntensityIter->second = atof(substrtemp.c_str());
//				mapPeptideAndIntensityIter++;
//			}
//			begin = end + 1;
//			end = strtemp.find(";", begin);
//		}
//		iPeptidesNumberOfOneProtein = mapPeptideAndIntensitys.size();
//		m_iPeptidesNumberOfProteins += iPeptidesNumberOfOneProtein;
//		mapPeptideAndIntensityIter = mapPeptideAndIntensitys.begin();
//		for (; mapPeptideAndIntensityIter != mapPeptideAndIntensitys.end(); mapPeptideAndIntensityIter++)
//		{
//			strWithFlankingRegionTemp = ProteinTemp.mf_GetPeptidesAdjacentSequence(log, mapPeptideAndIntensityIter->first);
//			if (strWithFlankingRegionTemp == "NULL")
//			{
//				cout << "can not find peptide " << mapPeptideAndIntensityIter->first << "'s adjacent Sequence!\n";
//				log << "can not find peptide " << mapPeptideAndIntensityIter->first << "'s adjacent Sequence!\n";
//			}
//			else
//			{
//				ProteinTemp.m_vPeptidesSequences.push_back(mapPeptideAndIntensityIter->first);
//				ProteinTemp.m_vPeptidesSequencesWithFlankingRegion.push_back(strWithFlankingRegionTemp);                                         //根据肽段ID取得相应肽段序列和左右15位氨基酸
//				peptideAttributeTemp.Clear();
//				peptideAttributeTemp.m_vecAttributes = PeptideAttributeWorkerTemp.mf_GetAttributeFromSequence(mapPeptideAndIntensityIter->first, strWithFlankingRegionTemp);
//				ProteinTemp.m_vPeptidesAttributes.push_back(peptideAttributeTemp);
//				ProteinTemp.m_vdPeptidesMW.push_back(peptideAttributeTemp.m_vecAttributes.at(2));  //第三个属性是MW
//
//				ProteinTemp.m_vPeptidesNativeIntensity.push_back(mapPeptideAndIntensityIter->second);
//				ProteinTemp.m_vPeptidesAdjustIntensity.push_back(mapPeptideAndIntensityIter->second); //assign the same value with m_vPeptidesNativeIntensity for iterration
//				ProteinTemp.m_vAdjustPeptidesIntensityFactor.push_back(0.0);
//			}
//
//		}
//
//		for (i = 0; i < m_iNumberOfExperiments - iExperimentIndex - 1; i++)
//		{
//			pstr = pstr1 + 1;
//			pstr1 = strstr(pstr, "\t");
//		}
//
//		pstr = pstr1 + 1;
//		pstr1 = strstr(pstr, "\t");
//		if (pstr1 != NULL)
//		{
//			*pstr1 = '\0';
//		}
//		else
//		{
//			log << "Error:\tThe format of proteins.txt is wrong.\n";
//			exit(0);
//		}
//		m_bIfIBAQIntensityExist = atoi(pstr);
//		pstr = pstr1 + 1;
//
//		for (i = 0; i < iExperimentIndex + 1; i++)
//		{
//			pstr = pstr1 + 1;
//			pstr1 = strstr(pstr, "\t");
//		}
//		pstr1 = strstr(pstr, "\t");
//		if (pstr1 != NULL)
//		{
//			*pstr1 = '\0';
//		}
//		else
//		{
//			log << "Error:\tThe format of proteins.txt is wrong.\n";
//			exit(0);
//		}
//		ProteinTemp.m_dMaxquantIBAQ = atof(pstr);
//
//
//		for (i = 0; i < m_iNumberOfExperiments - iExperimentIndex - 1; i++)
//		{
//			pstr = pstr1 + 1;
//			pstr1 = strstr(pstr, "\t");
//		}
//
//		pstr = pstr1 + 1;
//		pstr1 = strstr(pstr, "\t");
//		if (pstr1 != NULL)
//		{
//			*pstr1 = '\0';
//		}
//		else
//		{
//			log << "Error:\tThe format of proteins.txt is wrong.\n";
//			exit(0);
//		}
//		m_bIFLFQIntensityExist = atof(pstr);
//		pstr = pstr1 + 1;
//
//		for (i = 0; i < iExperimentIndex; i++)
//		{
//			pstr = pstr1 + 1;
//			pstr1 = strstr(pstr, "\t");
//		}
//		if (m_bIFLFQIntensityExist&&iExperimentIndex>0)
//		{
//			pstr1 = strstr(pstr, "\t");
//			if (pstr1 != NULL)
//			{
//				*pstr1 = '\0';
//			}
//			else
//			{
//				log << "Error:\tThe format of proteins.txt is wrong.\n";
//				exit(0);
//			}
//			ProteinTemp.m_dLFQ = atof(pstr);
//			pstr = pstr1 + 1;
//		}
//		else
//		{
//			ProteinTemp.m_dLFQ = 0.0;
//		}
//
//		ProteinTemp.m_iPeptidesNumber = ProteinTemp.m_vPeptidesSequences.size();
//
//		//divide the median of peptides' intensities of one protein
//
//		//ProteinTemp.m_dPeptidesIntensityMedian = GetMedianIntensity(ProteinTemp.m_vPeptidesNativeIntensity);
//		//if (ProteinTemp.m_dPeptidesIntensityMedian != 0.0)
//		//{
//		//	for (size_t i = 0; i < ProteinTemp.m_vAdjustPeptidesIntensityFactor.size(); i++)
//		//	ProteinTemp.m_vAdjustPeptidesIntensityFactor.at(i) = ProteinTemp.m_vPeptidesNativeIntensity.at(i) / ProteinTemp.m_dPeptidesIntensityMedian;
//		//}	
//		//divide the maximum of peptides' intensities of one protein	
//		ProteinTemp.m_dPeptidesIntensityMax = GetMaxIntensity(ProteinTemp.m_vPeptidesNativeIntensity);
//		if (ProteinTemp.m_dPeptidesIntensityMax != 0.0)
//		{
//			for (size_t i = 0; i < ProteinTemp.m_vAdjustPeptidesIntensityFactor.size(); i++)
//				ProteinTemp.m_vAdjustPeptidesIntensityFactor.at(i) = ProteinTemp.m_vPeptidesNativeIntensity.at(i) / ProteinTemp.m_dPeptidesIntensityMax;
//		}
//
//		if (iPeptidesNumberOfOneProtein > 0)
//		{
//			m_vecProteins.push_back(ProteinTemp);
//		}
//
//		fgets(Buffer, BUFFERLENGTH, pFile);
//		pstr = Buffer;
//	}
//	cout << "Loaded " << m_vecProteins.size() << " Proteins with " << m_iPeptidesNumberOfProteins << " peptides in experiment " << m_vecStrExperiments.at(iExperimentIndex) << endl;
//	log << "Loaded " << m_vecProteins.size() << " Proteins with " << m_iPeptidesNumberOfProteins << " peptides in experiment " << m_vecStrExperiments.at(iExperimentIndex) << endl;
//	return 1;
//
//}

//load protein fasta 
bool CProteinWorker::mf_LoadProteinFasta(std::ofstream& log, string strProteinFasterFilePath, FastaType fastatype)
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
	while ((Buffer[0] == '\n') && (!feof(pFile)))  //允许有空行
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
			cout << "Can not sparse the fasta file by the regular expression.\n";
			log << "Can not sparse the fasta file by the regular expression.\n";
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
			m_mapProteinsIdAndSequence.insert(pair<string, string>(strIDTemp, strSequenceTemp));
			strSequenceTemp.clear();
			strFirstRow = Buffer;
			valid = std::regex_search(strFirstRow, Matchresult, fastatype);
			if (valid == true)
			{
				strIDTemp = Matchresult[1];
			}
			else
			{
				cout << "Can not sparse the fasta file by the regular expression.\n";
				log << "Can not sparse the fasta file by the regular expression.\n";
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
	m_mapProteinsIdAndSequence.insert(pair<string, string>(strIDTemp, strSequenceTemp));

	fclose(pFile);
	return 1;
}

//calculate the IBAQ and SPEAQInPeptides by calling mf_CalculateProteinIBAQ and mf_CalculateProteinSPEAQInPeptides
bool CProteinWorker::mf_CalculateProteinIntensity(std::ofstream& log, int iExperimentIndex)
{
	mf_CalculateProteinIBAQ(log);
	mf_CalculateProteinSPEAQInPeptides(log, iExperimentIndex);
	return true;
}

// recalculate the IBAQ according to the protein' peptides intensity;
void  CProteinWorker::mf_CalculateProteinIBAQ(std::ofstream& log)
{
	CProteinDigestion proteinDigestion(m_predictParam);
	for (int i = 0; i < m_vecProteins.size(); i++)
	{
		for (size_t j = 0; j < m_vecProteins.at(i).m_vPeptidesNativeIntensity.size(); j++)
		{
			m_vecProteins.at(i).m_dPeptidesIntensitySum += m_vecProteins.at(i).m_vPeptidesNativeIntensity.at(j);
		}
		m_vecProteins.at(i).m_iPeptidesNumber = m_vecProteins.at(i).m_vPeptidesSequences.size();
		// enzyme number
		m_vecProteins.at(i).m_iNumberOfTheoreticEnzyme = proteinDigestion.mf_CalculateOptDigestionNumber(m_vecProteins.at(i).m_strProteinSequence);
		if (m_vecProteins.at(i).m_iNumberOfTheoreticEnzyme>0)
			m_vecProteins.at(i).m_dReCalculateIBAQ = m_vecProteins.at(i).m_dPeptidesIntensitySum / m_vecProteins.at(i).m_iNumberOfTheoreticEnzyme;
		else
		{
			cerr << "Protein " << m_vecProteins.at(i).m_strProteinID << " do not have Digension peptides." << endl;
			log << "Error:\tProtein " << m_vecProteins.at(i).m_strProteinID << " do not have Digension peptides." << endl;
			exit(0);
		}
	}
	
}
// recalculate the SPEAQInPeptides according to the protein' peptides intensity;
void CProteinWorker::mf_CalculateProteinSPEAQInPeptides(std::ofstream& log, int iExperimentIndex)
{
	CRegression regression;
	regression.mf_RegressionRun(log,m_vecProteins, m_trainParam, m_predictParam,m_vecStrExperiments.at(iExperimentIndex));

	//Load the correction correction factor calculated by other methods
	// maxquant 1..5.2.8 NotDeleteSharePeptide
	// UPS2_mouse	divided by the median 	stepwise realized by Matlab	
	 //string strPeptidesFactorFromOtherMethodPath="E:\\data\\BPRC\\regressionResult\\maxquant1.5.2.8\\NotDeleteSharePeptide\\UPS2_mouse\\DivideMid\\stepwiseMatlabVersion\\PredictPeptidesFactor.txt";

	//UPS2_mouse	divided by the Maximun	stepwise realized by Matlab	
	//string strPeptidesFactorFromOtherMethodPath="E:\\data\\BPRC\\regressionResult\\maxquant1.5.2.8\\NotDeleteSharePeptide\\UPS2_mouse\\DivideMax\\stepwiseMatlabVersion\\PredictPeptidesFactor.txt";
	// 
	// //UPS2_yeast	divided by the median 	stepwise realized by Matlab	
	 //string strPeptidesFactorFromOtherMethodPath="E:\\data\\BPRC\\regressionResult\\maxquant1.5.2.8\\NotDeleteSharePeptide\\UPS2_yeast\\DivideMid\\stepwiseMatlabVersion\\PredictPeptidesFactor.txt";
	// 
	// //UPS2_yeast	divided by the Maximun	stepwise realized by Matlab	
	 //string strPeptidesFactorFromOtherMethodPath="E:\\data\\BPRC\\regressionResult\\maxquant1.5.2.8\\NotDeleteSharePeptide\\UPS2_yeast\\DivideMax\\stepwiseMatlabVersion\\PredictPeptidesFactor.txt";
	// 
	// //UPS2_only	divided by the median 	stepwise realized by Matlab	
	//string strPeptidesFactorFromOtherMethodPath = "E:\\data\\BPRC\\regressionResult\\maxquant1.5.2.8\\NotDeleteSharePeptide\\UPS2_only\\DivideMid\\stepwiseMatlabVersion\\PredictPeptidesFactor.txt";
	// 
	// //UPS2_only	divided by the Maximun	stepwise realized by Matlab	
	//string strPeptidesFactorFromOtherMethodPath = "E:\\data\\BPRC\\regressionResult\\maxquant1.5.2.8\\NotDeleteSharePeptide\\UPS2_only\\DivideMax\\stepwiseMatlabVersion\\PredictPeptidesFactor.txt";
	

	//UPS2_mouse	divided by the median 	BART realized by R
	//string strPeptidesFactorFromOtherMethodPath = "E:\\data\\BPRC\\regressionResult\\maxquant1.5.2.8\\NotDeleteSharePeptide\\UPS2_mouse\\DivideMid\\BARTRVersion\\BARTPredictY.txt";

	// UPS2_mouse	divided by the Maximun	BART realized by R	
	//string strPeptidesFactorFromOtherMethodPath ="E:\\data\\BPRC\\regressionResult\\maxquant1.5.2.8\\NotDeleteSharePeptide\\UPS2_mouse\\DivideMax\\BARTRVersion\\BARTPredictY.txt";
	//
	//%UPS2_yeast	divided by the median 	BART realized by R	
	//string strPeptidesFactorFromOtherMethodPath ="E:\\data\\BPRC\\regressionResult\\maxquant1.5.2.8\\NotDeleteSharePeptide\\UPS2_yeast\\DivideMid\\BARTRVersion\\BARTPredictY.txt";
	//
	//%UPS2_yeast	divided by the Maximun	BART realized by R	
	//string strPeptidesFactorFromOtherMethodPath ="E:\\data\\BPRC\\regressionResult\\maxquant1.5.2.8\\NotDeleteSharePeptide\\UPS2_yeast\\DivideMax\\BARTRVersion\\BARTPredictY.txt";
	//
	//UPS2_only	divided by the median 	BART realized by R	
	//string strPeptidesFactorFromOtherMethodPath ="E:\\data\\BPRC\\regressionResult\\maxquant1.5.2.8\\NotDeleteSharePeptide\\UPS2_only\\DivideMid\\BARTRVersion\\BARTPredictY.txt";

	//UPS2_only	divided by the Maximun	BART realized by R	
	//string strPeptidesFactorFromOtherMethodPath = "E:\\data\\BPRC\\regressionResult\\maxquant1.5.2.8\\NotDeleteSharePeptide\\UPS2_only\\DivideMax\\BARTRVersion\\BARTPredictY.txt";

	//mf_LoadPeptidesFactorFromOtherMethod(strPeptidesFactorFromOtherMethodPath, m_vecProteins);

//maxquant 1..5.2.8  DeleteSharePeptides ***********************************************************
	// UPS2_mouse	divided by the median 	stepwise realized by Matlab	
	//string strPeptidesFactorFromOtherMethodPath="E:\\data\\BPRC\\regressionResult\\maxquant1.5.2.8\\DeleteSharePeptides\\UPS2_mouse\\DivideMid\\stepwiseMatlabVersion\\PredictPeptidesFactor.txt";

	//UPS2_mouse	divided by the Maximun	stepwise realized by Matlab	
	//string strPeptidesFactorFromOtherMethodPath="E:\\data\\BPRC\\regressionResult\\maxquant1.5.2.8\\DeleteSharePeptides\\UPS2_mouse\\DivideMax\\stepwiseMatlabVersion\\PredictPeptidesFactor.txt";
	// 
	// //UPS2_yeast	divided by the median 	stepwise realized by Matlab	
	//string strPeptidesFactorFromOtherMethodPath="E:\\data\\BPRC\\regressionResult\\maxquant1.5.2.8\\DeleteSharePeptides\\UPS2_yeast\\DivideMid\\stepwiseMatlabVersion\\PredictPeptidesFactor.txt";
	// 
	// //UPS2_yeast	divided by the Maximun	stepwise realized by Matlab	
	//string strPeptidesFactorFromOtherMethodPath="E:\\data\\BPRC\\regressionResult\\maxquant1.5.2.8\\DeleteSharePeptides\\UPS2_yeast\\DivideMax\\stepwiseMatlabVersion\\PredictPeptidesFactor.txt";
	// 
	// //UPS2_only	divided by the median 	stepwise realized by Matlab	
	//string strPeptidesFactorFromOtherMethodPath = "E:\\data\\BPRC\\regressionResult\\maxquant1.5.2.8\\DeleteSharePeptides\\UPS2_only\\DivideMid\\stepwiseMatlabVersion\\PredictPeptidesFactor.txt";
	// 
	// //UPS2_only	divided by the Maximun	stepwise realized by Matlab	
	//string strPeptidesFactorFromOtherMethodPath = "E:\\data\\BPRC\\regressionResult\\maxquant1.5.2.8\\DeleteSharePeptides\\UPS2_only\\DivideMax\\stepwiseMatlabVersion\\PredictPeptidesFactor.txt";


	//UPS2_mouse	divided by the median 	BART realized by R
	//string strPeptidesFactorFromOtherMethodPath = "E:\\data\\BPRC\\regressionResult\\maxquant1.5.2.8\\DeleteSharePeptides\\UPS2_mouse\\DivideMid\\BARTRVersion\\BARTPredictY.txt";

	// UPS2_mouse	divided by the Maximun	BART realized by R	
	//string strPeptidesFactorFromOtherMethodPath ="E:\\data\\BPRC\\regressionResult\\maxquant1.5.2.8\\DeleteSharePeptides\\UPS2_mouse\\DivideMax\\BARTRVersion\\BARTPredictY.txt";
	//
	//%UPS2_yeast	divided by the median 	BART realized by R	
	//string strPeptidesFactorFromOtherMethodPath ="E:\\data\\BPRC\\regressionResult\\maxquant1.5.2.8\\DeleteSharePeptides\\UPS2_yeast\\DivideMid\\BARTRVersion\\BARTPredictY.txt";
	//
	//%UPS2_yeast	divided by the Maximun	BART realized by R	
	//string strPeptidesFactorFromOtherMethodPath ="E:\\data\\BPRC\\regressionResult\\maxquant1.5.2.8\\DeleteSharePeptides\\UPS2_yeast\\DivideMax\\BARTRVersion\\BARTPredictY.txt";
	//
	//UPS2_only	divided by the median 	BART realized by R	
	//string strPeptidesFactorFromOtherMethodPath ="E:\\data\\BPRC\\regressionResult\\maxquant1.5.2.8\\DeleteSharePeptides\\UPS2_only\\DivideMid\\BARTRVersion\\BARTPredictY.txt";

	//UPS2_only	divided by the Maximun	BART realized by R	
	//string strPeptidesFactorFromOtherMethodPath = "E:\\data\\BPRC\\regressionResult\\maxquant1.5.2.8\\DeleteSharePeptides\\UPS2_only\\DivideMax\\BARTRVersion\\BARTPredictY.txt";

//maxquant 1.5.3.8  DeleteSharePeptides ***********************************************************
	// UPS2_mouse	divided by the median 	stepwise realized by Matlab	
	//string strPeptidesFactorFromOtherMethodPath="E:\\data\\BPRC\\regressionResult\\maxquant1.5.3.8\\DeleteSharePeptides\\UPS2_mouse\\DivideMid\\stepwiseMatlabVersion\\PredictPeptidesFactor.txt";

	//UPS2_mouse	divided by the Maximun	stepwise realized by Matlab	
	//string strPeptidesFactorFromOtherMethodPath="E:\\data\\BPRC\\regressionResult\\maxquant1.5.3.8\\DeleteSharePeptides\\UPS2_mouse\\DivideMax\\stepwiseMatlabVersion\\PredictPeptidesFactor.txt";
	// 
	// //UPS2_yeast	divided by the median 	stepwise realized by Matlab	
	//string strPeptidesFactorFromOtherMethodPath="E:\\data\\BPRC\\regressionResult\\maxquant1.5.3.8\\DeleteSharePeptides\\UPS2_yeast\\DivideMid\\stepwiseMatlabVersion\\PredictPeptidesFactor.txt";
	// 
	// //UPS2_yeast	divided by the Maximun	stepwise realized by Matlab	
	//string strPeptidesFactorFromOtherMethodPath="E:\\data\\BPRC\\regressionResult\\maxquant1.5.3.8\\DeleteSharePeptides\\UPS2_yeast\\DivideMax\\stepwiseMatlabVersion\\PredictPeptidesFactor.txt";
	// 
	// //UPS2_only	divided by the median 	stepwise realized by Matlab	
	//string strPeptidesFactorFromOtherMethodPath = "E:\\data\\BPRC\\regressionResult\\maxquant1.5.3.8\\DeleteSharePeptides\\UPS2_only\\DivideMid\\stepwiseMatlabVersion\\PredictPeptidesFactor.txt";
	// 
	// //UPS2_only	divided by the Maximun	stepwise realized by Matlab	
	//string strPeptidesFactorFromOtherMethodPath = "E:\\data\\BPRC\\regressionResult\\maxquant1.5.3.8\\DeleteSharePeptides\\UPS2_only\\DivideMax\\stepwiseMatlabVersion\\PredictPeptidesFactor.txt";


	//UPS2_mouse	divided by the median 	BART realized by R
	//string strPeptidesFactorFromOtherMethodPath = "E:\\data\\BPRC\\regressionResult\\maxquant1.5.3.8\\DeleteSharePeptides\\UPS2_mouse\\DivideMid\\BARTRVersion\\BARTPredictY.txt";

	// UPS2_mouse	divided by the Maximun	BART realized by R	
	//string strPeptidesFactorFromOtherMethodPath ="E:\\data\\BPRC\\regressionResult\\maxquant1.5.3.8\\DeleteSharePeptides\\UPS2_mouse\\DivideMax\\BARTRVersion\\BARTPredictY.txt";
	//
	//%UPS2_yeast	divided by the median 	BART realized by R	
	//string strPeptidesFactorFromOtherMethodPath ="E:\\data\\BPRC\\regressionResult\\maxquant1.5.3.8\\DeleteSharePeptides\\UPS2_yeast\\DivideMid\\BARTRVersion\\BARTPredictY.txt";
	//
	//%UPS2_yeast	divided by the Maximun	BART realized by R	
	//string strPeptidesFactorFromOtherMethodPath ="E:\\data\\BPRC\\regressionResult\\maxquant1.5.3.8\\DeleteSharePeptides\\UPS2_yeast\\DivideMax\\BARTRVersion\\BARTPredictY.txt";
	//
	//UPS2_only	divided by the median 	BART realized by R	
	//string strPeptidesFactorFromOtherMethodPath ="E:\\data\\BPRC\\regressionResult\\maxquant1.5.3.8\\DeleteSharePeptides\\UPS2_only\\DivideMid\\BARTRVersion\\BARTPredictY.txt";

	//UPS2_only	divided by the Maximun	BART realized by R	
	//string strPeptidesFactorFromOtherMethodPath = "E:\\data\\BPRC\\regressionResult\\maxquant1.5.3.8\\DeleteSharePeptides\\UPS2_only\\DivideMax\\BARTRVersion\\BARTPredictY.txt";

// maxquant 1.5.3.8 NotDeleteSharePeptide
// UPS2_mouse	divided by the median 	stepwise realized by Matlab	
//string strPeptidesFactorFromOtherMethodPath="E:\\data\\BPRC\\regressionResult\\maxquant1.5.3.8\\NotDeleteSharePeptide\\UPS2_mouse\\DivideMid\\stepwiseMatlabVersion\\PredictPeptidesFactor.txt";

//UPS2_mouse	divided by the Maximun	stepwise realized by Matlab	
//string strPeptidesFactorFromOtherMethodPath="E:\\data\\BPRC\\regressionResult\\maxquant1.5.3.8\\NotDeleteSharePeptide\\UPS2_mouse\\DivideMax\\stepwiseMatlabVersion\\PredictPeptidesFactor.txt";
// 
// //UPS2_yeast	divided by the median 	stepwise realized by Matlab	
//string strPeptidesFactorFromOtherMethodPath="E:\\data\\BPRC\\regressionResult\\maxquant1.5.3.8\\NotDeleteSharePeptide\\UPS2_yeast\\DivideMid\\stepwiseMatlabVersion\\PredictPeptidesFactor.txt";
// 
// //UPS2_yeast	divided by the Maximun	stepwise realized by Matlab	
//string strPeptidesFactorFromOtherMethodPath="E:\\data\\BPRC\\regressionResult\\maxquant1.5.3.8\\NotDeleteSharePeptide\\UPS2_yeast\\DivideMax\\stepwiseMatlabVersion\\PredictPeptidesFactor.txt";
// 
// //UPS2_only	divided by the median 	stepwise realized by Matlab	
//string strPeptidesFactorFromOtherMethodPath = "E:\\data\\BPRC\\regressionResult\\maxquant1.5.3.8\\NotDeleteSharePeptide\\UPS2_only\\DivideMid\\stepwiseMatlabVersion\\PredictPeptidesFactor.txt";
// 
// //UPS2_only	divided by the Maximun	stepwise realized by Matlab	
//string strPeptidesFactorFromOtherMethodPath = "E:\\data\\BPRC\\regressionResult\\maxquant1.5.3.8\\NotDeleteSharePeptide\\UPS2_only\\DivideMax\\stepwiseMatlabVersion\\PredictPeptidesFactor.txt";


//UPS2_mouse	divided by the median 	BART realized by R
//string strPeptidesFactorFromOtherMethodPath = "E:\\data\\BPRC\\regressionResult\\maxquant1.5.3.8\\NotDeleteSharePeptide\\UPS2_mouse\\DivideMid\\BARTRVersion\\BARTPredictY.txt";

// UPS2_mouse	divided by the Maximun	BART realized by R	
//string strPeptidesFactorFromOtherMethodPath ="E:\\data\\BPRC\\regressionResult\\maxquant1.5.3.8\\NotDeleteSharePeptide\\UPS2_mouse\\DivideMax\\BARTRVersion\\BARTPredictY.txt";
//
//%UPS2_yeast	divided by the median 	BART realized by R	
//string strPeptidesFactorFromOtherMethodPath ="E:\\data\\BPRC\\regressionResult\\maxquant1.5.3.8\\NotDeleteSharePeptide\\UPS2_yeast\\DivideMid\\BARTRVersion\\BARTPredictY.txt";
//
//%UPS2_yeast	divided by the Maximun	BART realized by R	
//string strPeptidesFactorFromOtherMethodPath ="E:\\data\\BPRC\\regressionResult\\maxquant1.5.3.8\\NotDeleteSharePeptide\\UPS2_yeast\\DivideMax\\BARTRVersion\\BARTPredictY.txt";
//
//UPS2_only	divided by the median 	BART realized by R	
//string strPeptidesFactorFromOtherMethodPath ="E:\\data\\BPRC\\regressionResult\\maxquant1.5.3.8\\NotDeleteSharePeptide\\UPS2_only\\DivideMid\\BARTRVersion\\BARTPredictY.txt";

//UPS2_only	divided by the Maximun	BART realized by R	
//string strPeptidesFactorFromOtherMethodPath = "E:\\data\\BPRC\\regressionResult\\maxquant1.5.3.8\\NotDeleteSharePeptide\\UPS2_only\\DivideMax\\BARTRVersion\\BARTPredictY.txt";

	vector<CProtein>::iterator proteinIter;
	vector<double>::iterator peptidesIter;
	vector<double>::iterator PeptidesCorrectionFactorIter;
	vector<double> AdjustPeptidesIntensity;
	vector<double> NativePeptidesIntensity;
	double dPeptideIntensityTemp;

	int i = 0;
	for (proteinIter = m_vecProteins.begin(); proteinIter != m_vecProteins.end(); proteinIter++)
	{
		i = 0;
		if (proteinIter->m_iPeptidesNumber >= UniquePeptidesTestThreshold)
		{
			proteinIter->m_bIfCalculateSPEAQInPeptides = true;
			if (proteinIter->m_iPeptidesNumber > UniquePeptidesCorrectThreshold)
			{
				proteinIter->m_dProteinSPEAQInPeptides = 0.0;
				proteinIter->m_dProteinSPEAQInProtein = 0.0;
				proteinIter->m_vPeptidesAdjustIntensity.clear();
				for (peptidesIter = proteinIter->m_vPeptidesNativeIntensity.begin(); peptidesIter != proteinIter->m_vPeptidesNativeIntensity.end(); peptidesIter++) //use vPeptidesNativeIntensity instead of vPeptidesIntensity! by CC 20150527
				{
					dPeptideIntensityTemp = *peptidesIter;
					proteinIter->m_dProteinSPEAQInProtein += dPeptideIntensityTemp;
					NativePeptidesIntensity.push_back(dPeptideIntensityTemp / proteinIter->m_vdPeptidesMW.at(i));   //change by gzhq 20151013

					if (proteinIter->m_vPeptidesCorrectionFactor.at(i) != 0.0)
						dPeptideIntensityTemp = dPeptideIntensityTemp / proteinIter->m_vPeptidesCorrectionFactor.at(i);
					else
						dPeptideIntensityTemp = 0.0;
					i++;
					AdjustPeptidesIntensity.push_back(dPeptideIntensityTemp);
					proteinIter->m_vPeptidesAdjustIntensity.push_back(dPeptideIntensityTemp);
					proteinIter->m_dProteinSPEAQInPeptides += dPeptideIntensityTemp;
				}
			}
			else
			{ //not enough peptides for correction
				proteinIter->m_dProteinSPEAQInPeptides = 0.0;
				proteinIter->m_dProteinSPEAQInProtein = 0.0;
				proteinIter->m_vPeptidesAdjustIntensity.clear();
				for (peptidesIter = proteinIter->m_vPeptidesNativeIntensity.begin(); peptidesIter != proteinIter->m_vPeptidesNativeIntensity.end(); peptidesIter++) //use vPeptidesNativeIntensity instead of vPeptidesIntensity! by CC 20150527
				{
					dPeptideIntensityTemp = *peptidesIter;
					proteinIter->m_dProteinSPEAQInProtein += dPeptideIntensityTemp;
					NativePeptidesIntensity.push_back(dPeptideIntensityTemp);
					AdjustPeptidesIntensity.push_back(dPeptideIntensityTemp);
					proteinIter->m_vPeptidesAdjustIntensity.push_back(dPeptideIntensityTemp);
					proteinIter->m_dProteinSPEAQInPeptides += dPeptideIntensityTemp;
				}
			}
			proteinIter->m_dCVNative = CalculateCV(NativePeptidesIntensity);
			proteinIter->m_dCVAfterAdjustAfter = CalculateCV(AdjustPeptidesIntensity);
			NativePeptidesIntensity.clear();
			AdjustPeptidesIntensity.clear();
			if (proteinIter->m_iPeptidesNumber != 0)
			{
				proteinIter->m_dProteinSPEAQInPeptides = proteinIter->m_dProteinSPEAQInPeptides / proteinIter->m_iPeptidesNumber;
				if (proteinIter->m_dLFQ != 0.0)
				{
					proteinIter->m_dLFQ = proteinIter->m_dLFQ / proteinIter->m_iPeptidesNumber;
				}
				PeptidesCorrectionFactorIter = proteinIter->m_vPeptidesCorrectionFactor.begin();
				proteinIter->m_dProteinCorrectionFactor = 0.0;
				for (; PeptidesCorrectionFactorIter != proteinIter->m_vPeptidesCorrectionFactor.end(); PeptidesCorrectionFactorIter++)
				{
					proteinIter->m_dProteinCorrectionFactor += *PeptidesCorrectionFactorIter;
				}

				proteinIter->m_dProteinSPEAQInProtein = proteinIter->m_dProteinSPEAQInProtein / proteinIter->m_dProteinCorrectionFactor;
			}
			else
			{
				proteinIter->m_dProteinSPEAQInPeptides = 0.0;
				proteinIter->m_dProteinSPEAQInProtein = 0.0;
				cout << "Protein " << proteinIter->m_strProteinID << " do not have unique peptides. ";
				cout << "So we make its SPEAQInPeptides and ProteinSPEAQInProtein  0\n";
				log << "Protein " << proteinIter->m_strProteinID << " do not have unique peptides. ";
				log << "So we make its SPEAQInPeptides and ProteinSPEAQInProtein 0\n";
			}
		}		
	}
}

void CProteinWorker::mf_SaveProteinsNativePeptides(std::ofstream& log, string path)
{
	ofstream ofile(path.c_str());
	if (!ofile)
	{
		cout << "can not open " << path << endl;
		log << "Error:\tCan not open " << path << endl;
		exit(0);
	}
	else
	{
		cout << "save Proteinpeptides in " << path << endl;
		log << "save Proteinpeptides in " << path << endl;
	}
		
	ofile << "PeptidesSequence\tProteinID\tPeptideCorrectionFactor\tpeptideNativeIntensity\tPeptidesAttributes\n";
	vector<CProtein>::iterator  ProteinIter;
	int i = 0;
	for (ProteinIter = m_vecProteins.begin(); ProteinIter != m_vecProteins.end(); ProteinIter++)
	{
		i = 0;

		for (i = 0; i < ProteinIter->m_iPeptidesNumber; i++)
		{
			
			ofile << ProteinIter->m_vPeptidesSequences.at(i) << "\t" ;
			ofile << ProteinIter->m_strProteinID << "\t";
			ofile << ProteinIter->m_vPeptidesCorrectionFactor.at(i) << "\t";
			ofile << ProteinIter->m_vPeptidesNativeIntensity.at(i) << "\t";
			//ofile << ProteinIter->m_dPeptidesIntensityMedian << "\t";
			ofile << ProteinIter->m_vPeptidesAttributes.at(i) << endl;
			//ofile << ProteinIter->m_dPeptidesIntensityMax << "\n";
		}

	}
	ofile.close();
}

void CProteinWorker::mf_saveProteinSPEAQInPeptides(std::ofstream& log, string ProteinSPEAQInPeptidesResultPath)
{
	ofstream ofile;
	ofile.open(ProteinSPEAQInPeptidesResultPath.c_str());
	if (!ofile)
	{
		cout << "Can not open " << ProteinSPEAQInPeptidesResultPath << endl;
		log << "Error:\tCan not open " << ProteinSPEAQInPeptidesResultPath << endl;

		exit(0);
	}
	else
	{
		cout << "Begin Save " << m_vecProteins.size() << " Proteins' SPEAQInPeptides in " << ProteinSPEAQInPeptidesResultPath << endl;
		log << "Begin Save " << m_vecProteins.size() << " Proteins' SPEAQInPeptides in " << ProteinSPEAQInPeptidesResultPath << endl;
	}

	vector<CProtein>::iterator proteinIter;
	vector<double>::iterator peptidesIter;
	if (m_bIfIBAQIntensityExist)
	{
		ofile << "Protein ID\tMajority protein IDs\tNumberOfUniquePeptides\tMaxquantIBAQ\t";
		ofile << "IfCalculateSPEAQInPeptides\tMyIBAQ\tLFQ\tSPEAQInPeptides\tdCVNative\tdCVAfterAdjust\t";
		ofile << "NumberOfTheoreticEnzyme\tm_dProteinCorrectionFactor\tm_dProteinQsocre2\t";
		ofile << "m_vPeptidesSequences\tm_vPeptidesNativeIntensity\tm_vPeptidesCorrectionFactor\t";
		ofile<<"m_vAdjustPeptidesIntensityFactor\tPeptidesAdjustIntensity\n";

	}
	int i = 0;
	for (proteinIter = m_vecProteins.begin(); proteinIter != m_vecProteins.end(); proteinIter++)
	{

		ofile << proteinIter->m_strProteinID << "\t" << proteinIter->m_strProteinFullName << "\t";
		ofile << proteinIter->m_iPeptidesNumber << "\t";
		ofile << proteinIter->m_dMaxquantIBAQ << "\t";
		ofile << proteinIter->m_bIfCalculateSPEAQInPeptides << "\t";
		ofile << proteinIter->m_dReCalculateIBAQ << "\t";
		ofile << proteinIter->m_dLFQ << "\t";
		ofile << proteinIter->m_dProteinSPEAQInPeptides << "\t" << proteinIter->m_dCVNative << "\t" << proteinIter->m_dCVAfterAdjustAfter << "\t";
		ofile << proteinIter->m_iNumberOfTheoreticEnzyme << "\t";
		ofile << proteinIter->m_dProteinCorrectionFactor << "\t";
		ofile << proteinIter->m_dProteinSPEAQInProtein << "\t";
		for (i = 0; i < proteinIter->m_vPeptidesSequences.size(); i++)
		{
			ofile << proteinIter->m_vPeptidesSequences.at(i) << ";";
		}
		ofile << "\t";

		for (i = 0; i < proteinIter->m_vPeptidesNativeIntensity.size() ; i++)
		{
			ofile << proteinIter->m_vPeptidesNativeIntensity.at(i) << ";";
		}
		ofile  << "\t";

		if (proteinIter->m_bIfCalculateSPEAQInPeptides)
		{
			for (i = 0; i < proteinIter->m_vPeptidesCorrectionFactor.size() ; i++)
			{
				ofile << proteinIter->m_vPeptidesCorrectionFactor.at(i) << ";";
			}
			ofile  << "\t";

			for (i = 0; i < proteinIter->m_vAdjustPeptidesIntensityFactor.size() ; i++)
			{
				ofile << proteinIter->m_vAdjustPeptidesIntensityFactor.at(i) << ";";
			}
			ofile  << "\t";

			for (i = 0; i < proteinIter->m_vPeptidesAdjustIntensity.size(); i++)
			{
				ofile << proteinIter->m_vPeptidesAdjustIntensity.at(i) << ";";
			}
			ofile  << "\n";
		}
		else
		{
			ofile << "\n";
		}

	}
	ofile.close();

}

bool CProteinWorker::mf_LoadProteinSPEAQInPeptides(std::ofstream& log, string path)
{
	cout << "\tLoading Protein!\n";
	log << "\tLoading Protein!\n";
	//open the file
	FILE * pFile;
	pFile = fopen(path.c_str(), "r");
	if (pFile == NULL)
	{
		cout << "Can not open " << path << endl;
		log << "Error:\tCan not open " << path << endl;
		exit(0);
		return false;
	}

	//temporary variable
	char Buffer[BUFFERLENGTH];
	char *pstr;
	char *pstr1;
	CProtein proteinTemp;


	fgets(Buffer, BUFFERLENGTH, pFile);    //jump the first row
	fgets(Buffer, BUFFERLENGTH, pFile);
	pstr = Buffer;

	while (!feof(pFile))
	{
		proteinTemp.Clear();
		pstr1 = strstr(pstr, "\t");
		if (pstr1 != NULL)
		{
			*pstr1 = '\0';
			proteinTemp.m_strProteinID = pstr;
			pstr = pstr1 + 1;
		}
		else
		{
			cout << "The format in " << path << " is wrong!\n";
			log << "Error:\tThe format in " << path << " is wrong!\n";
			exit(0);
		}

		pstr1 = strstr(pstr, "\t");
		if (pstr1 != NULL)
		{
			pstr = pstr1 + 1;
		}
		else
		{
			cout << "The format in " << path << " is wrong!\n";
			log << "Error:\tThe format in " << path << " is wrong!\n";
			exit(0);
		}

		pstr1 = strstr(pstr, "\t");
		if (pstr1 != NULL)
		{
			pstr = pstr1 + 1;
		}
		else
		{
			cout << "The format in " << path << " is wrong!\n";
			log << "Error:\tThe format in " << path << " is wrong!\n";
			exit(0);
		}

		pstr1 = strstr(pstr, "\t");
		if (pstr1 != NULL)
		{
			*pstr1 = '\0';
			proteinTemp.m_dMaxquantIBAQ = atof(pstr);
			pstr = pstr1 + 1;
		}
		else
		{
			cout << "The format in " << path << " is wrong!\n";
			log << "Error:\tThe format in " << path << " is wrong!\n";
			exit(0);
		}

		pstr1 = strstr(pstr, "\t");
		if (pstr1 != NULL)
		{
			pstr = pstr1 + 1;
		}
		else
		{
			cout << "The format in  " << path << " is wrong!\n";
			log << "Error:\tThe format in  " << path << " is wrong!\n";
			exit(0);
		}

		pstr1 = strstr(pstr, "\t");
		if (pstr1 != NULL)
		{
			pstr = pstr1 + 1;
		}
		else
		{
			cout << "The format in " << path << " is wrong!\n";
			log << "Error:\tThe format in " << path << " is wrong!\n";
			exit(0);
		}

		pstr1 = strstr(pstr, "\t");
		if (pstr1 != NULL)
		{
			*pstr1 = '\0';
			proteinTemp.m_dLFQ = atof(pstr);
			pstr = pstr1 + 1;
		}
		else
		{
			cout << "The format in " << path << " is wrong!\n";
			log << "Error:\tThe format in " << path << " is wrong!\n";
			exit(0);
		}
		pstr1 = strstr(pstr, "\t");
		if (pstr1 != NULL)
		{
			*pstr1 = '\0';
			proteinTemp.m_dProteinSPEAQInPeptides = atof(pstr);
			pstr = pstr1 + 1;
		}
		else
		{
			cout << "The format in " << path << " is wrong!\n";
			log << "Error:\tThe format in " << path << " is wrong!\n";
			exit(0);
		}

		for (int t = 0; t < 4; t++)
		{
			pstr1 = strstr(pstr, "\t");
			if (pstr1 != NULL)
			{
				pstr = pstr1 + 1;
			}
			else
			{
				cout << "The format in " << path << " is wrong!\n";
				log << "Error:\tThe format in " << path << " is wrong!\n";
				exit(0);
			}
		}

		pstr1 = strstr(pstr, "\t");
		if (pstr1 != NULL)
		{
			*pstr1 = '\0';
			proteinTemp.m_dProteinSPEAQInProtein = atof(pstr);
			pstr = pstr1 + 1;
		}
		else
		{
			cout << "The format in " << path << " is wrong!\n";
			log << "Error:\tThe format in " << path << " is wrong!\n";
			exit(0);
		}
		m_vecProteins.push_back(proteinTemp);
		fgets(Buffer, BUFFERLENGTH, pFile);
		pstr = Buffer;
	}
}

// merge all ProteinSPEAQInPeptides in different experiment;
void CProteinWorker::mf_MergeProteinSPEAQInPeptides(std::ofstream& log, string ProteinSPEAQInPeptidesResultPath)
{
	log << "Begin merge protein PepSPEAQ!\n";
	string strProteinSPEAQInPeptidesPath;
	FILE * Ofiles[10]; // at most ten experiments;
	char Buffer[BUFFERLENGTH];
	char *pstr;
	char *pstr1;
	string strIDTemp;
	string strTemp;
	

	CMergedProteins MergedProteins;
	CMergedProtein mergedProteinTemp;
	// assuming that every SPEAQInPeptides result is sorted, our result is OK;
	for (int i = 0; i < m_vecStrExperiments.size(); i++)
	{
		m_vecProteins.clear();
		strProteinSPEAQInPeptidesPath = m_predictParam.m_strProteinSPEAQInPeptidesPath + m_vecStrExperiments.at(i) + ".txt";
		mf_LoadProteinSPEAQInPeptides(log,strProteinSPEAQInPeptidesPath);
		if (MergedProteins.m_vecMergedProteins.size() == 0)
		{
			for (int j = 0; j < m_vecProteins.size(); j++)
			{
				mergedProteinTemp.Clear();
				mergedProteinTemp.m_strPrtoteinName = m_vecProteins[j].m_strProteinID;
				mergedProteinTemp.m_vecExperiments.push_back(m_vecStrExperiments.at(i));
				mergedProteinTemp.m_vecMaxquantIBAQOfExperiments.push_back(m_vecProteins[j].m_dMaxquantIBAQ);
				mergedProteinTemp.m_vecMaxquantLFQOfExperiments.push_back(m_vecProteins[j].m_dLFQ);
				mergedProteinTemp.m_vecSPEAQInPeptidesOfExperiments.push_back(m_vecProteins[j].m_dProteinSPEAQInPeptides);
				mergedProteinTemp.m_vecSPEAQInProteinOfExperiments.push_back(m_vecProteins[j].m_dProteinSPEAQInProtein);
				MergedProteins.m_vecMergedProteins.push_back(mergedProteinTemp);
			}

		}
		else
		{
			for (int j = 0; j < m_vecProteins.size(); j++)
			{
				int index = 0;
				while ((index<MergedProteins.m_vecMergedProteins.size())&&(MergedProteins.m_vecMergedProteins[index].m_strPrtoteinName != m_vecProteins[j].m_strProteinID))
				{
						index++;
				}
				if (index >= MergedProteins.m_vecMergedProteins.size())
				{
					//cout << "Can not find " << m_vecProteins[j].m_strProteinID <<" in "<<m_vecStrExperiments.at(i-1)<< endl;
					log << "Can not find " << m_vecProteins[j].m_strProteinID << " in " << m_vecStrExperiments.at(i - 1) << endl;
					mergedProteinTemp.Clear();
					mergedProteinTemp.m_strPrtoteinName = m_vecProteins[j].m_strProteinID;
					for (int t = 0; t <= i;t++)
						mergedProteinTemp.m_vecExperiments.push_back(m_vecStrExperiments.at(t));
					for (int t = 0; t < i; t++)
						mergedProteinTemp.m_vecMaxquantIBAQOfExperiments.push_back(0.0);
					mergedProteinTemp.m_vecMaxquantIBAQOfExperiments.push_back(m_vecProteins[j].m_dMaxquantIBAQ);
					for (int t = 0; t < i; t++)
						mergedProteinTemp.m_vecMaxquantLFQOfExperiments.push_back(0.0);
					mergedProteinTemp.m_vecMaxquantLFQOfExperiments.push_back(m_vecProteins[j].m_dLFQ);
					for (int t = 0; t < i; t++)
						mergedProteinTemp.m_vecSPEAQInPeptidesOfExperiments.push_back(0.0);
					mergedProteinTemp.m_vecSPEAQInPeptidesOfExperiments.push_back(m_vecProteins[j].m_dProteinSPEAQInPeptides);
					for (int t = 0; t < i; t++)
						mergedProteinTemp.m_vecSPEAQInProteinOfExperiments.push_back(0.0);
					mergedProteinTemp.m_vecSPEAQInProteinOfExperiments.push_back(m_vecProteins[j].m_dProteinSPEAQInProtein);
					MergedProteins.m_vecMergedProteins.push_back(mergedProteinTemp);

					j++;
					continue;
				}

				MergedProteins.m_vecMergedProteins[index].m_vecExperiments.push_back(m_vecStrExperiments.at(i));
				MergedProteins.m_vecMergedProteins[index].m_vecMaxquantIBAQOfExperiments.push_back(m_vecProteins[j].m_dMaxquantIBAQ);
				MergedProteins.m_vecMergedProteins[index].m_vecMaxquantLFQOfExperiments.push_back(m_vecProteins[j].m_dLFQ);
				MergedProteins.m_vecMergedProteins[index].m_vecSPEAQInProteinOfExperiments.push_back(m_vecProteins[j].m_dProteinSPEAQInProtein);
				MergedProteins.m_vecMergedProteins[index].m_vecSPEAQInPeptidesOfExperiments.push_back(m_vecProteins[j].m_dProteinSPEAQInPeptides);
			}

		}
		

	}

	MergedProteins.mf_ExperimentsAnalysis();

	ofstream oMergefile;
	oMergefile.open(ProteinSPEAQInPeptidesResultPath.c_str());
	if (!oMergefile)
	{
		cout << "Can not open " << ProteinSPEAQInPeptidesResultPath << endl;
		log << "Can not open " << ProteinSPEAQInPeptidesResultPath << endl;
		exit(0);
	}
	else
	{
		cout << "Begin Merge " << m_vecProteins.size() << " Proteins' SPEAQInPeptides in " << ProteinSPEAQInPeptidesResultPath << endl;
		log << "Begin Merge " << m_vecProteins.size() << " Proteins' SPEAQInPeptides in " << ProteinSPEAQInPeptidesResultPath << endl;
	}




	vector<CProtein>::iterator proteinIter;
	vector<double>::iterator peptidesIter;
	oMergefile << "Protein ID\tExperiments\tMaxquantIBAQ\tLFQ\tSPEAQInPeptides\tSPEAQInProtein\t";
	oMergefile<<"MaxquantIBAQCV\tMaxquantLFQCV\tSPEAQInPeptidesCV\tSPEAQInProteinCV\n";
	int i = 0;
	int j = 0;
	for (i = 0; i < MergedProteins.m_vecMergedProteins.size(); i++)
	{
		oMergefile << MergedProteins.m_vecMergedProteins[i].m_strPrtoteinName << "\t";
		for (j = 0; j < MergedProteins.m_vecMergedProteins[i].m_vecExperiments.size(); j++)
		{
			oMergefile << MergedProteins.m_vecMergedProteins[i].m_vecExperiments.at(j) << ";";
		}
		oMergefile << "\t";
		for (j = 0; j < MergedProteins.m_vecMergedProteins[i].m_vecMaxquantIBAQOfExperiments.size(); j++)
		{
			oMergefile << MergedProteins.m_vecMergedProteins[i].m_vecMaxquantIBAQOfExperiments.at(j) << ";";
		}
		oMergefile << "\t";
		for (j = 0; j < MergedProteins.m_vecMergedProteins[i].m_vecMaxquantLFQOfExperiments.size(); j++)
		{
			oMergefile << MergedProteins.m_vecMergedProteins[i].m_vecMaxquantLFQOfExperiments.at(j) << ";";
		}
		oMergefile << "\t";
		for (j = 0; j < MergedProteins.m_vecMergedProteins[i].m_vecSPEAQInPeptidesOfExperiments.size(); j++)
		{
			oMergefile << MergedProteins.m_vecMergedProteins[i].m_vecSPEAQInPeptidesOfExperiments.at(j) << ";";
		}
		oMergefile << "\t";
		for (j = 0; j < MergedProteins.m_vecMergedProteins[i].m_vecSPEAQInProteinOfExperiments.size(); j++)
		{
			oMergefile << MergedProteins.m_vecMergedProteins[i].m_vecSPEAQInProteinOfExperiments.at(j) << ";";
		}
		oMergefile << "\t";

		oMergefile << MergedProteins.m_vecMaxquantIBAQCV.at(i) << "\t";
		oMergefile << MergedProteins.m_vecLFQCV.at(i) << "\t";
		oMergefile << MergedProteins.m_vecSPEAQInPeptidesCV.at(i) << "\t";
		oMergefile << MergedProteins.m_vecSPEAQInProteinCV.at(i) << "\n";

	}



}

void CProteinWorker::mf_MergeRegressionResult(std::ofstream& log, string RegressionResultPath)
{
	string strRegressionResultOfOneExperimentPath;
	FILE * Ofile; 
	char Buffer[BUFFERLENGTH];
	char *pstr;
	char *pstr1;
	string strIDTemp;
	string strTemp;


	ofstream oMergefile;
	oMergefile.open(RegressionResultPath.c_str());
	if (!oMergefile)
	{
		cout << "Can not open " << RegressionResultPath << endl;
		log << "Error:\tCan not open " << RegressionResultPath << endl;
		exit(0);
	}
	else
	{
		log << "Begin Merge regression result to " << RegressionResultPath << endl;
	}
	for (int i = 0; i < m_vecStrExperiments.size(); i++)
	{
		strRegressionResultOfOneExperimentPath = m_predictParam.m_strRegressionResult + m_vecStrExperiments.at(i) + ".txt";
		Ofile = fopen(strRegressionResultOfOneExperimentPath.c_str(), "r");
		if (!Ofile)
		{
			log << "Error:\tCan not open " << strRegressionResultOfOneExperimentPath << endl;
			exit(0);
		}

		fgets(Buffer, BUFFERLENGTH, Ofile);
		pstr = Buffer;
		oMergefile << m_vecStrExperiments.at(i) << endl;
		while (!feof(Ofile))
		{
			oMergefile << pstr;
			fgets(Buffer, BUFFERLENGTH, Ofile);
			pstr = Buffer;
		}

		fclose(Ofile);
	}

}
//to validation the correctness of the program, calculate the pearson correlation 
// between  the UPS2 protein's quantifications and UPS2 protein's real moles.
void CProteinWorker::mf_ShowUPS2AnalysisResult(std::ofstream& log, int iExperimentIndex)
{
	map<string, double> mapUPS2Id2Mols;
	map<string, double> mapUPS2Id2MW;
	map<string, double> mapUPS2iBAQs;
	map<string, double>::iterator Ups2MolsIter;
	map<string, double>::iterator Ups2WMIter;
	map<string, double>::iterator UPS2IBAQIter;
	vector<CProtein>::iterator proteinIter;
	string strTemp;

	string UPS2path = "UPS2_mols.ini";
	mf_LoadUPS2Mols(log,UPS2path, mapUPS2Id2Mols, mapUPS2Id2MW);

	double dPearsonCorr = 0.0;
	vector<string> vecUPS2IDIdentified;
	vector<double> vecUPS2logSPEAQInPeptides;
	vector<double> vecUPS2logSPEAQInProtein;
	vector<double> vecUPS2logMolsIdentified;
	vector<double> vecUPS2logIBAQIdentified;
	vector<double> vecUPS2logLFQIdentified;
	int iNumberofUPS2Identified;
	string::size_type stStart,stPosition;

	ofstream ofile;
	int iPosTemp = m_predictParam.m_strProteinSPEAQInPeptidesPath.find_last_of("\\");
	string CorrelationPath = m_predictParam.m_strProteinSPEAQInPeptidesPath.substr(0, iPosTemp + 1) + "CheckCorrelationBetweenPredicted_realUPS2Identified"+m_vecStrExperiments.at(iExperimentIndex)+".txt";

	ofile.open(CorrelationPath);
	ofile << "UPS2 Protein\tPeptidesNumber\tlogUPS2Mols\tLogProteinSPEAQInPeptides\tLogProteinSPEAQInProtein\tMaxquantIBAQ\tMaxquantLFQ\n";
	for (proteinIter = m_vecProteins.begin(); proteinIter != m_vecProteins.end(); proteinIter++)
	{
		if (proteinIter->m_bIfCalculateSPEAQInPeptides == true)
		{
			strTemp = proteinIter->m_strProteinID;
			stStart = strTemp.find("|");
			stPosition = strTemp.find("ups");
			if (stPosition != strTemp.npos)
			{
				if (stStart != strTemp.npos)
				{
					strTemp = strTemp.substr(stStart + 1, stPosition - stStart - 1);
				}
				else
				{
					strTemp = strTemp.substr(0, stPosition);
				}
				//if (strTemp == "P07339")
				//	strTemp = "O76070";
				//if (strTemp == "P08311")
				//	strTemp = "P01579";
				Ups2MolsIter = mapUPS2Id2Mols.find(strTemp);
				if (Ups2MolsIter == mapUPS2Id2Mols.end())
				{
					cout << "Can not find " << strTemp << " in UPS2 File\n";
					log << " Can not find " << strTemp << " in UPS2 File\n";
					continue;
				}
				Ups2WMIter = mapUPS2Id2MW.find(strTemp);
				if (Ups2WMIter == mapUPS2Id2MW.end())
				{
					cout << " Can not find " << strTemp << " in UPS2 File\n";
					log << " Can not find " << strTemp << " in UPS2 File\n";
					continue;
				}

				vecUPS2IDIdentified.push_back(strTemp);
				ofile << strTemp << "\t";
				ofile << proteinIter->m_iPeptidesNumber << "\t";
				ofile << log10(Ups2MolsIter->second) << "\t";
				if (iExperimentIndex == 0)
				{
					if (proteinIter->m_dProteinSPEAQInPeptides != 0.0 && (proteinIter->m_dMaxquantIBAQ != 0.0 || (!m_bIfIBAQIntensityExist)))
					{
						vecUPS2logMolsIdentified.push_back(log10(Ups2MolsIter->second));
						vecUPS2logSPEAQInPeptides.push_back(log10(proteinIter->m_dProteinSPEAQInPeptides));
						vecUPS2logSPEAQInProtein.push_back(log10(proteinIter->m_dProteinSPEAQInProtein));
						ofile << log10(proteinIter->m_dProteinSPEAQInPeptides) << "\t";
						ofile << log10(proteinIter->m_dProteinSPEAQInProtein) << "\t";

						if (m_bIfIBAQIntensityExist)
						{
							vecUPS2logIBAQIdentified.push_back(log10(proteinIter->m_dMaxquantIBAQ));
							ofile << log10(proteinIter->m_dMaxquantIBAQ) << "\t";
						}

						if (m_bIFLFQIntensityExist)
						{
							vecUPS2logLFQIdentified.push_back(log10(proteinIter->m_dLFQ));
							ofile << log10(proteinIter->m_dLFQ) << "\t";
						}
						ofile << endl;
					} // end if 0.0
					else
					{
						ofile << log10(proteinIter->m_dProteinSPEAQInPeptides) << "\t";
						ofile << log10(proteinIter->m_dProteinSPEAQInProtein) << "\t";
						ofile << log10(proteinIter->m_dMaxquantIBAQ) << "\t";
						ofile << log10(proteinIter->m_dLFQ) << endl;
					}
				}
				else
				{
					if (proteinIter->m_dProteinSPEAQInPeptides != 0.0 && (proteinIter->m_dMaxquantIBAQ != 0.0 || (!m_bIfIBAQIntensityExist)) && (proteinIter->m_dLFQ != 0.0 || (!m_bIFLFQIntensityExist)))
					{
						vecUPS2logMolsIdentified.push_back(log10(Ups2MolsIter->second));
						vecUPS2logSPEAQInPeptides.push_back(log10(proteinIter->m_dProteinSPEAQInPeptides));
						vecUPS2logSPEAQInProtein.push_back(log10(proteinIter->m_dProteinSPEAQInProtein));
						ofile << log10(proteinIter->m_dProteinSPEAQInPeptides) << "\t";
						ofile << log10(proteinIter->m_dProteinSPEAQInProtein) << "\t";

						if (m_bIfIBAQIntensityExist)
						{
							vecUPS2logIBAQIdentified.push_back(log10(proteinIter->m_dMaxquantIBAQ));
							ofile << log10(proteinIter->m_dMaxquantIBAQ) << "\t";
						}

						if (m_bIFLFQIntensityExist)
						{
							vecUPS2logLFQIdentified.push_back(log10(proteinIter->m_dLFQ));
							ofile << log10(proteinIter->m_dLFQ) << "\t";
						}
						ofile << endl;
					} // end if 0.0
					else
					{
						ofile << log10(proteinIter->m_dProteinSPEAQInPeptides) << "\t";
						ofile << log10(proteinIter->m_dProteinSPEAQInProtein) << "\t";
						ofile << log10(proteinIter->m_dMaxquantIBAQ) << "\t";
						ofile << log10(proteinIter->m_dLFQ) << endl;
					}
				}
			}

		}

	}
	ofile.close();

	dPearsonCorr = spearsonCorrelation(vecUPS2logSPEAQInPeptides,vecUPS2logMolsIdentified);
	cout << "The Correlation between " << vecUPS2IDIdentified.size() << "  UPS2's predicted SPEAQinPeptides and Real UPS2 Mols is " << dPearsonCorr << endl;
	log << "The Correlation between " << vecUPS2IDIdentified.size() << " UPS2's predicted SPEAQinPeptides and Real UPS2 Mols is " << dPearsonCorr << endl;

	dPearsonCorr = spearsonCorrelation(vecUPS2logSPEAQInProtein, vecUPS2logMolsIdentified);
	cout << "The Correlation between " << vecUPS2IDIdentified.size() << "  UPS2's predicted SPEAQInProtein and Real UPS2 Mols is " << dPearsonCorr << endl;
	log << "The Correlation between " << vecUPS2IDIdentified.size() << " UPS2's predicted SPEAQInProtein and Real UPS2 Mols is " << dPearsonCorr << endl;

	dPearsonCorr = spearsonCorrelation(vecUPS2logIBAQIdentified, vecUPS2logMolsIdentified);
	cout << "The Correlation between " << vecUPS2logIBAQIdentified.size() << "  UPS2 IBAQ and Real UPS2 Mols is " << dPearsonCorr << endl;
	log << "The Correlation between " << vecUPS2logIBAQIdentified.size() << "  UPS2 IBAQ and Real UPS2 Mols is " << dPearsonCorr << endl;

	dPearsonCorr = spearsonCorrelation(vecUPS2logLFQIdentified, vecUPS2logMolsIdentified);
	cout << "The Correlation between " << vecUPS2logLFQIdentified.size() << "  UPS2 LFQ and Real UPS2 Mols is " << dPearsonCorr << endl;
	log << "The Correlation between " << vecUPS2logLFQIdentified.size() << "  UPS2 LFQ and Real UPS2 Mols is " << dPearsonCorr << endl;

	//calculate the CV of Peptides which belong to the proteins whose moles is same
	string PeptidesCVForSameUPS2MolsPath = m_predictParam.m_strProteinSPEAQInPeptidesPath.substr(0, iPosTemp + 1) + "PeptidesCVForSameUPS2MolsPath"+m_vecStrExperiments.at(iExperimentIndex)+".txt";
	mf_AnalysisPeptidesCVForSameUPS2Mols(log, PeptidesCVForSameUPS2MolsPath, m_vecProteins, mapUPS2Id2Mols);

}

int CProteinWorker::mf_GetExperimentNames(std::ofstream& log, string experimentalDesignPath)
{
	ifstream fin(experimentalDesignPath.c_str());
	if (!fin)
	{
		cout << "Can not open " << experimentalDesignPath << endl;
		log << "Error:\tCan not open " << experimentalDesignPath << endl;
		exit(0);
	}

	map<string, int> mapExperiments;
	map<string, int>::iterator mapExperimentsIter;
	string strTemp1, strTemp2;
	int iTemp1 = 0, iTemp2 = 0;
	m_iNumberOfExperiments = 0;
	getline(fin, strTemp1);//jump the first row
	strTemp2 = "";
	mapExperiments.insert(pair<string, int>(strTemp2, 0));
	while (getline(fin, strTemp1))
	{
		iTemp1 = strTemp1.find("\t");
		iTemp1 = strTemp1.find("\t", iTemp1 + 1);
		iTemp2 = strTemp1.find("\n", iTemp1 + 1);
		strTemp2 = strTemp1.substr(iTemp1 + 1, iTemp2 - iTemp1 - 1);
		mapExperiments.insert(pair<string, int>(strTemp2, 0));
	}

	m_iNumberOfExperiments = mapExperiments.size();
	mapExperimentsIter = mapExperiments.begin();
	for (; mapExperimentsIter != mapExperiments.end(); mapExperimentsIter++)
	{
		m_vecStrExperiments.push_back(mapExperimentsIter->first);
	}
	return m_iNumberOfExperiments;
}

bool CProteinWorker::mf_LoadUPS2Mols(std::ofstream& log, string path, map<string, double>&mapUPS2Id2Mols, map<string, double>&mapUPS2Id2MW)
{

	mapUPS2Id2MW.clear();
	mapUPS2Id2Mols.clear();

	char Buffer[BUFFERLENGTH];
	char *pstr;
	char *pstr1;
	string strTemp;
	FILE * pFile;

	cout << "\tLoading UPS2 mols\n";
	log << "\tLoading UPS2 mols\n";

	double dMolsTemp;
	double dWMTemp;


	pFile = fopen(path.c_str(), "r");
	if (pFile == NULL)
	{
		cout << "Can not open " << path<< endl;
		log << "Error:\tCan not open " << path << endl;
		exit(0);
	}

	fgets(Buffer, BUFFERLENGTH, pFile);
	fgets(Buffer, BUFFERLENGTH, pFile);

	while (!feof(pFile))
	{
		pstr = Buffer;
		pstr1 = strstr(pstr, "\t");
		if (pstr1 != NULL)
		{
			*pstr1 = '\0';
		}
		else
		{
			log << "Error:\tThe format of UPS2_mols.ini is wrong.\n";
			log << "Please do not modify the configuration file arbitrarily!\n";
			exit(0);
		}
		strTemp = pstr;

		pstr = pstr1 + 1;
		pstr1 = strstr(pstr, "\t");

		pstr = pstr1 + 1;
		pstr1 = strstr(pstr, "\t");
		if (pstr1 != NULL)
		{
			*pstr1 = '\0';
		}
		else
		{
			log << "Error:\tThe format of UPS2_mols.ini is wrong.\n";
			log << "Please do not modify the configuration file arbitrarily!\n";
			exit(0);
		}
		dMolsTemp = atof(pstr);
		mapUPS2Id2Mols.insert(pair<string, double>(strTemp, dMolsTemp));

		pstr = pstr1 + 1;
		pstr1 = strstr(pstr, "\t");

		pstr = pstr1 + 1;
		pstr1 = strstr(pstr, "\t");
		if (pstr1 != NULL)
		{
			*pstr1 = '\0';
		}
		else
		{
			log << "Error:\tThe format of UPS2_mols.ini is wrong.\n";
			log << "Please do not modify the configuration file arbitrarily!\n";
			exit(0);
		}
		dWMTemp = atof(pstr);
		mapUPS2Id2MW.insert(pair<string, double>(strTemp, dWMTemp));

		fgets(Buffer, BUFFERLENGTH, pFile);
	}
	fclose(pFile);
	return true;
}

bool CProteinWorker::mf_AnalysisPeptidesCVForSameUPS2Mols(std::ofstream& log, string path, vector<CProtein>& vecProteins, map<string, double>&mapUPS2Id2Mols)
{
	ofstream ofileCVForSameMols(path);
	if (!ofileCVForSameMols)
	{
		cout << "can not open " <<path<< endl;
		log << "can not open " << path << endl;
	}

	map<double, vector<double>> mapdUPS2Mols2PeptidesNativeIntensitys;
	map<double, vector<double>> mapdUPS2Mols2PeptidesAdjustIntensitys;
	map<double, vector<double>>::iterator mapdUPS2Mols2PeptidesNativeIntensitysIter;
	map<double, vector<double>>::iterator mapdUPS2Mols2PeptidesAdjustIntensitysIter;
	vector<double> vecdPeptideIntensityTemp;
	vector<CProtein>::iterator proteinIter;
	string::size_type stStart, stPosition;
	string strTemp;
	map<string, double>::iterator Ups2MolsIter;

	int i = 0; // just for index
	for (proteinIter = vecProteins.begin(); proteinIter != vecProteins.end(); proteinIter++)
	{
		if (proteinIter->m_bIfCalculateSPEAQInPeptides == true)
		{
			strTemp = proteinIter->m_strProteinID;
			stStart = strTemp.find("|");
			stPosition = strTemp.find("ups");
			if (stPosition != strTemp.npos)
			{
				if (stStart != strTemp.npos)
				{
					strTemp = strTemp.substr(stStart + 1, stPosition - stStart - 1);
				}
				else
				{
					strTemp = strTemp.substr(0, stPosition);
				}
				Ups2MolsIter = mapUPS2Id2Mols.find(strTemp);
				if (Ups2MolsIter == mapUPS2Id2Mols.end())
				{
					cout << "Error! Can not find " << strTemp << " in UPS2 File\n";
					log << "Error! Can not find " << strTemp << " in UPS2 File\n";
					continue;
				}
				else
				{
					mapdUPS2Mols2PeptidesNativeIntensitysIter = mapdUPS2Mols2PeptidesNativeIntensitys.find(Ups2MolsIter->second);
					mapdUPS2Mols2PeptidesAdjustIntensitysIter = mapdUPS2Mols2PeptidesAdjustIntensitys.find(Ups2MolsIter->second);
					if (mapdUPS2Mols2PeptidesNativeIntensitysIter == mapdUPS2Mols2PeptidesNativeIntensitys.end())
					{//new mols number group;
						vecdPeptideIntensityTemp.clear();
						vecdPeptideIntensityTemp = proteinIter->m_vPeptidesNativeIntensity;
						mapdUPS2Mols2PeptidesNativeIntensitys.insert(pair<double, vector<double>>(Ups2MolsIter->second, vecdPeptideIntensityTemp));

						vecdPeptideIntensityTemp.clear();
						vecdPeptideIntensityTemp = proteinIter->m_vPeptidesAdjustIntensity;
						mapdUPS2Mols2PeptidesAdjustIntensitys.insert(pair<double, vector<double>>(Ups2MolsIter->second, vecdPeptideIntensityTemp));

					}
					else
					{
						for (i = 0; i < proteinIter->m_vPeptidesNativeIntensity.size(); i++)
						{
							mapdUPS2Mols2PeptidesNativeIntensitysIter->second.push_back(proteinIter->m_vPeptidesNativeIntensity.at(i));
						}
						for (i = 0; i < proteinIter->m_vPeptidesAdjustIntensity.size(); i++)
						{
							mapdUPS2Mols2PeptidesAdjustIntensitysIter->second.push_back(proteinIter->m_vPeptidesAdjustIntensity.at(i));
						}
					}//end if new mols

				}
			}//end if UPS2 exist
		}// end if (proteinIter->m_bIfCalculateSPEAQInPeptides == true)
	} //end for proteins

	mapdUPS2Mols2PeptidesNativeIntensitysIter = mapdUPS2Mols2PeptidesNativeIntensitys.begin();
	mapdUPS2Mols2PeptidesAdjustIntensitysIter = mapdUPS2Mols2PeptidesAdjustIntensitys.begin();

	for (; mapdUPS2Mols2PeptidesNativeIntensitysIter != mapdUPS2Mols2PeptidesNativeIntensitys.end(); mapdUPS2Mols2PeptidesNativeIntensitysIter++, mapdUPS2Mols2PeptidesAdjustIntensitysIter++)
	{
		ofileCVForSameMols << mapdUPS2Mols2PeptidesNativeIntensitysIter->first << "\t";
		ofileCVForSameMols << CalculateCV(mapdUPS2Mols2PeptidesNativeIntensitysIter->second) << "\t";
		ofileCVForSameMols << CalculateCV(mapdUPS2Mols2PeptidesAdjustIntensitysIter->second) << endl;
	}

	ofileCVForSameMols.close();
	return true;
}