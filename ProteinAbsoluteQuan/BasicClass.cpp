#include"stdafx.h"
#include"BasicClass.h"

CProtein::CProtein()
{
	m_dPeptidesIntensityMax = 0.0;
	m_dPeptidesIntensityMedian = 0.0;
	m_dProteinSPEAQInPeptides = 0.0;
	m_dProteinSPEAQInProtein = 0.0;
	m_dProteinCorrectionFactor = 0.0;
	m_dMaxquantIBAQ = 0.0;
	m_dLFQ = 0.0;
	m_iPeptidesNumber = 0;
	m_dPeptidesIntensitySum = 0.0;
	m_dCVNative = 0.0;
	m_dCVAfterAdjustAfter = 0.0;
	m_bIfCalculateSPEAQInPeptides = false;
}
void CProtein::Clear()
{
	m_vPeptidesIDs.clear();
	m_vPeptidesSequencesWithFlankingRegion.clear();
	m_vPeptidesSequences.clear();
	m_vPeptidesAttributes.clear();
	m_vAdjustPeptidesIntensityFactor.clear();
	m_vPeptidesNativeIntensity.clear();
	m_vPeptidesAdjustIntensity.clear();
	m_dPeptidesIntensityMedian = 0.0;
	m_dPeptidesIntensityMax = 0.0;
	m_mapPeptidesSequenceAndCorrespondRawFiles.clear();
	m_vPeptidesCorrectionFactor.clear();
	m_dMaxquantIBAQ = 0.0;
	m_vdPeptidesMW.clear();
	m_dReCalculateIBAQ = 0.0;
	m_dLFQ = 0.0;
}

string CProtein::mf_GetPeptidesAdjacentSequence(std::ofstream& log, string PeptideSequence)
{
	string::size_type peptideLocation = 0, begin = 0, Getwidth = 0, end;
	end = PeptideSequence.find("{");
	PeptideSequence = PeptideSequence.substr(0, end);
	peptideLocation = m_strProteinSequence.find(PeptideSequence);
	if (peptideLocation != m_strProteinSequence.npos)
	{

		if ((peptideLocation >= FlankingRegionWidth) && (peptideLocation <= m_strProteinSequence.size() - PeptideSequence.size() - FlankingRegionWidth))
		{
			begin = peptideLocation - FlankingRegionWidth;
			Getwidth = PeptideSequence.size() + 2 * FlankingRegionWidth;
		}
		else if (peptideLocation < FlankingRegionWidth)
		{
			begin = 0;
			Getwidth = peptideLocation + PeptideSequence.size() + FlankingRegionWidth;
		}
		else if (peptideLocation > m_strProteinSequence.size() - PeptideSequence.size() - FlankingRegionWidth)
		{
			begin = peptideLocation - FlankingRegionWidth;
			Getwidth = m_strProteinSequence.size() - peptideLocation - 1 + PeptideSequence.size() + FlankingRegionWidth;
		}
		return m_strProteinSequence.substr(begin, Getwidth);
	}
	else
	{
		cout << "Can not find peptide " << PeptideSequence << " in protein " << m_strProteinID << endl;
		cout << "The protein sequence is " << m_strProteinSequence << endl; 
		log << "Can not find peptide " << PeptideSequence << " in protein " << m_strProteinID << endl;
		log << "The protein sequence is " << m_strProteinSequence << endl;
		return "NULL";
	}

}
void Cpeptides::mf_Clear()
{
	m_mapPeptideIDAndSequence.clear();
	m_mapPeptideIDAndIntensity.clear();
	m_mapPeptideSequenceAndIntensity.clear();
	m_mapPeptideSequencuAndIfUse.clear();
	m_mapPeptideSequenceAndAttribute.clear();
}

void CMergedProtein::Clear()
{
	m_vecExperiments.clear();
	m_vecMaxquantIBAQOfExperiments.clear();
	m_vecMaxquantLFQOfExperiments.clear();
	m_vecSPEAQInProteinOfExperiments.clear();
	m_vecSPEAQInPeptidesOfExperiments.clear();
}

void CMergedProteins::mf_ExperimentsAnalysis()
{
	double dCVTemp;
	for (int i = 0; i < m_vecMergedProteins.size(); i++)
	{
		dCVTemp = CalculateCV(m_vecMergedProteins[i].m_vecMaxquantIBAQOfExperiments);
		m_vecMaxquantIBAQCV.push_back(dCVTemp);
		dCVTemp = CalculateCV(m_vecMergedProteins[i].m_vecMaxquantLFQOfExperiments);
		m_vecLFQCV.push_back(dCVTemp);
		dCVTemp = CalculateCV(m_vecMergedProteins[i].m_vecSPEAQInPeptidesOfExperiments);
		m_vecSPEAQInPeptidesCV.push_back(dCVTemp);
		dCVTemp = CalculateCV(m_vecMergedProteins[i].m_vecSPEAQInProteinOfExperiments);
		m_vecSPEAQInProteinCV.push_back(dCVTemp);
	}
}

double Average(double *Array, int Number)
{
	double ave = 0.0;
	for (int i = 0; i < Number; i++)
	{
		ave += Array[i];
	}
	return ave / Number;

}
double Average(vector<double> Array)
{
	double ave = 0.0;
	vector<double>::iterator doubleIter;
	for (doubleIter = Array.begin(); doubleIter != Array.end(); doubleIter++)
	{
		ave += *doubleIter;
	}
	double len = Array.size();
	return ave / len;

}

void quickSort(double v[], int indexTemp[], int left, int right)
{
	if (left < right)
	{
		double key = v[left];
		int indexKey = indexTemp[left];
		int low = left;
		int high = right;
		while (low < high)
		{
			while (low < high && v[high] <= key)
			{
				high--;
			}
			if (low < high)
			{
				v[low] = v[high];
				indexTemp[low] = indexTemp[high];
				low++;
			}

			while (low < high && v[low] > key)
			{
				low++;
			}
			if (low < high)
			{
				v[high] = v[low];
				indexTemp[high] = indexTemp[low];
				high--;
			}
		}
		v[low] = key;
		indexTemp[low] = indexKey;
		quickSort(v, indexTemp, left, low - 1);
		quickSort(v, indexTemp, low + 1, right);
	}
}

void quickSort(vector<size_t> &v, vector<size_t> &indexTemp, int left, int right)
{
	if (left < right)
	{
		size_t key = v[left];
		size_t indexKey = indexTemp[left];
		int low = left;
		int high = right;
		while (low < high)
		{
			while (low < high && v[high] <= key)
			{
				high--;
			}
			if (low < high)
			{
				v[low] = v[high];
				indexTemp[low] = indexTemp[high];
				low++;
			}

			while (low < high && v[low] > key)
			{
				low++;
			}
			if (low < high)
			{
				v[high] = v[low];
				indexTemp[high] = indexTemp[low];
				high--;
			}
		}
		v[low] = key;
		indexTemp[low] = indexKey;
		quickSort(v, indexTemp, left, low - 1);
		quickSort(v, indexTemp, low + 1, right);
	}
}
// calculat Gene abundance from proteins abundance
void CGene::mf_CaluculateIntensityFromProteinsIntensity1()
{
	m_dExpressionIntensity = 0.0;
	int DigestionNumber = 0;
	map<string, vector<double>>::iterator PRoteinIDAndPeptidesIntensitysIter;
	int i = 0;
	for (PRoteinIDAndPeptidesIntensitysIter = m_mapPRoteinIDAndPeptidesIntensitys.begin(); PRoteinIDAndPeptidesIntensitysIter != m_mapPRoteinIDAndPeptidesIntensitys.end(); PRoteinIDAndPeptidesIntensitysIter++)
	{
		for (i = 0; i < PRoteinIDAndPeptidesIntensitysIter->second.size(); i++)
			m_dExpressionIntensity += PRoteinIDAndPeptidesIntensitysIter->second.at(i);
	}
	// get the number of peptides of theory enzyme of proteins
	map<string, int>::iterator ProteinIdAndDigestNumberIter;
	for (ProteinIdAndDigestNumberIter = m_mapProteinIdAndDigestNumber.begin(); ProteinIdAndDigestNumberIter != m_mapProteinIdAndDigestNumber.end(); ProteinIdAndDigestNumberIter++)
	{
		DigestionNumber += ProteinIdAndDigestNumberIter->second;
	}

	m_dExpressionIntensity = m_dExpressionIntensity / DigestionNumber;
}

void CGene::mf_CaluculateIntensityFromProteinsIntensity2()
{
	m_dExpressionIntensity = 0.0;
	int DigestionNumber = 0;
	map<string, double>::iterator ProteinIdAndIntensityIter;
	int i = 0;
	for (ProteinIdAndIntensityIter = m_mapProteinIDAndIntensity.begin(); ProteinIdAndIntensityIter != m_mapProteinIDAndIntensity.end(); ProteinIdAndIntensityIter++)
	{
		m_dExpressionIntensity += ProteinIdAndIntensityIter->second;
	}
}

double spearsonCorrelation(double s1[], double s2[], int NumberOfSamples)
{   //change m_iTestSampleNumber to NumberOfSamples by gzhq 20150529
	double dCorrelation = 0.0;
	double means2 = 0.0;
	double means1 = 0.0;
	double sd1 = 0.0, sd2 = 0.0;
	means1 = Average(s1, NumberOfSamples);
	for (int i = 0; i < NumberOfSamples; i++)
		means2 += s2[i];
	means2 = means2 / NumberOfSamples;
	for (int i = 0; i < NumberOfSamples; i++)
		dCorrelation += (s2[i] - means2)*(s1[i] - means1);
	for (int i = 0; i < NumberOfSamples; i++)
	{
		sd1 += (s1[i] - means1)*(s1[i] - means1);
		sd2 += (s2[i] - means2)*(s2[i] - means2);
	}
	dCorrelation = dCorrelation / sqrt(sd1*sd2);
	if (sd1*sd2 <= 0) //addded by CC 20150527
		return -1;
	else	return dCorrelation;
}

double spearsonCorrelation(vector<double> s1, double s2[])
{
	double dCorrelation = 0.0;
	double means2 = 0.0;
	double means1 = 0.0;
	double sd1 = 0.0, sd2 = 0.0;
	int NumberOfSamples = s1.size();
	means1 = Average(s1);
	means2 = Average(s2, NumberOfSamples);
	for (int i = 0; i < NumberOfSamples; i++)
		dCorrelation += (s2[i] - means2)*(s1[i] - means1);
	for (int i = 0; i < NumberOfSamples; i++)
	{
		sd1 += (s1[i] - means1)*(s1[i] - means1);
		sd2 += (s2[i] - means2)*(s2[i] - means2);
	}
	dCorrelation = dCorrelation / sqrt(sd1*sd2);
	if (sd1*sd2 <= 0) //addded by CC 20150527
		return -1;
	else	return dCorrelation;
}

double spearsonCorrelation(vector<double> s1, vector<double> s2)
{
	//delete the nan in s1 and s2;
	/*vector<double> s1, s2;
	for (int i = 0; i < os1.size(); i++)
	{
		if (_finite(os1[i]) &&( _finite(os2[i])))
		{
			s1.push_back(os1[i]);
			s2.push_back(os2[i]);
		}
	}*/

	double dCorrelation = 0.0;
	double means2 = 0.0;
	double means1 = 0.0;
	double sd1 = 0.0, sd2 = 0.0;
	int NumberOfSamples = s1.size();
	means1 = Average(s1);
	means2 = Average(s2);
	for (int i = 0; i < NumberOfSamples; i++)
		dCorrelation += (s2[i] - means2)*(s1[i] - means1);
	for (int i = 0; i < NumberOfSamples; i++)
	{
		sd1 += (s1[i] - means1)*(s1[i] - means1);
		sd2 += (s2[i] - means2)*(s2[i] - means2);
	}
	dCorrelation = dCorrelation / sqrt(sd1*sd2);
	if (sd1*sd2 <= 0) //addded by CC 20150527
		return -1;
	else	return dCorrelation;
}

double CalculateCV(vector<double>s)
{
	double mean = 0.0;
	double sd = 0.0, sd2 = 0.0;
	int len = s.size();


	for (int i = 0; i < len; i++)
	{
		//s[i] = log10(s[i]); //change intensity to log10(intensity) by CC 20150527
		mean += s[i];
	}

	if (len <= 1)  //addded by CC 20150527
		return 0.0;
	else 
	{
		mean = mean / len;
		for (int i = 0; i < len; i++)
		{
			sd += (s[i] - mean)*(s[i] - mean);
		}
	}

	sd = sqrt(sd / len);  //change by gzhq 20151118

	if (mean <= 0)  //addded by CC 20150527
		return -1;
	else	return sd / mean;
}

void Discretization(double *TestY, int m_iTestSampleNumber)
{

	vector <pair<double, int>> vctemp;
	int bin_num = m_iTestSampleNumber / 10;
	for (int i = 0; i < m_iTestSampleNumber; i++)
	{
		vctemp.push_back(make_pair(TestY[i], i));
	}
	sort(vctemp.begin(), vctemp.end());
	for (int i = 0; i < m_iTestSampleNumber; i++)
	{
		TestY[vctemp[i].second] = int(i / bin_num + 1)*0.1;
		if (TestY[vctemp[i].second]>1) //added by CC 20150527
			TestY[vctemp[i].second] = 1;
	}

}

void Discretization(vector<double> &TestY)
{

	vector <pair<double, int>> vctemp;
	int m_iTestSampleNumber = TestY.size();
	int bin_num = m_iTestSampleNumber / 10;
	for (int i = 0; i < m_iTestSampleNumber; i++)
	{
		vctemp.push_back(make_pair(TestY[i], i));
	}
	sort(vctemp.begin(), vctemp.end());
	for (int i = 0; i < m_iTestSampleNumber; i++)
	{
		TestY[vctemp[i].second] = int(i / bin_num + 1)*0.1;
		if (TestY[vctemp[i].second]>1) //added by CC 20150527
			TestY[vctemp[i].second] = 1;
	}

}
double GetMaxIntensity(const vector<double> PeptidesIntensity)
{
	double dMaxIntensitytemp = 0.0;
	for (int i = 0; i < PeptidesIntensity.size(); i++)
	{
		if (dMaxIntensitytemp < PeptidesIntensity.at(i))
		{
			dMaxIntensitytemp = PeptidesIntensity.at(i);
		}

	}
	return dMaxIntensitytemp;
}

double GetMedianIntensity(vector<double> PeptidesIntensity)
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
void Assertion(vector<double> &TestY)
{
	for (size_t i = 0; i < TestY.size(); i++)
	{
		if (TestY.at(i) <= 0.0)
			TestY.at(i) = 0.1;
	}
}

void Assertion(double *TestY, int m_iTestSampleNumber)
{

	for (int i = 0; i < m_iTestSampleNumber; i++)
	{
		if (TestY[i] <= 0.0)
		{
			TestY[i] = 0.1;
		}
	}

}
//integral using trapzoid formular
double IntegralFromVector(std::ofstream& log, vector<double> Values, vector<double> times)
{
	double result = 0.0;
	if (Values.size() != times.size())
	{
		cout << "The size of Values does not correspond to times'" << endl;
		log << "The size of Values does not correspond to times'" << endl;
		exit(0);
	}
	for (int i = 1; i < times.size(); i++)
	{
		result = result + Values.at(i)*(times.at(i) - times.at(i - 1)) / 2;
	}
	return result;
}

void GetQuantiles(vector<double> PeptidesIntensity, double &FirstQ, double &median, double &ThirdQ)
{
	sort(PeptidesIntensity.begin(), PeptidesIntensity.end());
	int len = PeptidesIntensity.size();
	if (len ==0)             //added by gzhq 20150529
	{
		FirstQ = 0.0;
		median = 0.0;
		ThirdQ = 0.0;
	}
	else if (len == 1)
	{
		FirstQ = PeptidesIntensity.front();
		median = PeptidesIntensity.front();
		ThirdQ = PeptidesIntensity.front();
	}
	else if (len == 2)
	{
		FirstQ = PeptidesIntensity.at(0);
		median = (PeptidesIntensity.at(0) + PeptidesIntensity.at(1)) / 2;
		ThirdQ = PeptidesIntensity.at(1);
	}
	else if (len == 3)
	{
		FirstQ = PeptidesIntensity.at(0);
		median = PeptidesIntensity.at(1);
		ThirdQ = PeptidesIntensity.at(2);
	}
	else if (len % 4 == 0)
	{
		FirstQ = (PeptidesIntensity.at(len/4)+PeptidesIntensity.at(len/4-1))/2;
		median = (PeptidesIntensity.at(len / 2) + PeptidesIntensity.at(len /2-1))/ 2;
		ThirdQ = (PeptidesIntensity.at(3 * (len / 4)) + PeptidesIntensity.at(3 * (len / 4) - 1)) / 2;
	} 
	else if (len % 4 == 1)
	{
		FirstQ = (PeptidesIntensity.at(floor(len / 4)) + PeptidesIntensity.at(ceil(len / 4))) / 2;
		median = PeptidesIntensity.at(floor(len / 2));
		ThirdQ = (PeptidesIntensity.at(floor(3 * (len / 4))) + PeptidesIntensity.at(ceil(3 * (len / 4) ))) / 2;
	}
	else if (len % 4 == 2)
	{
		FirstQ = PeptidesIntensity.at(floor(len / 4)) ;
		median = (PeptidesIntensity.at(len / 2) + PeptidesIntensity.at(len / 2 - 1)) / 2;
		ThirdQ = PeptidesIntensity.at(floor(3 * (len / 4))) ;

	}
	else if (len % 4 == 3)
	{
		FirstQ = PeptidesIntensity.at(floor(len / 4));
		median = PeptidesIntensity.at(ceil(len / 2));
		ThirdQ = PeptidesIntensity.at(floor(3 * (len / 4)));
	}
}

void GetAttributesFromFirstRow(char * pstr, map<string, int>& mapAttributesAndColumns)
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
			cout << "Error:\tThe format  is wrong.\n";
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
		cout << "Error:\tThe format is wrong.\n";
		exit(0);
	}
	mapAttributesAndColumns.insert(pair<string, int>(pstr, icolumns));
}

bool fStringToBool(string str)
{
	transform(str.begin(), str.end(), str.begin(), ::tolower);
	if (str == "true")
		return true;
	else if (str == "false")
		return false;
	else
	{
		cout << "Can not convert " << str << endl;
	}
}
