#include"stdafx.h"
#include"BasicClass.h"
CProtein::CProtein()
{
	//m_dPeptidesIntensityMedian = 0.0;
	m_dProteinSPEAQInPeptides = 0.0;
	//m_dMaxquantIBAQ = 0.0;
	//m_dMaxquantLFQ = 0.0;
	m_iPeptidesNumber = 0;
	m_dPeptidesIntensitySum = 0.0;
}
void CProtein::Clear()
{
	m_vPeptidesIDs.clear();
	//m_vPeptidesSequencesWithFlankingRegion.clear();
	m_vPeptidesSequences.clear();
	m_mapExperimentAndPeptidesIntensity.clear();
	m_vecdMaxquantIBAQ.clear();
	m_vecdMaxquantLFQ.clear();
	//m_vPeptidesAttributes.clear();
	//m_vPeptidesIntensity.clear();
//	m_vPeptidesNativeIntensity.clear();
	//m_dPeptidesIntensityMedian = 0.0;
	m_mapPeptidesSequenceAndCorrespondRawFiles.clear();
}

string CProtein::mf_GetPeptidesAdjacentSequence(std::ofstream& log, string PeptideSequence)
{
	string::size_type peptideLocation = 0, begin = 0, Getwidth = 0,end;
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
	m_mapPeptideIDAndIntensitys.clear();
	m_mapPeptideSequenceAndIntensitys.clear();
	m_mapPeptideSequencuAndIfUse.clear();
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


double spearsonCorrelation(double s1[], double s2[], int NumberOfSamples)
{   //change m_NumberOfSamples to NumberOfSamples by gzhq 20150529
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

double CalculateCV(vector<double>s)
{
	double mean = 0.0;
	double sd = 0.0, sd2 = 0.0;
	int len = s.size();
	for (int i = 0; i < len; i++)
	{
		s[i] = log10(s[i]); //change intensity to log10(intensity) by CC 20150527
		mean += s[i];
	}
	mean = mean / len;
	for (int i = 0; i < len; i++)
	{
		sd += (s[i] - mean)*(s[i] - mean);
	}
	sd = sqrt(sd / len);
	if (mean <= 0)  //addded by CC 20150527
		return -1;
	else	return sd / mean;
}

double IntegralFromVector(std::ofstream& log, vector<double> Values, vector<double> times)
{
	double result = 0.0;
	if (Values.size() != times.size())
	{
		cout << "The size of Values does not correspond to times'" << endl;
		log << "Error:\tThe size of Values does not correspond to times'" << endl;
		exit(0);
	}
	for (int i = 1; i < times.size(); i++)
	{
		result = result + Values.at(i)*(times.at(i) - times.at(i - 1)) / 2;
	}
	return result;
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
//remove the whitespace at the end of the string or at the begin of the string
string strim(string &str) 
{
	int ilocation = str.find_last_of(' ');
	if (ilocation != -1)
	{
		str.erase(str.find_last_not_of(' ') + 1);    //remove the whitespace at the end of the string
	}
	ilocation = str.find_first_not_of(' ');
	str.erase(0, ilocation);    //remove the whitespace  at the begin of the string
	return str;
}