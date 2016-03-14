// ProteinAbsoluteQuan.cpp : the main file of the program

#include "stdafx.h"
#include"ProteinWorker.h"
#include<windows.h>
#include<time.h>

int _tmain(int argc, _TCHAR* argv[]) 
{
	clock_t begin, end;
	begin = clock();

	if (argc <= 1)
	{
		cout << "ERROR:Input parameters is not enough" << endl;
		cout << "Input parameters\'s Format:" << endl;
		cout << "xxx.txt" << endl;
		system("pause");
		return 1;
	}

	string ParaFilePath;
	DWORD dwNum = WideCharToMultiByte(CP_OEMCP, NULL, argv[1], -1, NULL, 0, NULL, FALSE);
	char *psText;
	psText = new char[dwNum];
	if (!psText)
	{
		delete[]psText;
	}
	WideCharToMultiByte(CP_OEMCP, NULL, argv[1], -1, psText, dwNum, NULL, FALSE);
	ParaFilePath = psText;
	delete[]psText;

	CProteinWorker ProteinWorker;
	string NativePeptidesPath;
	string strProteinSPEAQInPeptidesPath;
	
	ofstream flog("log.txt",ios::app);
	
	//set parameters
	ProteinWorker.m_trainParam.mf_setparameters(flog, ParaFilePath);
	ProteinWorker.m_predictParam.mf_setparameters(flog,ParaFilePath);

	//load data
	ProteinWorker.mf_LoadProteinFasta(flog,ProteinWorker.m_predictParam.m_strProteinFastaPath, ProteinWorker.m_predictParam.m_eFastaType);
	int iNumberOfExperiments=ProteinWorker.mf_GetExperimentNames(flog,ProteinWorker.m_trainParam.m_strExperimentDesignPath);
	for (int i = 0; i < iNumberOfExperiments; i++)
	{
		flog << "In experiment " << ProteinWorker.m_vecStrExperiments.at(i) << endl;
		ProteinWorker.mf_LoadProteins(flog,ProteinWorker.m_trainParam,i);

		
		//quantify the proteins  
		ProteinWorker.mf_CalculateProteinIntensity(flog,i);

		// save results
		NativePeptidesPath = ProteinWorker.m_predictParam.m_strQuantiResultPath + "\\PeptidesAttribute"+ProteinWorker.m_vecStrExperiments.at(i)+".txt";
		ProteinWorker.mf_SaveProteinsNativePeptides(flog, NativePeptidesPath);
		strProteinSPEAQInPeptidesPath = ProteinWorker.m_predictParam.m_strProteinSPEAQInPeptidesPath+ProteinWorker.m_vecStrExperiments.at(i)+".txt";
		ProteinWorker.mf_saveProteinSPEAQInPeptides(flog, strProteinSPEAQInPeptidesPath);
		if (ProteinWorker.m_trainParam.m_bIfCotainStardProtein)
		{
			ProteinWorker.mf_ShowUPS2AnalysisResult(flog, i);
		}
		
	}
	string AllProteinSPEAQInPeptidesPath = ProteinWorker.m_predictParam.m_strQuantiResultPath + "\\MergedProteinSPEAQs.txt";
	ProteinWorker.mf_MergeProteinSPEAQInPeptides(flog, AllProteinSPEAQInPeptidesPath);

	string strAllRegressionResultPath = ProteinWorker.m_predictParam.m_strQuantiResultPath + "\\MergedRegressionResult.txt";
	ProteinWorker.mf_MergeRegressionResult(flog, strAllRegressionResultPath);


	end = clock();
	cout << "The Program run " << (end - begin) / CLOCKS_PER_SEC << " seconds." << endl;
	flog << "The Program run " << (end - begin) / CLOCKS_PER_SEC << " seconds." << endl;
	flog.close();
	return 0;
}

