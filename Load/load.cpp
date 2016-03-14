// load.cpp : Define the entry point of the console application
#include "stdafx.h"
#include"loadParam.h"
#include "BasicClass.h"
#include"DataIO.h"
#include"time.h"
#include<windows.h>
#include<iostream>
using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{
	clock_t begin, end;
	begin = clock();

	ofstream flog("log.txt",ios::out);
	cout << "Load\t" << endl;
	flog <<  "Load\t" << endl;
	CLoadParam param;
	if (argc <= 1)
	{
		cout << "ERROR:Input parameters is not enough" << endl;
		cout << "Input parameters\'s Format:" << endl;
		cout << "xxx.config" << endl;
		system("pause");
		return 0;
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
	psText = NULL;

	param.mf_setparameters(flog, ParaFilePath);
	if (param.strLoadFileType == "pfind")
	{	//read pfind results

	}
	else if (param.strLoadFileType == "maxquant")
	{	//read maxquant results
		CLoadMaxquantIO Proteinsio(flog, param);
		Proteinsio.mf_LoadProteinFasta(flog, param.m_strProteinFasterFilePath, param.m_eFastaType);
		//Proteinsio.mf_LoadPeptides(param.m_strPeptideFilePath);
		Proteinsio.mf_LoadPeptidesEithoutRedunPeptides(flog, param.m_strPeptideFilePath);
		Proteinsio.mf_LoadProteins(flog, param.m_strProteinFilePath);
		Proteinsio.mf_saveProteins(flog, param.m_strProteinsPath);
	}


	end = clock();
	cout << "\tLoading work using " << (end - begin) / CLOCKS_PER_SEC << " seconds." << endl;
	flog << "\tLoading work using " << (end - begin) / CLOCKS_PER_SEC << " seconds." << endl;
	flog.close();
	return 0;
}

