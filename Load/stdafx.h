// stdafx.h : ��׼ϵͳ�����ļ��İ����ļ���
// ���Ǿ���ʹ�õ��������ĵ�
// �ض�����Ŀ�İ����ļ�
//

#pragma once

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>
#include<fstream>
#include<string>
#include<vector>
#include<map>
#include <algorithm> 
#include<time.h>
#include<regex>
#include<iostream>
using namespace std;

//#define BUFFERLENGTH 2048; 


// TODO: �ڴ˴����ó�����Ҫ������ͷ�ļ�
//#include"AbsoluteQuan.h"
//#include"stepwise.h"
//#include "BasicClass.h"
//#include "DataIO.h"
enum DataType{ MaxquantTpye = 0, pfindType = 1 };
typedef std::regex FastaType;
const int BUFFERLENGTH = 100000;
const int UniquePeptidesThreshold = 5;
