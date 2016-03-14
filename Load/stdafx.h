// stdafx.h : 标准系统包含文件的包含文件，
// 或是经常使用但不常更改的
// 特定于项目的包含文件
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


// TODO: 在此处引用程序需要的其他头文件
//#include"AbsoluteQuan.h"
//#include"stepwise.h"
//#include "BasicClass.h"
//#include "DataIO.h"
enum DataType{ MaxquantTpye = 0, pfindType = 1 };
typedef std::regex FastaType;
const int BUFFERLENGTH = 100000;
const int UniquePeptidesThreshold = 5;
