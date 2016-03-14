#pragma once

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>
#include<time.h>
#include<fstream>
#include<string>
#include<vector>
#include <algorithm>
#include<regex>
#include<iostream>
using namespace std;

//#define BUFFERLENGTH 2048; 
#define BUF_LENGTH 200000
const int FlankingRegionWidth = 15;
const int NumberOfFeature = 587;

// TODO: 
#include"AbsoluteQuan.h"
#include"stepwise.h"
typedef std::regex FastaType;
const int BUFFERLENGTH=20000; 
const int UniquePeptidesTrainThreshold=5;
const int UniquePeptidesTestThreshold = 0; //Thr threshold used to determine if calculate SPEAQInPeptides or not
const int UniquePeptidesCorrectThreshold = 0;// Thr threshold used to determine if correct peptide or not
