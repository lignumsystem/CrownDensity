#ifndef MODECHANGEFUNCTIONS_H
#define MODECHANGEFUNCTIONS_H
#include<glob.h>
#include<iostream>
#include<string>
#include<sstream>
#include<deque>
using namespace std;
///\brief Insert MetaFiles from the command line. Glob expressions are accepted.
///\param globexpr The Glo expression denoting MetaFiles needed
///\param metafiles The queue of MetaFiles to used in the simulations  
void InsertMetaFiles(const string& globexpr, deque<string>& metafiles);
///\brief Insert comma separated list of mode change years from the command line
///\param years The comma separated list of years
///\modeyears The queue of Mode Change Years
void InsertModeChangeYears(const string& years, deque<int>& modeyears);
///\brief Check the consistency between number of MetaFiles and Growth Mode Years
///\return true If the number of MetaFiles the is number of GrowthModeYears+1, false otherwise
///\note The first MetaFile is consumed in the first tree initialization.
bool ModeChangeConsistency(deque<string>& metafiles,deque<int>& modeyears);
#endif
