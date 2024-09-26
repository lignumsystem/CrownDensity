#include<ModeChangeFunctions.h>

void InsertMetaFiles(const string& globexpr, deque<string>& metafiles)
{
  glob_t glob_result;
  glob(globexpr.c_str(),GLOB_TILDE|GLOB_BRACE,NULL,&glob_result);
  cout << "Inserting MetaFiles"<<endl;
  for(unsigned int i=0;i<glob_result.gl_pathc;++i){
    string fname = glob_result.gl_pathv[i];
    cout << fname << endl;
    metafiles.push_back(fname);
  }
  cout << "Done" <<endl;
  sort(metafiles.begin(),metafiles.end(),less<string>());
  globfree(&glob_result);
}

void InsertModeChangeYears(const string& mode_years,deque<int>& modeyears)
{
  stringstream year_ss(mode_years);
  string year_s;
  cout << "Inserting mode change years" <<endl;
  while(getline(year_ss,year_s,',')){
    int year = atoi(year_s.c_str());
    cout << year <<endl;
    modeyears.push_back(year);
  }
  cout << "Done" <<endl;
  sort(modeyears.begin(),modeyears.end(),less<int>());
}

bool ModeChangeConsistency(deque<string>& metafiles,deque<int>& modeyears)
{
  if (metafiles.size() != (modeyears.size() + 1)){
    cout << "Growth Mode Change constistency error" <<endl;
    cout << "The number of MetaFiles " << metafiles.size()
	 << " must be the number of Mode Change Years " << modeyears.size() << " plus 1" <<endl;
    return false;
  }
  return true;
}
