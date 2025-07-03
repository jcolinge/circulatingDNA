#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<regex>
#include<cmath>
#include<cstring>
#include<sys/types.h>
#include<dirent.h>

using namespace std;


void readDirectory(const string name, vector<string>& v, const regex& pattern)
{
    
    DIR* dirp = opendir(name.c_str());
    struct dirent * dp;
    while ((dp = readdir(dirp)) != NULL)
      if (regex_search(dp->d_name, pattern))
        v.push_back(dp->d_name);
    closedir(dirp);

} // readDirectory


int main(int argc, char* argv[]){
    
  bool w_startend = false;

  if ((strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "-help") == 0)){
    cerr << "Usage: blockCoverage [-h] block_size input_folder input_file_regexp output_folder" << endl;
    return 0;
  }
    
  unsigned int bs = atoi(argv[1]);

  // parse files and compute coverage
  vector<string> beds;
  regex pattern(argv[3]);
  readDirectory(argv[2], beds, pattern);
  vector<int> block;
  for (int i = 0; i < beds.size(); i++){
    ifstream file;
    cerr << "  Processing " << beds[i] << endl;
    file.open(string(argv[2]) + "/" + beds[i]);
    if (!file.is_open()) {
      cerr << "Failed to open the TSV file." << endl;
      return 1;
    }
      
    string outFN = string(argv[4]) + "/" + "blocks-" + beds[i];
    cerr << "   output in " << outFN << endl;
    ofstream outFile(outFN);
    if (!outFile.is_open()) {
      cerr << "Failed to create the " << outFN << " file." << endl;
      return 1;
    }
   
    unsigned int c;
    unsigned int pos = 0;
    block.resize(0);
    while (file >> c){
      if (block.size() == bs){
        unsigned int tot = 0;
        unsigned int tot2 = 0;
        for (unsigned int i = 0; i < bs; i++){
          tot += block[i];
          tot2 += block[i]*block[i];
        }
        double avg = 1.0*tot/bs;
        double sd = sqrt(1.0*tot2/bs - avg*avg);
        double cv = 0.0;
        if (avg > 0)
          cv = sd/avg;
        outFile << pos-bs << "\t" << avg << "\t" << cv << endl;
        block.resize(0);
      }
      block.push_back(c);
      pos++;
    }
    if (block.size() > 0){
      unsigned int tot = 0;
      unsigned int tot2 = 0;
      for (unsigned int i = 0; i < block.size(); i++){
        tot += block[i];
        tot2 += block[1]*block[i];
      }
      double avg = 1.0*tot/block.size();
      double sd = sqrt(1.0*tot2/bs - avg*avg);
      double cv = 0.0;
      if (avg > 0)
        cv = sd/avg;
      outFile << pos-block.size() << "\t" << avg << "\t" << cv << endl;
    }
    outFile.close();
  }
  
  return 0;

} // main
