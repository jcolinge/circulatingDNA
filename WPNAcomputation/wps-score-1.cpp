#include<stdlib.h>
#include<cstring>
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<cmath>
#include<algorithm>
#include<regex>
#include<sys/types.h>
#include<dirent.h>

using namespace std;

bool verbose = true;


class FragmentPosition
{
public:
  FragmentPosition():first(-1),last(-1),length(-1) {}
  FragmentPosition(int f, int l, int len){
    first = f;
    last = l;
    length = len;
  }
  ~FragmentPosition(){}

  int first, last, length;

}; // class FragmentPosition


bool operator<(const FragmentPosition& a, const FragmentPosition& b){

  return (a.first < b.first) || ((a.first == b.first) && (a.last < b.last));

} // operator<(const FragmentPosition& a, const FragmentPosition& b)    


void readBEDFile(ifstream& file, vector<FragmentPosition>& frag, int minQual,
                 string chromosome, unsigned int from, unsigned int to,
                 int minLen, int maxSize){

  char dir;
  string line, chr;
  unsigned int start, end, qual;
    
  // read
  unsigned int nLines = 0;
  unsigned int good = 0;
  while (file >> chr >> start >> end >> qual >> dir) {
//   string strand;
//   while (file >> chr >> start >> end >> strand >> qual) { // DSP format
//     chr = "chr"+chr;
    nLines++;
    int len = end-start;
    if ((qual >= minQual) && (len >= minLen) && (len <= maxSize) &&
        (chr == chromosome) && (end <= to) && (start >= from)){
      frag.push_back(FragmentPosition(start+1,end,len));
      good++;
    }
  }
  if (verbose){
    cerr << nLines << " lines processed.\n";
    cerr << good << " lines with QUAL >= " << minQual << " on " << chromosome << 
                    " and fragment length >= " << minLen << " and <= " << maxSize << endl;
  }
  
} // readBEDFile


void computeWPS(vector<FragmentPosition> &frag, unsigned int from,
                unsigned int to, int winSize){

  if (to == 0)
    for (auto c = frag.begin(); c != frag.end(); ++c)
      if ((*c).last > to)
        to = (*c).last;
  int whl = trunc(0.5*winSize);
  

  auto c = frag.begin();
  for (int pos = from; pos <= to; pos++){
    while (((*c).first+winSize <= pos-whl) && (c != frag.end()))
      ++c;
    int wps = 0;
    for (auto check = c; (check != frag.end()) && ((*check).first <= pos+whl); ++check){
      if ((*check).first <= pos-whl)
        if ((*check).last >= pos+whl)
          wps++;
        else
          wps--;
      else
        if (((*check).first <= pos+whl) && ((*check).last >= pos+whl))
          wps--;
    }
    printf("%i\t%i\n", pos, wps);
  }
    
} // computeWPS


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

  if ((strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "-help") == 0)){
    cerr << "Usage: wps-score [-h] chromosome start_pos end_pos input_folder input_file_regexp" << endl;
    return 0;
  }
    
  string chromosome(argv[1]);
  unsigned int from = atoi(argv[2]);
  unsigned int to = atoi(argv[3]);
  if (verbose)
    cerr << "Computing WPS from positions " << from << " to " << to << " on chromosome " << chromosome << endl;

  // loop over BED files
  vector<FragmentPosition> fPos;
  vector<string> beds;
  regex pattern(argv[5]);
  readDirectory(argv[4], beds, pattern);
  for (int i = 0; i < beds.size(); i++){
    // open BED file
    ifstream file;
    if (verbose)
      cerr << "  Parsing " << beds[i] << endl;
    file.open(string(argv[4]) + "/" + beds[i]);
    if (!file.is_open()) {
      cerr << "Failed to open the TSV file." << endl;
      return 1;
    }
  
    // read BED file
    readBEDFile(file, fPos, 30, chromosome, from, to, 120, 180);
    file.close();
  }

  // sort fragments
  sort(fPos.begin(), fPos.end());


  // compute WPS scores
  if (verbose)
    cerr << "Computing WPS ...\n";
  if (fPos.size() >= 100)
    computeWPS(fPos, from, to, 120);
  else
    if (verbose)
      cerr << "Not enough data\n";

  return 0;
    
} // main
