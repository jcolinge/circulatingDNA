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

bool verbose = false;


vector<int> slidingMax(const vector<int> &array, const int window) {

	vector<int> result;
    int max = array[0];
    deque<int> deque;
	int countdown = window - 1;
	typename vector<int>::const_iterator tail = array.cbegin();
	for (const int &val : array) {
		while (!deque.empty() && (val > deque.back()))
			deque.pop_back();
		deque.push_back(val);
		
		if (countdown > 0)
			countdown--;
		else {
			result.push_back(deque.front());
			if (*tail == deque.front())
				deque.pop_front();
			++tail;
		}
	}
	
	return result;
}


vector<int> slidingMin(const vector<int> &array, const int window) {

	vector<int> result;    
	deque<int> deque;
	int countdown = window - 1;
	typename vector<int>::const_iterator tail = array.cbegin();
	for (const int &val : array) {
		while (!deque.empty() && (val < deque.back()))
			deque.pop_back();
		deque.push_back(val);
		
		if (countdown > 0)
			countdown--;
		else {
			result.push_back(deque.front());
			if (*tail == deque.front())
				deque.pop_front();
			++tail;
		}
	}
	
	return result;
}


vector<double> normalize(const vector<int> &v, const vector<int> &mins,
                         const vector<int> &maxes, const unsigned int winSize){
  
    unsigned int hwl = trunc(0.5*winSize);
    vector<double> nv;
    double n, delta;
    for (unsigned int i = 0; (i <= hwl) && (i < v.size()); i++){
        delta = maxes[0]-mins[0];
        if (delta > 0)
            nv.push_back(1.0*(v[i]-mins[0])/delta);
        else
            nv.push_back(0.0);
    }
    for (unsigned int i = hwl+1; i < v.size()-hwl; i++){
        delta = maxes[i-hwl]-mins[i-hwl];
        if (delta > 0)
            nv.push_back(1.0*(v[i]-mins[i-hwl])/delta);
        else
            nv.push_back(0.0);
    }
    for (unsigned int i = v.size()-hwl; i < v.size(); i++){
        delta = maxes[maxes.size()-1]-mins[mins.size()-1];
        if (delta > 0)
            nv.push_back(1.0*(v[i]-mins[mins.size()-1])/delta);
        else
            nv.push_back(0.0);
    }
    
    return nv;
}


int main(int argc, char* argv[]){

  if ((strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "-help") == 0)){
      cerr << "Usage: normalize [-h] start_pos end_pos win_size input_file" << endl;
      return 0;
  }
    
  unsigned int from = atoi(argv[1]);
  int tmp = atoi(argv[2]);
  unsigned int to;
  if (tmp <= 0)
      to = 1e9;
  else
      to = tmp;
  unsigned int winSize = atoi(argv[3]);
  if ((winSize % 2) == 0){
      cerr << "winSize must be odd." << endl;
      return 1;
  }

  // read data
  if (verbose)
      cerr << "Reading data ..." << endl;
  ifstream file;
  file.open(string(argv[4]));
  if (!file.is_open()) {
      cerr << "Failed to open the file." << endl;
      return 1;
  }
  int val;
  unsigned int line = 0;
  vector<int> c;
  while (file >> val) {
      if ((line >= from) && (line <= to))
          c.push_back(val);
      line++;
      if (line > to)
          break;
  }
  file.close();

  // local normalization
  if (verbose)
      cerr << "Normalizing ..." << endl;
  vector<int> mins = slidingMin(c, winSize);
  vector<int> maxes = slidingMax(c, winSize);
  vector<double> nc = normalize(c, mins, maxes, winSize);
  
  // output
  if (verbose)
      cerr << "Writing ..." << endl;
  for (unsigned int i = 0; i < c.size(); i++)
      printf("%i\t%.3f\n", from+i, nc[i]);
  
  return 0;
    
} // main
