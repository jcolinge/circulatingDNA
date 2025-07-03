#include<stdlib.h>
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


string chrSizeFileName = "pos_chromosomes.tsv";
vector<vector<unsigned int>> genome, randGenome;

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


bool readChrSizes(string fn) {

    ifstream sf;
    sf.open(fn);
    if (!sf.is_open()) {
      cerr << "Failed to open " << fn << endl;
      return false;
    }

    unsigned int indice, n1, n2, size;
    string chr;
    genome.resize(24);
    while (sf >> indice >> chr >> size >> n1 >> n2){
      genome[indice].resize(size+1);
      fill(genome[indice].begin(),genome[indice].end(),0);
    }
    sf.close();
    randGenome = genome;
    return true;

} // readChrSizes


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

    int tol = 20;
    int maxLength = 1000;
    
    if ((strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "-help") == 0)){
        cerr << "Usage: computeLengthDistrib [-h] k-stack_file input_folder input_file_regexp baseName" << endl;
        return 0;
    }
    
    // get the chromosome sizes to allocate memory
    cerr << "Allocating memory ...\n";
    if (!readChrSizes(chrSizeFileName))
        return 1;
    
    // read the centromere coordinates (genome coordinate convention)
    cerr << "Read centromere positions." << endl;
    vector<FragmentPosition> centromere;
    centromere.resize(23);
    ifstream centroFile("centromere-region.txt");
    string chr;
    int start, end;
    unsigned int indice;
    while (centroFile >> chr >> start >> end)
        if ((chr != "Y") & (chr != "chr")){
            if (chr == "X")
                indice = 22;
            else
                indice = stoul(chr)-1;
            int len = end-start+1;
            centromere[indice] = FragmentPosition(start,end,len);
        }
    centroFile.close();

    // read the k-stack coordinates (genome coordinate convention)
    ifstream atlasFile(argv[1]);
    cerr << "  Parsing the atlas..." << endl;
    int good = 0;
    string strand;
    unsigned int dummy, kst;
    while (atlasFile >> chr >> start >> end >> dummy >> strand >> kst)
        if (chr != "Y"){
            if (chr == "X")
                indice = 22;
            else
                indice = stoul(chr)-1;
            if ((end <= centromere[indice].first) || (start >= centromere[indice].last)){
                //correct posisition
                for (unsigned int i = start; (i < genome[indice].size()) && (i <= end); i++)
                    genome[indice][i] = 1;
                //generate a random position
                int rstart, rend;
                do{
                    rstart = start+rand()%100000-50000;
                    rend = rstart+end-start;
                } while ((rstart >= centromere[indice].first && rend <= centromere[indice].last) ||
                    (rstart < centromere[indice].first && rend > centromere[indice].first) ||
                    (rstart < centromere[indice].last && rend > centromere[indice].last));
                for (unsigned int i = rstart; (i < randGenome[indice].size()) && (i <= rend); i++)
                    randGenome[indice][i] = 1;                
            }
            good++;
        }
    atlasFile.close();
    cerr << good << " used entries." << endl;
  
    // read BED files and compute length distributions (BED file coordinate convention)
    vector<string> beds;
    regex pattern(argv[3]);
    readDirectory(argv[2], beds, pattern);
    vector<unsigned int> lenAll, lenIn, lenOut, lenRandIn;
    lenAll.resize(maxLength+1);
    fill(lenAll.begin(),lenAll.end(),0);
    lenIn = lenAll;
    lenOut = lenAll;
    lenRandIn = lenAll;
    for (int i = 0; i < beds.size(); i++){
        ifstream file;
        cerr << "  Parsing " << beds[i] << endl;
        file.open(string(argv[2]) + "/" + beds[i]);
        if (!file.is_open()) {
            cerr << "Failed to open the TSV file." << endl;
            return 1;
        }
        
        unsigned int nLines = 0;
        string strand;
        int qual;
//        while (file >> chr >> start >> end >> strand >> qual) { // DSP format
//            chr = "chr"+chr;
        while (file >> chr >> start >> end >> qual >> strand) { // FinaleDB format
            nLines++;
            unsigned len = end-start;
            if ((qual > 30) && (chr != "chrY") && (chr != "chrM") && 
                (chr != "chrEBV") && (chr.find("_") == -1) && (chr.find("KI") == -1) &&
                (chr.find("MT") == -1) && (chr.find("GL") == -1) && (len <= maxLength)){
                if (chr == "chrX")
                    indice = 22;
                else
                    indice = stoul(chr.substr(3,chr.length()))-1;
                lenAll[len]++;
                int sum = 0;
                for (unsigned int j = start+1; (j < genome[indice].size()) && (j < end); j++)
                    sum += genome[indice][j];
                if ((sum >= 167) || (sum >= end-start-tol))
                    lenIn[len]++;
                else
                    lenOut[len]++;
                sum = 0;
                for (unsigned int j = start+1; (j < randGenome[indice].size()) && (j < end); j++)
                    sum += randGenome[indice][j];
                if ((sum >= 167) || (sum >= end-start-tol))
                    lenRandIn[len]++;
            }
        }
        file.close();
        cerr << "    " << nLines << " lines processed\n";
    }
    
    // output results
    string baseName(argv[4]);
    ofstream outAll(baseName + "_all.txt");
    outAll << "length\tfreq\n";
    for (unsigned int i = 1; i < lenAll.size(); i++)
        outAll << i << "\t" << lenAll[i] << endl;
    outAll.close();
    ofstream outIn(baseName + "_in.txt");
    outIn << "length\tfreq\n";
    for (unsigned int i = 1; i < lenIn.size(); i++)
        outIn << i << "\t" << lenIn[i] << endl;
    outIn.close();
    ofstream outOut(baseName + "_out.txt");
    outOut << "length\tfreq\n";
    for (unsigned int i = 1; i < lenOut.size(); i++)
        outOut << i << "\t" << lenOut[i] << endl;
    outOut.close();
    ofstream outRandIn(baseName + "_randin.txt");
    outRandIn << "length\tfreq\n";
    for (unsigned int i = 1; i < lenRandIn.size(); i++)
        outRandIn << i << "\t" << lenRandIn[i] << endl;
    outRandIn.close();

    return 0;

} // main
