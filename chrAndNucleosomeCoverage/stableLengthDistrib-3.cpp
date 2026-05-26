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
vector<vector<unsigned int>> stable, randStable;

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
    stable.resize(24);
    while (sf >> indice >> chr >> size >> n1 >> n2)
        stable[indice].resize(size+1,0);
    sf.close();
    randStable = stable;
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

    int maxLength = 1000;
    int tol = 20;
    int hw = 100;
    
    if ((strcmp(argv[1], "-h") == 0) || (strcmp(argv[1], "-help") == 0)){
        cerr << "Usage: stableLengthDistrib [-h] reference_file input_folder input_file_regexp baseName" << endl;
        return 0;
    }
    
    // get the chromosome sizes to allocate memory
    cerr << "Allocating memory ...\n";
    if (!readChrSizes(chrSizeFileName))
        return 1;
    
    // read the centromere coordinates (genome coordinate convention, chr Y regions are hardcoded)
    cerr << "Read centromere positions" << endl;
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

    // read the reference table (stable chromatosome center positions)
    ifstream atlasFile(argv[1]);
    cerr << "Parsing the reference position table" << endl;
    int pos;
    double cover;
    int nLines = 0;
    string header;
    atlasFile >> header >> header; // skip column headers
    while (atlasFile >> chr >> pos){
        nLines++;
        if (chr == "X")
            indice = 22;
        else if (chr == "Y")
            indice = 23;
        else
            indice = stoul(chr)-1;
            
        // stable posisition
        int from = pos-hw;
        if (from < 0)
            from = 0;
        for (unsigned int i = from; (i < stable[indice].size()) && (i <= pos+hw); i++)
            stable[indice][i] = 1;
            
        // generate a random stable position
		int rpos;
		if (chr == "Y"){
            do{
                rpos = pos+rand()%1000000-500000;
            } while ((rpos-hw <= 2781000) || (rpos+hw > 26673000) ||
            (rpos+hw < 18866000 && rpos-hw > 17457000) ||
            (rpos+hw < 26310000 && rpos-hw > 22378000) ||
            (rpos < 0) || (rpos >= randStable[indice].size()));
		}
		else{
            do{
                rpos = pos+rand()%1000000-500000;
            } while ((rpos-hw >= centromere[indice].first && rpos+hw <= centromere[indice].last) ||
            (rpos-hw < centromere[indice].first && rpos+hw > centromere[indice].first) ||
            (rpos-hw < centromere[indice].last && rpos+hw > centromere[indice].last) ||
            (rpos < 0) || (rpos >= randStable[indice].size()));
		}
        from = rpos-hw;
        if (from < 0)
            from = 0;
        for (unsigned int i = from; (i < randStable[indice].size()) && (i <= rpos+hw); i++)
            randStable[indice][i] = 1;
    }
    atlasFile.close();    
    cerr << "  loaded " << nLines << " lines\n";
  
    // read BED files and compute length distributions (BED file coordinate convention)
    vector<string> beds;
    regex pattern(argv[3]);
    readDirectory(argv[2], beds, pattern);
    vector<unsigned int> lenStable, lenRandStable;
    lenStable.resize(maxLength+1);
    fill(lenStable.begin(),lenStable.end(),0);
    lenRandStable = lenStable;
    for (int i = 0; i < beds.size(); i++){
        ifstream file;
        cerr << "Parsing " << beds[i] << endl;
        file.open(string(argv[2]) + "/" + beds[i]);
        if (!file.is_open()) {
            cerr << "Failed to open the TSV file" << endl;
            return 1;
        }
        
        unsigned int nLines = 0;
        int qual;
        string strand;
//         while (file >> chr >> start >> end >> strand >> qual) { // DSP format
//             chr = "chr"+chr;
        while (file >> chr >> start >> end >> qual >> strand) {
            nLines++;
            unsigned int len = end-start;
            if ((qual > 30) && (chr != "chrM") && 
                (chr != "chrEBV") && (chr.find("_") == -1) && (chr.find("KI") == -1) &&
                (chr.find("MT") == -1) && (chr.find("GL") == -1) && (len <= maxLength)){
                if (chr == "chrX")
                    indice = 22;
				else if (chr == "chrY")
					indice = 23;
                else
                    indice = stoul(chr.substr(3,chr.length()))-1;
                
                // overlap with a stable chromatosome
                int sum = 0;
                for (unsigned int j = start+1; (j < stable[indice].size()) && (j < end); j++)
                    sum += stable[indice][j];
                if ((sum >= 2*hw+1) || (sum >= end-start-tol))
                    lenStable[len]++;
                
                // overlap with a random position
                sum = 0;
                for (unsigned int j = start+1; (j < randStable[indice].size()) && (j < end); j++)
                    sum += randStable[indice][j];
                if ((sum >= 2*hw+1) || (sum >= end-start-tol))
                    lenRandStable[len]++;
            }
        }
        file.close();
        cerr << "  " << nLines << " lines processed\n";
    }
    
    // output results
    string baseName(argv[4]);
    ofstream outIn(baseName + "_stable.txt");
    outIn << "length\tfreq\n";
    for (unsigned int i = 1; i < lenStable.size(); i++)
        outIn << i << "\t" << lenStable[i] << endl;
    outIn.close();
    ofstream outRandIn(baseName + "_rndstable.txt");
    outRandIn << "length\tfreq\n";
    for (unsigned int i = 1; i < lenRandStable.size(); i++)
        outRandIn << i << "\t" << lenRandStable[i] << endl;
    outRandIn.close();

    return 0;

} // main
