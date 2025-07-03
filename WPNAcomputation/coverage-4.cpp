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


bool verbose = false;
string chrSizeFileName = "pos_chromosomes.tsv";


bool readChrSizesAndAllocateMemory(string fn, bool w_startend, vector<vector<unsigned int>>& genome,
                                   vector<vector<unsigned int>>& fragStart, vector<vector<unsigned int>>& fragEnd)
{
    ifstream sf;
    sf.open(fn);
    if (!sf.is_open()) {
      cerr << "Failed to open " << fn << endl;
      return false;
    }

    unsigned int indice, n1, n2, size;
    string chr;
    genome.resize(24);
    if (w_startend){
      fragStart.resize(24);
      fragEnd.resize(24);
    }
    while (sf >> indice >> chr >> size >> n1 >> n2){
      genome[indice].resize(size+1);
      if (w_startend){
        fragStart[indice] = genome[indice];
        fragEnd[indice] = genome[indice];
      }
    }
    sf.close();
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


void usage()
{
  cerr << "Usage: coverage [options] -infolder <filename> -fileregexp <string> -outfolder <filename>\n\nOptions are:\n";
  cerr << "-help\n";
  cerr << "-verbose                         prints more details about the computation\n";
  cerr << "-sizerange <string>              a token indicating the fragment size range (chroma|mono|tf|all), default=chroma\n";
  cerr << "-withstartend                    also output fragment start and end position coverage\n";
  cerr << "-chromosome <string>             limits the computation and output to a specific chromosome\n";
  cerr << "-chrsizes <filename>             file containing the chromosome sizes, default is " << chrSizeFileName << " in the current folder\n";
  exit(0);

} // usage


int main(int argc, char* argv[]){
    
    bool w_startend = false;
    string sizeRange = "chroma";
    string chromosome, inFolder, outFolder, fileRegexp;

    for (int i = 1; i < argc; i++){
        if (strcmp(argv[i], "-help") == 0)
            usage();
        else if (strcmp(argv[i], "-withstartend") == 0)
            w_startend = true;
        else if (strcmp(argv[i], "-verbose") == 0)
            verbose = true;
        if ((strcmp(argv[i], "-sizerange") == 0) && (i < argc-1))
            sizeRange = string(argv[++i]);
        else if ((strcmp(argv[i], "-chromosome") == 0) && (i < argc-1))
            chromosome = string(argv[++i]);
        else if ((strcmp(argv[i], "-infolder") == 0) && (i < argc-1))
            inFolder = string(argv[++i]);
        else if ((strcmp(argv[i], "-outfolder") == 0) && (i < argc-1))
            outFolder = string(argv[++i]);
        else if ((strcmp(argv[i], "-fileregexp") == 0) && (i < argc-1))
            fileRegexp = string(argv[++i]);
    }
    
    // get fragment size range
    unsigned int minFragSize, maxFragSize;   
    if (sizeRange == "chroma"){
        minFragSize = 120;
        maxFragSize = 180;
        if (verbose)
            cerr << "Chromatosome fragment size range" << endl;
    }
    else if (sizeRange == "tf"){
        minFragSize = 35;
        maxFragSize = 119;
        if (verbose)
            cerr << "TF fragment size range" << endl;
    }
    else if (sizeRange == "all"){
        minFragSize = 0;
        maxFragSize = 10000;
        if (verbose)
            cerr << "All fragment size range" << endl;
    }
    else if (sizeRange == "mono"){
        minFragSize = 120;
        maxFragSize = 150;
        if (verbose)
            cerr << "Mononucleosome fragment size range" << endl;
    }
    else{
        cerr << "Fragment size range must be specified (chroma|mono|tf|all)." << endl;
        return 1;
    }
    
    // get the chromosome sizes to allocate memory
    if (verbose)
        cerr << "Allocating memory\n";
    vector<vector<unsigned int>> genome;
    vector<vector<unsigned int>> fragStart;
    vector<vector<unsigned int>> fragEnd;
    if (!readChrSizesAndAllocateMemory(chrSizeFileName, w_startend, genome, fragStart, fragEnd))
        return 1;

    // parse files and compute coverage
    string chr, strand;
    unsigned int start, end, qual;   
    vector<string> beds;
    regex pattern(fileRegexp);
    readDirectory(inFolder, beds, pattern);
    if (verbose)
        cerr<<"Directory read\n";
    for (int i = 0; i < beds.size(); i++){
        ifstream file;
        if (verbose)
            cerr << "  Parsing " << beds[i] << endl;
        file.open(inFolder + "/" + beds[i]);
        if (!file.is_open()) {
            cerr << "Failed to open the TSV file." << endl;
            return 1;
        }
        
        unsigned int nLines = 0;
        while (file >> chr >> start >> end >> qual >> strand) {
//         while (file >> chr >> start >> end >> strand >> qual) { // DSP format
            if (chr.length() < 4)
                chr = "chr"+chr;
            nLines++;
            unsigned int len = end-start;
            if ((qual > 30) && (len <= maxFragSize) && (len >= minFragSize)){
                if ((qual > 30) && (chr != "chrM") && 
                    (chr != "chrEBV") && (chr.find("_") == -1) && (chr.find("KI") == -1) &&
                    (chr.find("MT") == -1) && (chr.find("GL") == -1)){
                    unsigned int indice;
                    
                    if (chr == "chrX")
                        indice = 22;
                    else
                        if (chr == "chrY")
                            indice = 23;
                        else
                            indice = stoul(chr.substr(3,chr.length()))-1;
                        if ((indice < 0) || (indice > 23))
                         cerr << "Error: wrong indice " << indice << " line " << nLines << endl;
                        else
                        for (unsigned int j = start+1; j <= end; j++)
                            genome[indice][j]++;
                        if (w_startend){
                            fragStart[indice][start+1]++;
                            fragEnd[indice][end]++;
                        }
                }
            }
        }
        
        file.close();
        if (verbose)
            cerr << "    " << nLines << " lines processed\n";
    }
    
    // output one coverage file per chromosome
    for (int j=0; j < 24; j++){
        string fname;
        if (j < 22)
            chr = "chr" + to_string((j+1));
        else
            if (j == 22)
                chr = "chrX";
            else
                chr = "chrY";
        if ((chromosome.length() == 0) || (chromosome == chr)){
            fname = outFolder + "/" + chr;
            ofstream outFile(fname+"_coverage.tsv");
            for (const auto & k : genome[j]) outFile << k << "\n";
            outFile.close();
            if (w_startend){
                ofstream outFileS(fname+"_starts.tsv");
                for (const auto & k : fragStart[j]) outFileS << k << "\n";
                outFileS.close();
                ofstream outFileE(fname+"_ends.tsv");
                for (const auto & k : fragEnd[j]) outFileE << k << "\n";
                outFileE.close();
            }
            if (verbose)
                cerr << "  Finished chromosome " << j+1 << endl;
        }
    }
    
    return 0;
    
} // main
