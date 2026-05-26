#include<cmath>
#include<stdlib.h>
#include<cstring>
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<algorithm>


using namespace std;

bool verbose = false;


class Extremum
{
    public:
    Extremum():pos(0),height(0.0),max(true) {}
    Extremum(unsigned int p, double h, bool m){
        pos = p;
        height = h;
        max = m;
    }
    ~Extremum(){}

    unsigned int pos;
    double height;
    bool max;

}; // class Extremum


bool operator<(const Extremum& a, const Extremum& b){

    return a.pos < b.pos;

} // operator<(const Extremum& a, const Extremum& b)


bool operator<(const Extremum& a, const int x){

    return a.pos < x;

} // operator<(const Extremum& a, const Extremum& b)


void findMinMax(const vector<unsigned int>& pos, vector<Extremum>& ex, const vector<double>& v, const int len, const int tol){
    
    ex.resize(0);
    for (unsigned int i = 500; i < v.size()-len-1; i++){
        if ((v[i] < v[i-2]) && (v[i] < v[i+2])){
            int num = 0;
            for (unsigned int j = i-len; j <= i+len; j++)
                if (v[j] > v[i])
                    num++;
            if (num >= 2*len-tol)
                ex.push_back(Extremum(pos[i], v[i], false));
        }
        if ((v[i] > v[i-2]) && (v[i] > v[i+2])){
            int num = 0;
            for (unsigned int j = i-len; j <= i+len; j++)
                if (v[j] < v[i])
                    num++;
            if (num >= 2*len-tol)
                ex.push_back(Extremum(pos[i], v[i], true));
        }            
    }
    sort(ex.begin(), ex.end());

}  // findMinMax


int maximumInRange(const vector<Extremum>& sorted_array, int from, int to) {

    const double lowest = -1.0E10;
    auto iter_geq = lower_bound(sorted_array.begin(), sorted_array.end(), to);

    Extremum best = *iter_geq;
    if (best.pos > to)
        best.height = lowest;
    int bestIndex = iter_geq - sorted_array.begin();
    while ((iter_geq >= sorted_array.begin()) && ((*iter_geq).pos >= from)){
        if (((*iter_geq).pos <= to) && (*iter_geq).max && ((*iter_geq).height > best.height)){
            best = *iter_geq;
            bestIndex = iter_geq - sorted_array.begin();
        }
        --iter_geq;
    }
    
    if ((best.pos <= to) && best.max && (best.height > lowest))
        return bestIndex;
    else
        return -1;

} // maximumInRange


void findChromatosomes(vector<unsigned int>& pos, vector<double>& sc, vector<double>& rsc,
                       vector<double>& sw, vector<double>& ss, vector<double>& se, unsigned int len,
                       unsigned int minDist, unsigned int minDistSE, unsigned int maxDistSE,
                       double thresHeight, int tol, bool strictSE, string& baseName){
    
    // extrema finding
    if (verbose)
        cerr << "Finding extrema ...\n";
    vector<Extremum> ex_sw, ex_sc, ex_st, ex_en;
    findMinMax(pos, ex_sw, sw, len, 0);
    findMinMax(pos, ex_sc, sc, len, 0);
    bool w_se;
    if ((ss.size() == 0) || (se.size() == 0))
        w_se = false;
    else{
        findMinMax(pos, ex_st, ss, len, tol);
        findMinMax(pos, ex_en, se, len, tol);
        w_se = true;
    }

    // candidate coverage peaks
    if (verbose){
        cerr << "Finding candidate peaks ...\n";
		int nWPSmax = 0;
		for (unsigned int i = 1; i < ex_sw.size()-1; i++)
			if (ex_sw[i].max)
				nWPSmax++;
		cerr << "  #WPS local maxima: " << nWPSmax << endl;
	}
    ofstream outCand(baseName + "_candidates_cover.txt");
    outCand << "position\theight\tnorm_coverage\tcoverage\tdiameter\n";
    vector<unsigned int> cand;
    vector<int> diameters;
    for (unsigned int i = 1; i < ex_sw.size()-1; i++){
        if (!ex_sw[i-1].max && ex_sw[i].max && !ex_sw[i+1].max &&
            (ex_sw[i].pos-ex_sw[i-1].pos>=minDist) && (ex_sw[i+1].pos-ex_sw[i].pos>=minDist)){

            // found bona fide WPS peak, check sequence coverage
            double coverage = sc[ex_sw[i].pos-pos[0]];
            if (coverage >= thresHeight){
                // estimate diameter
                auto iter_geq = lower_bound(ex_sc.begin(), ex_sc.end(), ex_sw[i].pos);
                auto iter = iter_geq;
                unsigned int right;
                if (!(*iter).max)
                    right = (*iter).pos;
                else{
                    right = ex_sw[i].pos+1000;
                    while ((iter < ex_sc.end()) && (*iter).max)
                        ++iter;
                    if ((iter < ex_sc.end()) && !(*iter).max)
                        right = (*iter).pos;
                }
                unsigned int left = ex_sw[i].pos-1000;
                iter = iter_geq-1;
                if (iter >= ex_sc.begin()){
                    if (!(*iter).max)
                        left = (*iter).pos;
                    else{
                        while ((iter >= ex_sc.begin()) && (*iter).max)
                            --iter;
                        if ((iter >= ex_sc.begin()) && !(*iter).max)
                            left = (*iter).pos;
                    }
                }
                int diameter = right-left;
                
                // output candidate peak and store diameter for subsequent selection
                outCand << ex_sw[i].pos << "\t" << ex_sw[i].height << "\t" << coverage << "\t" << rsc[ex_sw[i].pos-pos[0]] <<
                           "\t" << diameter << endl;
                cand.push_back(i);
                diameters.push_back(diameter);
            }
        }
    }
    outCand.close();
    if (verbose){
        cerr << "  Found " << cand.size() << " candidate WPS peaks." << endl;
        cerr << "Selecting ...\n";
    }
    
    // peak selection
    ofstream outSel(baseName + "_selection_cover.txt");
    outSel<< "position\theight\tnorm_coverage\tcoverage\tdiameter\n";
    unsigned int nSel = 0;
    for (unsigned j = 0; j < cand.size(); j++){
        unsigned int i = cand[j];
        int diameter = diameters[j];
        unsigned int center = ex_sw[i].pos;
        if (w_se){
            // check for adjacent maxima in starts and ends
            int ks = maximumInRange(ex_st, center-maxDistSE, center-minDistSE);
            int ke = maximumInRange(ex_en, center+minDistSE, center+maxDistSE);
            if ((ks >= 0) && (ke >= 0)){
                if (strictSE){
                    // also check starts and ends on the other sides, must be absent or 10% lower
                    int sks = maximumInRange(ex_st, center-minDistSE+1, center+minDistSE);
                    int ske = maximumInRange(ex_en, center-minDistSE, center+minDistSE-1);
                    if (((sks < 0) || (ex_st[ks].height > 1.1*ex_st[sks].height)) &&
                        ((ske < 0) || (ex_en[ks].height > 1.1*ex_en[ske].height))){
                        // good for strict SE criteria, write the peak location
                        outSel << center << "\t" << ex_sw[i].height << "\t" << sc[center-pos[0]] << "\t" << rsc[center-pos[0]] <<
                                  "\t" << diameter << endl;
                        nSel++;
                    }
                }
                else{
                    // all good for non strict SE criteria, write the peak location
                    outSel << center << "\t" << ex_sw[i].height << "\t" << sc[center-pos[0]] << "\t" << rsc[center-pos[0]] <<
                              "\t" << diameter << endl;
                    nSel++;
                }
            }
        }
        else{
            // every candidate WPS peak is selected
            outSel << center << "\t" << ex_sw[i].height << "\t" << sc[center-pos[0]] << "\t" << rsc[center-pos[0]] <<
                      "\t" << diameter << endl;
            nSel++;
        }
    }
    outSel.close();
    if (verbose){
        cerr << "Selected " << nSel << " candidates." << endl;
    
        // write the WPS, starts, and ends local maxima for information
        ofstream outWPS(baseName + "_wps_max.txt");
        outWPS << "position\theight\n";
        for (const auto & wpsm : ex_sw)
            if (wpsm.max)
                outWPS << wpsm.pos << "\t" << sw[wpsm.pos-pos[0]] << endl;
        outWPS.close();
        if (w_se){
            ofstream outS(baseName + "_starts_max.txt");
            outS << "position\theight\n";
            for (const auto & stm : ex_st)
                if (stm.max)
                    outS << stm.pos << "\t" << ss[stm.pos-pos[0]] << endl;
            outS.close();
            ofstream outE(baseName + "_ends_max.txt");
            outE << "position\theight\n";
            for (const auto & enm : ex_en)
                if (enm.max)
                    outE << enm.pos << "\t" << se[enm.pos-pos[0]] << endl;
            outE.close();
        }
    }
    
} // findChromatosomes


int main(int argc, char* argv[]){

    unsigned int len = 30;
    unsigned int minDist = 40;
    unsigned int minDistSE = 25;
    unsigned int maxDistSE = 120;
    int tol = 2;
    double thresHeight = 0.15;
    bool strictSE = false;
    string scFile, rscFile, swFile, ssFile, seFile, baseName;

    for (int i = 1; i < argc; i++){
        if ((strcmp(argv[i], "-help") == 0) || (strcmp(argv[i], "-h") == 0)){
            cerr << "Usage: findChromatosomes [options] -sc <filename> -rsc <filename> -sw <filename> -basename <filename>\n\nOptions are:\n";
            cerr << "-help\n";
            cerr << "-verbose                  print more details about the computation\n";
            cerr << "-strictse                 activate strict mode to check starts on the right of WPS peaks and end on the left as well\n";
            cerr << "-basename <filename>      base name for output files\n";
            cerr << "-sc <filename>            file with smoothed normalized position coverage\n";
            cerr << "-sw <filename>            file with smoothed WPS\n";
            cerr << "-ss <filename>            file with start position coverage\n";
            cerr << "-se <filename>            file with end position coverage\n";
            cerr << "-rsc <filename>           file with smoothed raw (not normalized) coverage\n";
            cerr << "-len <int>                number of bp befeore/after a positions to find an extremum, default=" << len << endl;
            cerr << "-tol <int>                extrema tolerance (more extreme than all but tol before/after positions), default=" << tol << endl;
            cerr << "-mindist <int>            minimum distance to the next/previous extremum, default=" << minDist << endl;
            cerr << "-mindistse <int>          minimum distance to closest start or end position coverage maximum, default=" << minDistSE << endl;
            cerr << "-maxdistse <int>          maximum distance to closest start or end position coverage maximum, default=" << maxDistSE << endl;
            cerr << "-thresheight <float>      height in the smoothed normalized coverage to determine chromatosome presence, default=" << thresHeight << endl;
            exit(0);
        }
        else if (strcmp(argv[i], "-verbose") == 0)
            verbose = true;
        else if (strcmp(argv[i], "-strictse") == 0)
            strictSE = true;
        else if ((strcmp(argv[i], "-sc") == 0) && (i < argc-1))
            scFile = string(argv[++i]);
        else if ((strcmp(argv[i], "-sw") == 0) && (i < argc-1))
            swFile = string(argv[++i]);
        else if ((strcmp(argv[i], "-ss") == 0) && (i < argc-1))
            ssFile = string(argv[++i]);
        else if ((strcmp(argv[i], "-se") == 0) && (i < argc-1))
            seFile = string(argv[++i]);
        else if ((strcmp(argv[i], "-rsc") == 0) && (i < argc-1))
            rscFile = string(argv[++i]);
        else if ((strcmp(argv[i], "-basename") == 0) && (i < argc-1))
            baseName = string(argv[++i]);
        if ((strcmp(argv[i], "-len") == 0) && (i < argc-1))
            len = atoi(argv[++i]);
        if ((strcmp(argv[i], "-tol") == 0) && (i < argc-1))
            tol = atoi(argv[++i]);
        if ((strcmp(argv[i], "-mindist") == 0) && (i < argc-1))
            minDist = atoi(argv[++i]);
        if ((strcmp(argv[i], "-mindistse") == 0) && (i < argc-1))
            minDistSE = atoi(argv[++i]);
        if ((strcmp(argv[i], "-maxdistse") == 0) && (i < argc-1))
            maxDistSE = atoi(argv[++i]);
        if ((strcmp(argv[i], "-thresheight") == 0) && (i < argc-1))
            thresHeight = atof(argv[++i]);
    }
    
    // read data
    vector<double> sc, sw, ss, se, rsc;
    vector<unsigned int> position;
    if (verbose)
        cerr << "Reading data ..." << endl;
    ifstream file;
    double pos;
    double val;
    file.open(scFile);
    while (file >> pos >> val){
        sc.push_back(val);
        position.push_back(static_cast<unsigned int>(pos));
    }
    file.close();
    if (verbose)
        cerr << "  read " << position.size() << " rows, last position=" << position[position.size()-1] << endl;
    file.open(rscFile);
    while (file >> pos >> val)
        rsc.push_back(val);
    file.close();
    if (verbose)
        cerr << "  read " << rsc.size() << " rows, last position=" << pos << endl;
    file.open(swFile);
    while (file >> pos >> val)
        sw.push_back(val);
    file.close();
    if (verbose)
        cerr << "  read " << sw.size() << " rows, last position=" << pos << endl;
    if ((ssFile.length() > 0) && (seFile.length() > 0)){
        if (verbose)
            cerr << "Reading starts and ends ..." << endl;
        file.open(ssFile);
        while (file >> pos >> val)
            ss.push_back(val);
        file.close();
        if (verbose)
            cerr << "  read " << ss.size() << " rows, last position=" << pos << endl;
        file.open(seFile);
        while (file >> pos >> val)
            se.push_back(val);
        file.close();
        if (verbose)
            cerr << "  read " << se.size() << " rows, last position=" << pos << endl;
    }
    
    // local normalization
    if (verbose)
        cerr << "Peak picking ..." << endl;
    findChromatosomes(position, sc, rsc, sw, ss, se, len, minDist, minDistSE, maxDistSE,
                      thresHeight, tol, strictSE, baseName);
    
    return 0;
    
} // main
