/*
*
* Command-line client for the KM-config
*
*
* An algorithm for finding multiple core-periphery pairs in networks
*
*
* Core-periphery structure requires something else in the network
* Sadamori Kojaku and Naoki Masuda
* Preprint arXiv:1710.07076
* 
*
* Please do not distribute without contacting the authors.
*
*
* AUTHOR - Sadamori Kojaku
*
*
* DATE - 11 Oct, 2017
*/

#include "../lib/km_config.h"
#include "../lib/km_config.cpp"
#include <string>
#include <stdio.h>
#include <unistd.h>

char delimiter = ' '; // the delimiter for the [input-file] and [output-file]

char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}


void split(const string& s, char c,
    vector<string>& v);

void readEdgeTable(string filename, vector<vector<int> >& A, vector<vector<double>>& W, int& N, double& M, int& edgeNum, bool& isunweighted, bool& isSelfLoop, int startIndex);
//void readEdgeTable(string filename, vector<vector<int> >& A, int& N);

void writeLabels(const string filename, const vector<int>& c, const vector<bool>& x, const vector<double> p_values, const double pval, int startIndex);

void usage();

std::string myreplace(std::string &s,
                      const std::string &toReplace,
                      const std::string &replaceWith)
{
    return(s.replace(s.find(toReplace), toReplace.length(), replaceWith));
}


int main(int argc, char* argv[])
{
    if (argc <= 2) {
        usage();
        return 1;
    }

    string linkfile = argv[1];
    string outputfile = argv[2];

    double alpha = 1.0;
    int num_of_runs = 10;
    int num_of_rand_nets = 500;
    int startIndex = 1;
   
    int opt;
    opterr = 0;
    string tmp;
   
    if(cmdOptionExists(argv, argv + argc, "-h")){
        usage();
    }
    if(cmdOptionExists(argv, argv + argc, "-a")){
        alpha = atof( getCmdOption(argv, argv + argc, "-a") );
    }
    if(cmdOptionExists(argv, argv + argc, "-r")){
        num_of_runs = atoi( getCmdOption(argv, argv + argc, "-r") );
    }
    if(cmdOptionExists(argv, argv + argc, "-l")){
        num_of_rand_nets = atoi( getCmdOption(argv, argv + argc, "-l") );
    }
    if(cmdOptionExists(argv, argv + argc, "-i")){
        startIndex = atoi( getCmdOption(argv, argv + argc, "-i") );
    }
    if(cmdOptionExists(argv, argv + argc, "-d")){
        tmp = getCmdOption(argv, argv + argc, "-d");
	if(!tmp.compare("\\t")){
		tmp = myreplace( tmp, "\\t", "\t");
		delimiter = *tmp.c_str();
	}else{
		delimiter = *tmp.c_str();
	}
    }

    cout << "===================" << endl;
    cout << "# input file:  "<< linkfile<< endl;
    cout << "# output file: "<< outputfile<< endl;
    cout << "" << endl;
  
   
    /* Read file */
    cout << "# Status: "<< endl;
    cout << "   Reading input file..."<<endl;
    
    int N;
    vector<vector<int>> A;
    vector<vector<double>> W;
    bool isSelfLoop; 
    bool isUnweighted;
    double M = 0;
    int edgeNum = 0;
    readEdgeTable(linkfile, A, W, N, M, edgeNum, isUnweighted, isSelfLoop, startIndex);

    cout << "   Number of nodes: "<< N << endl;
    cout << "   Number of edges: "<< edgeNum << endl;
    cout << "" << endl;

    /* Run the KM algorithm */
    cout << "   Seeking core-periphery pairs..."<<endl;
    cout << "      - Number of runs: "<< num_of_runs<< endl;
    vector<int> c(N);
    vector<bool> x(N);
    double Q;
    vector<double> q;
    srand(time(NULL));
    mt19937_64 mtrnd = init_random_number_generator();
    km_config_label_switching(A, W, num_of_runs, c, x, Q, q, mtrnd);
    cout <<"   end"<<endl<<endl;

    /* Statistical test */
    int K = q.size();
    vector<double> p_values(K);
    fill(p_values.begin(), p_values.end(), 0.0);
    double corrected_alpha =1.0 - pow(1.0 - alpha, 1.0 / (double)K); // Sidak correction.
    if (alpha < 1.0) {
    	cout << "   Significance test..."<<endl;
    	cout << "      - Number of random networks: "<< num_of_rand_nets<< endl;
    	cout << "      - Significance level: "<< alpha<< endl;
    	cout << "      - Corrected-significance level: "<< corrected_alpha<< endl;
    	cout << "      - Number of core-periphery pairs under testing: "<< K<< endl;
        estimate_statistical_significance(A, W, c, x, num_of_runs, num_of_rand_nets, p_values);
       	cout <<"   end"<<endl<<endl;
    }
    cout << "   "<<endl;
    int Ksig = 0;
    for(int i = 0; i < K; i++){
	if(p_values[i]<=corrected_alpha) Ksig++;
    }
    
    cout << "# Results: "<< endl;
    cout << "   Number of significant core-periphery pairs: "<< Ksig<<endl;
    cout << "   Number of insignificant core-periphery pairs: "<< K-Ksig<<endl;
    cout << "   Save to "<<outputfile<<"..."<<endl;
    cout <<"   end"<<endl;
    cout << "===================" << endl;

    /* Save results */
    writeLabels(outputfile, c, x, p_values, corrected_alpha, startIndex);
    return 0;
}

void usage()
{

    cout << endl
         << "\e[1mNAME:\e[0m" << endl
         << endl;
    cout << "  km_config - Implementation of the KM-config algorithm." << endl
         << endl
         << endl;

    cout << "\e[1mUSAGE:\e[0m" << endl
         << endl;
    cout << "  \e[1mkm_config\e[0m [input-file] [output-file] [options]" << endl
         << endl
         << endl;
    cout << "\e[1mDESCRIPTION:\e[0m" << endl
         << endl;
    cout << "  \e[1mkm_config\e[0m seeks multiple core-periphery pairs in the network given by [input-file] and" << endl;
    cout << "  saves the detected core-periphery pairs in [output-file]." << endl
         << endl;

    cout << "  \e[1m[input-file]\e[0m" << endl;
    cout << "    The file should contain a list of edges (space-separated)." << endl;
    cout << "    The first and second columns represent the IDs of the two nodes forming an edge." << endl;
    cout << "    The third column represent the weight of edgde between the two nodes." << endl;
    cout << "    If the third column is not provided, then the weight is set to 1" << endl;
    cout << "    The node's ID is assumed to start from 1. (you can change this by setting option, e.g., -i 0)" << endl
         << endl;

    cout << "  \e[1m[output_file]\e[0m" << endl;
    cout << "    The first column represents the node's ID." << endl;
    cout << "    The second column represents the index of the core-periphery pair to which each node belongs." << endl;
    cout << "    The third column indicates whether each node is a core node (= 1) or a peripheral node (= 0)." << endl;
    cout << "    The fourth column indicates whether each node belongs to a significant core-periphery pair (= 1) or not (= 0)." << endl
         << endl
         << endl;
    cout << "\e[1mOPTIONS:\e[0m" << endl
         << endl;
    cout << "  \e[1m-r R\e[0m  Run the KM algorithm R times. (Default: 10)" << endl
         << endl
         << "  \e[1m-a ALPHA\e[0m  Set the significance level before the Šidák correction to ALPHA. (Default: 1)."<< endl
         << "            If this option is not set, the statistical test is not carried out."<< endl
         << endl
         << "  \e[1m-l NUM\e[0m  Set the number of randomised networks to NUM. (Default: 500)" << endl
         << endl
         << "  \e[1m-i NUM\e[0m  The node's ID starts from NUM. (Default: 1)" << endl
         << endl
         << "  \e[1m-d \"D\"\e[0m  Change the delimiter for [input-file] and [output-file] to D. (Default: space)" << endl
         << endl
         << endl;
}

void split(const string& s, char c,
    vector<string>& v)
{
    string::size_type i = 0;
    string::size_type j = s.find(c);

    while (j != string::npos) {
        v.push_back(s.substr(i, j - i));
        i = ++j;
        j = s.find(c, j);

        if (j == string::npos)
            v.push_back(s.substr(i, s.length()));
    }
}

void readEdgeTable(string filename, vector<vector<int> >& A, vector<vector<double>>& W, int& N, double& M, int& edgeNum, bool& isUnweighted, bool& isSelfLoop, int startIndex)
{

    std::ifstream ifs(filename);
    vector<int> edgeList;
    vector<double> wList;
    string str;
    N = 0;
    edgeList.clear();
    while (getline(ifs, str)) {
        vector<string> v;
        split(str, delimiter, v);

        int sid = stoi(v[0]) - startIndex;
        int did = stoi(v[1]) - startIndex;
	double w = 1;
	
	if(v.size()>2){
        	w = stof(v[2]);
	}

        if (N < sid)
            N = sid;
        if (N < did)
            N = did;
        
        edgeList.push_back(sid);
        edgeList.push_back(did);
        wList.push_back(w);

    }
    N = N + 1;
   
    vector<vector<int>>tmp(N, vector<int>(0));
    vector<vector<double>>tmp2(N, vector<double>(0));
    A = tmp; 
    W = tmp2; 

   
    int wid = 0;
    int edgeListSize = edgeList.size();
    M = 0; 
    isSelfLoop = false; 
    isUnweighted = true; 
    for (int i = 0; i < edgeListSize; i += 2) {
        int sid = edgeList[i];
        int did = edgeList[i + 1];
	double w = wList[wid];

	wid++;
	
	if(sid == did){
        	A[sid].push_back(did);
        	W[sid].push_back(w);
		M+=w;
		isSelfLoop = true;
	}else{
        	A[sid].push_back(did);
        	A[did].push_back(sid);
        	W[sid].push_back(w);
        	W[did].push_back(w);
		M+=w;
	}
        isUnweighted = isUnweighted & (abs(w-1.0)<=1e-20);	
    }
    edgeNum = wid;
}

void writeLabels(const string filename, const vector<int>& c, const vector<bool>& x, const vector<double> p_values, const double pval, int startIndex)
{
    FILE* fid = fopen(filename.c_str(), "w");
    fprintf(fid, "id%cc%cx%csig\n", delimiter, delimiter, delimiter);
    int csize = c.size();
    for (int i = 0; i < csize; i++) {
        fprintf(fid, "%d%c%d%c%d%c%d\n", i + startIndex, delimiter, c[i] + startIndex, delimiter, !!(x[i]), delimiter, !!(p_values[c[i]] <= pval));
    }
    fclose(fid);
}
