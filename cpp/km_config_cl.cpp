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

#include "km_config.h"
#include "km_config.cpp"
#include <string>
#include <stdio.h>
#include <unistd.h>

char delimiter = ' '; // the delimiter for the [input-file] and [output-file]

void split(const string& s, char c,
    vector<string>& v);

void readEdgeTable(string filename, vector<vector<int> >& A, int& N);

void writeLabels(const string filename, const vector<int>& c, const vector<bool>& x, const vector<double> p_values, const double pval);

void usage();

std::string myreplace(std::string &s,
                      const std::string &toReplace,
                      const std::string &replaceWith)
{
    return(s.replace(s.find(toReplace), toReplace.length(), replaceWith));
}

// initialise mtrnd
void init_random_number_generator()
{
    /* If you can initialise mtrnd with random_device, use the following codes.*/
    random_device r;
    seed_seq seed{ r(), r(), r(), r(), r(), r(), r(), r() };
    mtrnd.seed(seed);

    /* Otherwise use the following codes */
    /*
	int seeds[624];
	size_t size = 624 * 4; // Declare size of data
	std::ifstream urandom("/dev/urandom", std::ios::in | std::ios::binary); // Open stream
	if (urandom) // Check if stream is open
	{
	    urandom.read(reinterpret_cast<char*>(seeds), size); // Read from urandom
	    urandom.close();// Close stream
	}
	else // Open failed
	{
	    std::cerr << "Failed to open /dev/urandom" << std::endl;
	}
	std::seed_seq seed(&seeds[0], &seeds[624]);
	mtrnd.seed(seed);
    */
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

    int opt;
    opterr = 0;
    string tmp;
    while ((opt = getopt(argc, argv, "hr:a:l:d:")) != -1) {
        switch (opt) {
        case 'h':
            usage();
            break;
        case 'a':
            tmp.assign(optarg);
            alpha = atof(tmp.c_str());
            break;
        case 'r':
            tmp.assign(optarg);
            num_of_runs = atoi(tmp.c_str());
            break;
        case 'l':
            tmp.assign(optarg);
            num_of_rand_nets = atoi(tmp.c_str());
            break;
        case 'd':
            tmp.assign(optarg);
	    if(!tmp.compare("\\t")){
	    	tmp = myreplace( tmp, "\\t", "\t");
	    	delimiter = *tmp.c_str();
	    }else{
	    	delimiter = *tmp.c_str();
	    }
            break;
        default: /* '?' */
            usage();
            break;
        }
    }
  
   
    /* Read file */
    int N;
    vector<vector<int> > A;
    readEdgeTable(linkfile, A, N);

    /* Run the KM algorithm */
    vector<int> c(N);
    vector<bool> x(N);
    double Q;
    vector<double> q;
    srand(time(NULL));
    init_random_number_generator();
    km_config_label_switching(A, num_of_runs, c, x, Q, q);

    /* Statistical test */
    int K = q.size();
    vector<double> p_values(K);
    fill(p_values.begin(), p_values.end(), 0.0);
    if (alpha < 1.0) {
        estimate_statistical_significance(A, c, x, num_of_runs, num_of_rand_nets, p_values);
    }

    /* Save results */
    double corrected_pval = 1.0 - pow(1.0 - alpha, 1.0 / (double)K);
    writeLabels(outputfile, c, x, p_values, corrected_pval);

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
    cout << "    The node's ID is assumed to start from 1." << endl
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

void readEdgeTable(string filename, vector<vector<int> >& A, int& N)
{

    std::ifstream ifs(filename);
    vector<int> edgeList;
    string str;
    N = 0;
    edgeList.clear();
    while (getline(ifs, str)) {
        vector<string> v;
        split(str, delimiter, v);

        int sid = stoi(v[0]) - 1;
        int did = stoi(v[1]) - 1;

        if (sid == did)
            continue;

        if (N < sid)
            N = sid;
        if (N < did)
            N = did;
        vector<int> tmp;
        edgeList.push_back(sid);
        edgeList.push_back(did);
    }
    N = N + 1;

    vector<int> eids;
    for (int i = 0; i < edgeList.size(); i += 2) {
        int sid = edgeList[i];
        int did = edgeList[i + 1];
        eids.push_back(MIN(sid, did) + N * MAX(sid, did));
    }
    sort(eids.begin(), eids.end());
    eids.erase(unique(eids.begin(), eids.end()), eids.end());

    for (int i = 0; i < N; i++) {
        vector<int> tmp;
        A.push_back(tmp);
    }

    for (int i = 0; i < eids.size(); i++) {
        int sid = eids[i] % N;
        int did = (eids[i] - sid) / N;
        A[sid].push_back(did);
        A[did].push_back(sid);
    }
}

void writeLabels(const string filename, const vector<int>& c, const vector<bool>& x, const vector<double> p_values, const double pval)
{
    FILE* fid = fopen(filename.c_str(), "w");
    fprintf(fid, "id%cc%cx%csig\n", delimiter, delimiter, delimiter);
    int cid = 0;
    for (int i = 0; i < c.size(); i++) {
        fprintf(fid, "%d%c%d%c%d%c%d\n", i + 1, delimiter, c[i] + 1, delimiter, !!(x[i]), delimiter, !!(p_values[c[i]] <= pval));
    }
    fclose(fid);
}
