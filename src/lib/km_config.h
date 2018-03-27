/*
*
* Header file of the KM-config algorithm (C++ version)
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
* DATE - 11 Oct 2017
*/

#include <algorithm>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <numeric>
#include <cmath>
#include <omp.h>

#if !defined(MAX)
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#endif
using namespace std;

/* Global variables */
//std::mt19937_64 mtrnd;
uniform_real_distribution<double> udist(0.0, 1.0);


/* 
*
* This function initialises mtrnd. 
*
*/
std::mt19937_64 init_random_number_generator(){
	mt19937_64 mtrnd;
    	random_device r;
    	seed_seq seed{ r(), r(), r(), r(), r(), r(), r(), r() };
    	mtrnd.seed(seed);
	return mtrnd;
}

/* 
* Use this function instead if you have problem in initialising mtrnd
* This sometimes happens when the version of gcc compiler is different from the one required by mex compiler.
*/ 
/*
std::mt19937_64 init_random_number_generator(){
	int seeds[624];
	size_t size = 624*4; //Declare size of data
	std::ifstream urandom("/dev/urandom", std::ios::in | std::ios::binary); //Open stream
	if (urandom) //Check if stream is open
	{
	    urandom.read(reinterpret_cast<char*>(seeds), size); //Read from urandom
	    urandom.close(); //close stream
	}
	else //Open failed
	{
	    		std::cerr << "Failed to open /dev/urandom" << std::endl;
	}
	std::mt19937_64 mtrnd;
	std::seed_seq seed(&seeds[0], &seeds[624]);
	mtrnd.seed(seed);
	return mtrnd;
}
*/


/* ---- KM-config algorithm ----
* INPUT:
*
*   A[i][j] - Adjacency matrix. A[i][j] is the node id of the j-th neighbour of node i.
*
*   W[i][j] - Weight matrix. W[i][j] is the weight of edges between node i and the jth neighbour.
*
*   num_of_runs - number of times we run the KM-config algorithm.
*
*
* OUTPUT
*
*   c - N-dimensional vector. c[i] is the index of the core-periphery pair to which node i belongs.
*  
*   x - N-dimensional vector. If node i is a core node, x[i] = 1. If node i is a periphery node, x[i] = 0.
*  
*   Q - The quality value of the detected core-periphery structure.
*  
*   q - C-dimensional vector. q[i] represents the contribution of the i-th core-periphery pair to Q.
*/
void km_config_label_switching(
    const vector<vector<int>>& A,
    const vector<vector<double>>& W,
    const int num_of_runs,
    vector<int>& c,
    vector<bool>& x,
    double& Q,
    vector<double>& q);


/* ---- Quaity function Qconf ----
* INPUT:
*
*   A[i][j] - Adjacency matrix. A[i][j] is the node id of the j-th neighbour of node i.
*
*   W[i][j] - Weight matrix. W[i][j] is the weight of edges between node i and the jth neighbour.
*
*   c - N-dimensional vector. c[i] is the index of the core-periphery pair to which node i belongs.
*  
*   x - N-dimensional vector. If node i is a core node, x[i] = 1. If node i is a periphery node, x[i] = 0.
*
*
* OUTPUT:
*
*   Q - The quality value of the detected core-periphery structure.
*  
*   q - C-dimensional vector. q[i] represents the contribution of the i-th core-periphery pair to Q.
*
*/
void calc_Qconf(
    const vector<vector<int>>& A,
    const vector<vector<double>>& W,
    const vector<int>& c,
    const vector<bool>& x,
    double& Q,
    vector<double>& q);


/* ---- Statistical significance of core-periphery pairs ----
* INPUT
*
*   A[i][j] - Adjacency matrix. A[i][j] is the node id of the j-th neighbour of node i.
*
*   W[i][j] - Weight matrix. W[i][j] is the weight of edges between node i and the jth neighbour.
*
*   c - N-dimensional vector. c[i] is the index of the core-periphery pair to which node i belongs.
*  
*   x - N-dimensional vector. If node i is a core node, x[i] = 1. If node i is a periphery node, x[i] = 0.
*
*   num_of_runs - number of times we run the KM-config algorithm.
*
*   num_of_rand_nets - number of randomised networks.
*
*
* OUTPUT
*
*   p_values - C-dimensional vector. p_values[i] represents the p-value of the i-th core-periphery pair. 
*
*/
void estimate_statistical_significance(
    const vector<vector<int> >& A,
    const vector<vector<double> >& W,
    const vector<int>& c,
    const vector<bool>& x,
    const int num_of_runs,
    const int num_of_rand_nets,
    vector<double>& p_values);


/*
 * Chung, F. & Lu, L. 
 * Connected Components in Random Graphs with Given Expected Degree Sequences. 
 * Ann. Comb. 6, 125–145 (2002). 
 *
 * Miller, J. C. & Hagberg, A. 
 * Efficient Generation of Networks with Given Expected Degrees. 
 * in Algorithms and Models for the Web Graph (eds. Frieze, A., Horn, P. & Prałat, P.) 
 * 6732 LNCS, 115–126 (Springer Berlin Heidelberg, 2011). 
 *
 * */
void Chung_Lu_Algorithm(const vector<double>& deg, const vector<int>& nodes, vector<vector<int>>& A, vector<vector<double>>&W, bool noSelfloop, bool isunweighted, mt19937_64& mtrnd)
{
    uniform_real_distribution<double> udist(0.0, 1.0);
    int N = deg.size();
    double M = accumulate(deg.begin(), deg.end(), 0.0);
    M /= 2;
    A.clear();
    vector<vector<int>> tmp(N, vector<int>(0));
    A =tmp;
    W.clear();
    vector<vector<double>> tmp2(N, vector<double>(0));
    W =tmp2;
	
    // sort degree
/*
    vector<int> nodes(N);
    iota(nodes.begin(), nodes.end(), 0);
    sort(
        nodes.begin(),
        nodes.end(),
        [&](int x, int y){return deg[x] > deg[y];}
    );
*/

    for(int u = 0; u < N-1; u++){
	int v = u;

	if(noSelfloop){
		v = v + 1;
	}	
	
	double p = MIN(1, deg[ nodes[u] ] * deg[ nodes[v] ] / (2.0*M) );
	while(v < N && p > 0){
		if(p!=1){
    			//geometric_distribution<int> gdist(p);
			//v = v + gdist(mtrnd);
			double r = udist(mtrnd);
			v = v + floor( log(r) / log(1-p) );
		}
		if(v < N ){
			double q = MIN(deg[ nodes[u] ] * deg[ nodes[v] ] / (2.0*M), 1);
			double w = 1;
			bool addEdge = false;
			if(isunweighted){
				double r = udist(mtrnd);
				addEdge = r < q / p;	
			}else{
	    			poisson_distribution<int> distribution(q / p);
	    			w = distribution(mtrnd);
				addEdge = w>0; 	
			}
			if(addEdge){
				if(u!=v){
	            			A[nodes[u]].push_back(nodes[v]);
	            			A[nodes[v]].push_back(nodes[u]);
	            			W[nodes[u]].push_back(w);
	            			W[nodes[v]].push_back(w);
				}else{
	            			A[nodes[u]].push_back(nodes[u]);
	            			W[nodes[u]].push_back(2*w);
				}
			}
			p = q;
			v = v + 1;
		}
	}
    }
}



void generate_randomised_nets(
    const vector<double>& deg,
    const int num_of_rand_nets,
    vector<vector<vector<int>>>& Alist,
    vector<vector<vector<double>>>& Wlist,
    bool noSelfLoop,
    bool isunweighted
){
 
    Alist.clear();
    int N = deg.size();
    
    vector<int> deg_rank(N); // deg_rank[k] is the id of the node with the kth largest degree. 
    iota(deg_rank.begin(), deg_rank.end(), 0);
    sort(
        deg_rank.begin(),
        deg_rank.end(),
        [&](int x, int y){return deg[x] > deg[y];}
    );

    // create random number generator per each thread
    int numthread;
    # pragma omp parallel
    {
    numthread = omp_get_num_threads();
    }
    vector<mt19937_64> mtrnd_list(numthread);
    for(int i = 0; i < numthread; i++){
	mtrnd_list[i] = init_random_number_generator();
    }
   
    #ifdef _OPENMP
    #pragma omp parallel for shared(Alist, Wlist, mtrnd_list)
    #endif
    for (int it = 0; it < num_of_rand_nets; it++) {
	
	vector<vector<int>> A_rand;
	vector<vector<double>> W_rand;
       
        int tid = omp_get_thread_num();
        mt19937_64 mtrnd = mtrnd_list[tid];

	Chung_Lu_Algorithm(deg, deg_rank, A_rand, W_rand, noSelfLoop, isunweighted, mtrnd);
	
        #ifdef _OPENMP
        #pragma omp critical
        #endif
	{
		Alist.push_back(A_rand);
		Wlist.push_back(W_rand);
	}
    }
}
