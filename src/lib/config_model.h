#include <array>
#include <random>
#include <vector>
#include <fstream>
#include "math.h"
#include <random>
#include <algorithm>
#include <iostream>     // std::cout, std::end
#include <numeric>      // std::iota
#include <omp.h>

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif
using namespace std;

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
	mt19937_64 mtrnd;
    	random_device r;
    	seed_seq seed{ r(), r(), r(), r(), r(), r(), r(), r() };
    	mtrnd.seed(seed);
	mtrnd_list[i] = mtrnd;
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
