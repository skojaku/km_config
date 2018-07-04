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

#include <graph.h>

#if !defined(MAX)
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#endif
using namespace std;

class KMAlgorithm{
public:	
	// Constructor 
	KMAlgorithm();
	//KMAlgorithm(int num_runs, double significance_level);
	
	// Getter
	vector<int> get_c () const;
	vector<bool> get_x () const;
	vector<double> get_p_values () const;
	
	// Setter
	void set_significance_level (double s);
	void set_num_rand_nets (double r);
	
	// Detect significant CP pairs in networks 
	void detect(const Graph& G);
	
protected:
	vector<int> _c; // _c[i] indicates the index of the CP pair of node i 
	vector<bool> _x;  // x[i]=1 or x[i]=0 indicates a core or a periphery, respectively.  
	vector<double> _p_values; // p_values 
	double _Q; // quality value 
	vector<double> _q;  // quality values
	
	double _significance_level; // statistical significance level
	int _num_rand_nets; // number of randomised networks to be generated
	int _num_runs; // Number of runs of the algorithm 
	mt19937_64 _mtrnd; // random number generator

	/* --------------------------
	 Functions needed to be implemented	
	-------------------------- */
	
	// Detect CP structure and compute their quality
	virtual void _detect_( //
	    const Graph& G,
	    vector<int>& c,
	    vector<bool>& x,
	    double& Q,
	    vector<double>& q,
            mt19937_64& mtrnd) = 0;
	
	// Initialise parameter of randomised networks generator 
	virtual void _init_randomised_network_generator(const Graph& G)= 0;	
	
	// Generate randomised networks 
        virtual void _generate_randomised_network(Graph& G, mt19937_64& mtrnd) = 0;
	
	// Compute the quality of CP pairs 
	virtual void _calc_Q(
	    const Graph& G,
	    const vector<int>& c,
	    const vector<bool>& x,
	    double& Q,
	    vector<double>& q) = 0;
	
	/* --------------------------
	 Statistical test	
	-------------------------- */
	void _estimate_statistical_significance(
	    const Graph& G,
	    const vector<int>& c,
	    const vector<bool>& x,
	    const int num_of_rand_nets,
    	    vector<double>& p_values
	    );


	/* --------------------------
	Other utility functions 
	-------------------------- */
	double _normcdf(double value);

	mt19937_64 _init_random_number_generator();
};

KMAlgorithm::KMAlgorithm(){
	_mtrnd = _init_random_number_generator(); 
	_num_rand_nets = 500;
	_num_runs = 10;
}

mt19937_64 KMAlgorithm::_init_random_number_generator(){
	mt19937_64 mtrnd;
	random_device r;
	seed_seq seed{ r(), r(), r(), r(), r(), r(), r(), r() };
	mtrnd.seed(seed);
	return mtrnd;
}

// Getter
vector<int> KMAlgorithm::get_c() const{
	return _c;
}
vector<bool> KMAlgorithm::get_x() const{
	return _x;
}

vector<double> KMAlgorithm::get_p_values() const{
	return _p_values;
}

void KMAlgorithm::set_significance_level (double s){
	_significance_level = s;
}

void KMAlgorithm::set_num_rand_nets (double r){
	_num_rand_nets = r;
}

void KMAlgorithm::detect(const Graph& G){
	_detect_(G, _c, _x, _Q, _q, _mtrnd);
	
	/* Statistical test */
	int K = _q.size();
	vector<double> tmp(K,0.0);
	_p_values = tmp;
	if (_significance_level < 1.0) {
	    _estimate_statistical_significance(G, _c, _x, _num_rand_nets, _p_values);
	}
}

void KMAlgorithm::_estimate_statistical_significance(
    const Graph& G,
    const vector<int>& c,
    const vector<bool>& x,
    const int num_of_rand_nets,
    vector<double>& p_values
)
{

    /* Initialise variables */
    bool noSelfloop = false;
    bool isunweighted = false;
    int K = *max_element(c.begin(), c.end()) + 1;
    int N = G.get_num_nodes();
 
    double Q;
    vector<double> q;
    vector<int> n(K);
    vector<double> deg(N);
    fill(n.begin(), n.end(), 0);
    _calc_Q(G, c, x, Q, q);
    for (int i = 0; i < N; i++) {
        n[c[i]]++;
    };
	
    /* Generate \hat q^{(s)} and \hat n^{(s)} (1 \leq s \leq S) */
    int numthread;// create random number generator per each thread
    _init_randomised_network_generator(G);

    # pragma omp parallel
    {
    	numthread = omp_get_num_threads();
    }
    vector<mt19937_64> mtrnd_list(numthread);
    for(int i = 0; i < numthread; i++){
	mtrnd_list[i] = _init_random_number_generator();
    }
    
    vector<int> nhat;
    vector<double> qhat;
    #ifdef _OPENMP
    #pragma omp parallel for shared(nhat, qhat, mtrnd_list)
    #endif
    for (int it = 0; it < num_of_rand_nets; it++) {

        // Generate a randomised network using the configuration model.
        Graph G_rand(N);
        int tid = omp_get_thread_num();
        mt19937_64 mtrnd = mtrnd_list[tid];
		
	_generate_randomised_network(G_rand, mtrnd);
	//_Chung_Lu_Algorithm(deg, deg_rank, G_rand, noSelfloop, isunweighted, mtrnd);

        // Detect core-periphery pairs using the KM--config algorithm
        vector<int> c_rand;
        vector<bool> x_rand;
        vector<double> q_rand;
        double Q_rand;
        _detect_(G_rand, c_rand, x_rand, Q_rand, q_rand, mtrnd);
	
        // Save the quality and size of core-periphery pairs in the randomised network.
        int K_rand = q_rand.size();
        vector<int> nsr(K_rand, 0);
        for (int i = 0; i < N; i++) {
            nsr[c_rand[i]]++;
        }

        #ifdef _OPENMP
        #pragma omp critical
        #endif
        for (int k = 0; k < K_rand; k++) {
            nhat.push_back(nsr[k]);
            qhat.push_back(q_rand[k]);
        }
    }

    /* Compute the mean and variance of the quality and size */
    int S = nhat.size();
    double mu_n = (double)accumulate(nhat.begin(), nhat.end(), 0.0) / (double)S;
    double mu_q = (double)accumulate(qhat.begin(), qhat.end(), 0.0) / (double)S;
    double sig_nn = 0;
    double sig_qq = 0;
    double sig_nq = 0;
    for (int s = 0; s < S; s++) {
        sig_nn += pow((double)nhat[s] - mu_n, 2) / (double)(S - 1);
        sig_qq += pow(qhat[s] - mu_q, 2) / (double)(S - 1);
        sig_nq += ((double)nhat[s] - mu_n) * (qhat[s] - mu_q) / (double)(S - 1);
    }

    /* Compute p-values using the Gaussian kernel density estimator */
    double h = MAX(pow((double)S, -1.0 / 6.0), 1e-32);
    p_values.clear();
    p_values.assign(K, 1.0);
    for (int k = 0; k < K; k++) {
        double numer = 0.0;
        double denom = 0.0;
        for (int s = 0; s < S; s++) {
            double qbar = qhat[s] +  sig_nq / sig_nn * (double)(n[k] - nhat[s]);

            double t = sig_nn * (q[k] - qbar) / (sqrt(sig_nn * sig_qq - sig_nq * sig_nq) * h);
            double cum = _normcdf(t);
            double w = exp(- (double)pow(n[k] - nhat[s], 2) / (2.0 * h * h * sig_nn)) + 1e-33;

            numer += cum * w;
            denom += w;
        }
        p_values[k] = 1.0 - numer / denom;
    }
}

double KMAlgorithm::_normcdf(double value)
{
    return 0.5 + 0.5 * erf(value * M_SQRT1_2);
}
