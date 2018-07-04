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
#ifndef KM_ALGORITHM 
#define KM_ALGORITHM
	#include "kmalgorithm.h" 
#endif

class KM_modmat: public KMAlgorithm{
public:
	// Constructor 
	KM_modmat();
	KM_modmat(int num_runs, double significance_level);
	
protected: // function needed to be implemented
	int _num_runs; 

	void _init_randomised_network_generator(const Graph& G);

        void _generate_randomised_network(Graph& G, mt19937_64& mtrnd);

	void _detect_(
	    const Graph& G,
	    vector<int>& c,
	    vector<bool>& x,
	    double& Q,
	    vector<double>& q,
            mt19937_64& mtrnd);
	
	void _calc_Q(
	    const Graph& G,
	    const vector<int>& c,
	    const vector<bool>& x,
	    double& Q,
	    vector<double>& q);

private:
	
	void _km_modmat_louvain(
	    const vector<vector<double>>& M,
	    const int num_of_runs,
	    vector<int>& c,
	    vector<bool>& x,
	    double& Q,
	    vector<double>& q,
	    mt19937_64& mtrnd
	    );

	void _km_modmat_louvain_core(
		const vector<vector<double>>& M, 
	    	vector<int>& c,
	    	vector<bool>& x,
	        mt19937_64& mtrnd
		);

	void _km_modmat_label_switching(
	    const vector<vector<double>>& M,
	    const int num_of_runs,
	    vector<int>& c,
	    vector<bool>& x,
	    double& Q,
	    vector<double>& q,
	    mt19937_64& mtrnd
	    );

	void _km_modmat_label_switching_core(
	    const vector<vector<double>>& M,
	    vector<int>& c,
	    vector<bool>& x,
	    mt19937_64& mtrnd
	    );

	void _propose_new_label_modmat(
	    const vector<vector<double>>& M,
	    const vector<int>& c,
	    const vector<bool>& x,
	    const int node_id,
	    int& cprime,
	    bool& xprime,
	    double& dQ,
	    mt19937_64& mtrnd
	    );

	void _calc_Q_modmat(
	    const vector<vector<double>>& M,
	    const vector<int>& c,
	    const vector<bool>& x,
	    double& Q,
	    vector<double>& q);

	void _coarsing(
	    	const vector<vector<double>>& M,
	    	const vector<int>& c,
	    	const vector<bool>& x,
	    	vector<vector<double>>& newM,
	    	vector<int>& toLayerId 
		);

	int _count_non_empty_block(
    		vector<int>& c,
    		vector<bool>& x
		);

	void _relabeling(vector<int>& c);
};



/*-----------------------------
Constructor
-----------------------------*/
KM_modmat::KM_modmat(int num_runs, double significance_level):KMAlgorithm(){
	KM_modmat();
	_num_runs = num_runs;
	_significance_level = significance_level;
};

KM_modmat::KM_modmat(): KMAlgorithm(){
	_num_runs = 10;
	_significance_level = 0.05;
	_num_rand_nets = 500;
};


/*-----------------------------
Functions inherited from the super class (KMAlgorithm)
-----------------------------*/
void KM_modmat::_detect_(
	    const Graph& G,
	    vector<int>& c,
	    vector<bool>& x,
	    double& Q,
	    vector<double>& q,
            mt19937_64& mtrnd){

	vector<vector<double>>M = G.to_matrix();
	
	_km_modmat_louvain(M, _num_runs, c, x, Q, q, mtrnd);
}

// not implemented
void KM_modmat::_init_randomised_network_generator(const Graph& G){

}

// not implemented
void KM_modmat::_generate_randomised_network(Graph& G_rand, mt19937_64& mtrnd){

}

void KM_modmat::_calc_Q(
    const Graph& G,
    const vector<int>& c,
    const vector<bool>& x,
    double& Q,
    vector<double>& q)
{
	vector<vector<double>>M = G.to_matrix();
	_calc_Q_modmat(M,c,x,Q,q);
}

/*-----------------------------
Private functions (internal use only)
-----------------------------*/
void KM_modmat::_km_modmat_label_switching(
    const vector<vector<double>>& M,
    const int num_of_runs,
    vector<int>& c,
    vector<bool>& x,
    double& Q,
    vector<double>& q,
    mt19937_64& mtrnd
    )
{
    /* Generate \hat q^{(s)} and \hat n^{(s)} (1 \leq s \leq S) */
    // create random number generator per each thread
    int numthread = 1;
    #ifdef _OPENMP
    # pragma omp parallel
    {
    	numthread = omp_get_num_threads();
    }
    #endif
    vector<mt19937_64> mtrnd_list(numthread);
    for(int i = 0; i < numthread; i++){
	mt19937_64 mtrnd = _init_random_number_generator();
	mtrnd_list[i] = mtrnd;
    }

    Q = -1;
    int N = M.size();
    #ifdef _OPENMP
    #pragma omp parallel for shared(c, x, Q, q, N, mtrnd_list)
    #endif
    for (int i = 0; i < num_of_runs; i++) {
        vector<int> ci;
        vector<bool> xi;
        vector<double> qi;
        double Qi = 0.0;

        int tid = 0;
    	#ifdef _OPENMP
        	tid = omp_get_thread_num();
    	#endif

        mt19937_64 mtrnd = mtrnd_list[tid];
        _km_modmat_label_switching_core(M, ci, xi, mtrnd);

        _calc_Q_modmat(M, ci, xi, Qi, qi);
	#ifdef _OPENMP
        #pragma omp critical
        #endif
        {
        	if (Qi > Q) {
		    for(int i = 0; i < N; i++){
			c[i] = ci[i];
			x[i] = xi[i];
		    }
		    q.clear();
		    int K = qi.size();
	            vector<double> tmp(K,0.0);	
		    q = tmp;
		    for(int k = 0; k < K; k++){
			q[k] = qi[k];
		    }
        	    Q = Qi;
        	}
	}
    }
}

void KM_modmat::_calc_Q_modmat(
    const vector<vector<double>>& M,
    const vector<int>& c,
    const vector<bool>& x,
    double& Q,
    vector<double>& q)
{
    int N = M.size();
    int K = *max_element(c.begin(), c.end()) + 1;
    
    q.assign(K, 0.0);
    for (int i = 0; i < N; i++) {
    	for (int j = 0; j < N; j++) {
		if( c[i] != c[j]) continue;
		if(x[i] | x[j]){
			q[c[i]]+=M[i][j];
		}	
    	}
    }
    Q = 0;
    for (int k = 0; k < K; k++) {
        Q += q[k];
    }
}


void KM_modmat::_propose_new_label_modmat(
    const vector<vector<double>>& M,
    const vector<int>& c,
    const vector<bool>& x,
    const int node_id,
    int& cprime,
    bool& xprime,
    double& dQ,
    mt19937_64& mtrnd
    )
{
        int N = M.size();
    	int K = *max_element(c.begin(), c.end()) + 1;
	vector<double> dq_core(K,0.0);
	vector<double> dq_peri(K,0.0);
	double dq_old = 0;
	
	for(int i = 0; i < N; i++){
		if(i == node_id) continue;
		dq_core[c[i]]+= M[node_id][i];
		dq_peri[c[i]]+= !!(x[i]) * M[node_id][i];
		dq_old+= !!( (x[i] | x[node_id]) && c[node_id] == c[i]) * M[node_id][i];
	}

	double dqmax = 0;
	dq_old+=!!(x[node_id])*M[node_id][node_id]/2; // add quality induced by self-edges
	for(int k = 0; k < K; k++){
		dq_core[k]+=M[node_id][node_id]/2; // add quality induced by self-edges
		if(dq_core[k] > dq_peri[k]){
			if( dq_core[k]-dq_old >0 && dq_core[k] > dqmax ){
				xprime = true;
				cprime = k;	
				dqmax = dq_core[k];	
				dQ = dq_core[k] - dq_old;	
			}
		}else{
			if( dq_peri[k]-dq_old >0 && dq_peri[k] > dqmax ){
				xprime = false;
				cprime = k;		
				dqmax = dq_peri[k];
				dQ = dq_peri[k] - dq_old;	
			}
		}
	};	
}

void KM_modmat::_relabeling(
    	vector<int>& c
	){

    int N = c.size(); 
    std::vector<int> labs;
    for (int i = 0; i < N; i++) {
        int cid = -1;
	int labsize = labs.size();
        for (int j = 0; j < labsize; j++) {
            if (labs[j] == c[i]) {
                cid = j;
                break;
            }
        }

        if (cid < 0) {
            labs.push_back(c[i]);
            cid = labs.size() - 1;
        }
        c[i] = cid;
    }
}

void KM_modmat::_km_modmat_label_switching_core(
    const vector<vector<double>>& M,
    vector<int>& c,
    vector<bool>& x,
    mt19937_64& mtrnd
    )
{
    /* Variable declarations */
    int N = M.size();
    vector<int> order(N);
    vector<double> degs(N);
    bool isupdated = false;
    c.clear();
    x.clear();
    c.assign(N, 0);
    x.assign(N, true);
    for (int i = 0; i < N; i++) {
        order[i] = i;
        c[i] = i;
    };

    /* Label switching algorithm */
    do {
        isupdated = false;
        shuffle(order.begin(), order.end(), mtrnd);

        for (int scan_count = 0; scan_count < N; scan_count++) {
            int i = order[scan_count];

            int cprime = c[i]; // c'
            bool xprime = x[i]; // x'

            double dQ = 0;
            _propose_new_label_modmat(M, c, x, i, cprime, xprime, dQ, mtrnd);

            if (dQ <= 0)
                continue;

            if ( (c[i] == cprime) & (x[i] == xprime) )
                continue;

            c[i] = cprime;
            x[i] = xprime;

            isupdated = true;
        }

    } while (isupdated == true);

    /* Remove empty core-periphery pairs */
    _relabeling(c);
}


void KM_modmat::_coarsing(
    	const vector<vector<double>>& M,
    	const vector<int>& c,
    	const vector<bool>& x,
    	vector<vector<double>>& newM,
    	vector<int>& toLayerId 
	){
		
        int N = c.size();
	vector<int> ids(N,0);
    	int maxid = 0;
	for(int i = 0;i<N;i++){
		ids[i] = 2 * c[i] + x[i];
		maxid = MAX(maxid, ids[i]);
	}
	_relabeling(ids);
	toLayerId.clear();
	toLayerId.assign(maxid+1,0);
	for(int i = 0;i<N;i++){
		toLayerId[2 * c[i] + x[i]] = ids[i];
	}
	
	
    	int K = *max_element(ids.begin(), ids.end()) + 1;
	vector<vector<double>> tmp(K, vector<double>(K,0));
	newM.clear();
	newM = tmp;

	for(int i = 0;i<N;i++){
		int mi = 2 * c[i] + x[i];
		for(int j = 0;j<N;j++){
			int mj = 2 * c[j] + x[j];
			newM[ toLayerId[mi] ][ toLayerId[mj] ]+=M[i][j];
		}
	}
}

int KM_modmat::_count_non_empty_block(
    	vector<int>& c,
    	vector<bool>& x
	){
	int N = c.size();
	vector<int> ids(N,0);
	for(int i = 0; i< N; i++){
		ids[i] = 2 * c[i] + x[i];
	}
	sort(ids.begin(), ids.end());
	return unique(ids.begin(), ids.end()) - ids.begin();
}

void KM_modmat::_km_modmat_louvain_core(
	const vector<vector<double>>& M, 
    	vector<int>& c,
    	vector<bool>& x,
        mt19937_64& mtrnd
	){

	// Intiialise variables	
	int N = M.size();
	c.clear();
	x.clear();
    	c.assign(N, 0); 
    	x.assign(N, true);
    	for (int i = 0; i < N; i++) c[i] = i;
	
	vector<int>ct = c; // label of each node at tth iteration
	vector<bool>xt = x; // label of each node at tth iteration. 
	vector<vector<double>> cnet_M; // coarse network
	vector<int> toLayerId; //toLayerId[i] maps 2*c[i] + x[i] to the id of node in the coarse network 
	_coarsing(M, ct, xt, cnet_M, toLayerId); // Initialise toLayerId

	double Qbest = 0; // quality of the current partition

	int cnet_N;
	do{
		cnet_N = cnet_M.size();
		
		// Core-periphery detection	
		vector<int> cnet_c; // label of node in the coarse network, Mt 
		vector<bool> cnet_x; // label of node in the coarse network, Mt 
		_km_modmat_label_switching_core(cnet_M, cnet_c, cnet_x, mtrnd);
	
		// Update the label of node in the original network, ct and xt.	
		for(int i = 0; i< N; i++){
			int cnet_id = toLayerId[2 * ct[i] + xt[i]];
			ct[i] = cnet_c[ cnet_id ];
			xt[i] = cnet_x[ cnet_id ];
		}
 		
		// Compute the quality       	
		double Qt = 0; vector<double> qt;
		_calc_Q_modmat(cnet_M, cnet_c, cnet_x, Qt, qt);

		if(Qt>=Qbest){ // if the quality is the highest among those detected so far
			c = ct;
			x = xt;
			Qbest = Qt;
		}
	
		// Coarsing	
		vector<vector<double>> new_cnet_M; 
		_coarsing(cnet_M, cnet_c, cnet_x, new_cnet_M, toLayerId);
		cnet_M = new_cnet_M;
		
		int sz = cnet_M.size();
		if(sz == cnet_N) break;	
			
	}while( true );
	_relabeling(c);
}

void KM_modmat::_km_modmat_louvain(
    const vector<vector<double>>& M,
    const int num_of_runs,
    vector<int>& c,
    vector<bool>& x,
    double& Q,
    vector<double>& q,
    mt19937_64& mtrnd
    )
{

    int N = M.size();
    c.clear();
    x.clear();
    c.assign(N, 0);
    x.assign(N, true);

    /* Generate \hat q^{(s)} and \hat n^{(s)} (1 \leq s \leq S) */
    // create random number generator per each thread
    int numthread = 1;
    #ifdef _OPENMP
    	# pragma omp parallel
    	{
    		numthread = omp_get_num_threads();
    	}
    #endif
    vector<mt19937_64> mtrnd_list(numthread);
    for(int i = 0; i < numthread; i++){
	mt19937_64 mtrnd = _init_random_number_generator();
	mtrnd_list[i] = mtrnd;
    }

    Q = -1;
    #ifdef _OPENMP
    #pragma omp parallel for shared(c, x, Q, q, N, M, mtrnd_list)
    #endif
    for (int i = 0; i < num_of_runs; i++) {
        vector<int> ci;
        vector<bool> xi;
        vector<double> qi;
        double Qi = 0.0;

        int tid = 0;
    	#ifdef _OPENMP
        	tid = omp_get_thread_num();
    	#endif
	
        mt19937_64 mtrnd = mtrnd_list[tid];
        _km_modmat_louvain_core(M, ci, xi, mtrnd);

        _calc_Q_modmat(M, ci, xi, Qi, qi);

        #pragma omp critical
        {
        	if (Qi > Q) {
		    for(int i = 0; i < N; i++){
			c[i] = ci[i];
			x[i] = xi[i];
		    }
		    q.clear();
		    int K = qi.size();
	            vector<double> tmp(K,0.0);	
		    q = tmp;
		    for(int k = 0; k < K; k++){
			q[k] = qi[k];
		    }
        	    Q = Qi;
        	}
	}
    }
}
