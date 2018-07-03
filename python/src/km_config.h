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

class KM_config: public KMAlgorithm{
public:
	// Constructor 
	KM_config();
	KM_config(int num_runs, double significance_level);
	
protected: // function needed to be implemented

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
	int _num_runs;
	uniform_real_distribution<double> _udist;
	
	// variables for statistical test
	vector<double> _deg;
	vector<int> _deg_rank;
        bool _noSelfloop;
        bool _isunweighted;
        int _N;

	void _km_config_label_switching(
	    const Graph& G,
	    const int num_of_runs,
	    vector<int>& c,
	    vector<bool>& x,
	    double& Q,
	    vector<double>& q,
            mt19937_64& mtrnd
		);

	void _Chung_Lu_Algorithm(const vector<double>& deg, const vector<int>& nodes, Graph& G, bool noSelfloop, bool isunweighted, mt19937_64& mtrnd);

	double _calc_dQ_conf(double d_i_c,
	    double d_i_p,
	    double d_i,
	    double D_c,
	    double D_p,
	    double selfloop,
	    bool x,
	    const double M);
	
	void _propose_new_label(
	    const Graph& G,
	    const vector<int>& c,
	    const vector<bool>& x,
	    const vector<double>& sum_of_deg_core,
	    const vector<double>& sum_of_deg_peri,
	    const double M,
	    const int node_id,
	    int& cprime,
	    bool& xprime,
	    double& dQ,
	    mt19937_64& mtrnd
	    );
	
	
	void _km_config_label_switching_core(
	    const Graph& G,
	    vector<int>& c,
	    vector<bool>& x,
	    mt19937_64& mtrnd
	    );
};


/*-----------------------------
Constructor
-----------------------------*/
KM_config::KM_config(int num_runs, double significance_level):KMAlgorithm(){
	KM_config();
	_num_runs = num_runs;
	_significance_level = significance_level;
};

KM_config::KM_config():KMAlgorithm(){
	_num_runs = 10;
	_significance_level = 0.05;
	_num_rand_nets = 500;
	//_null_model = "config"; 
	//_algorithm = "louvain"; 
	
	uniform_real_distribution<double> tmp(0.0, 1.0);
	_udist = tmp; 
	_mtrnd = _init_random_number_generator();
};


/*-----------------------------
Functions inherited from the super class (KMAlgorithm)
-----------------------------*/
void KM_config::_detect_(
	    const Graph& G,
	    vector<int>& c,
	    vector<bool>& x,
	    double& Q,
	    vector<double>& q,
            mt19937_64& mtrnd){

	_km_config_label_switching(G, _num_runs, c, x, Q, q, mtrnd);
}

void KM_config::_init_randomised_network_generator(const Graph& G){

    /* Initialise variables */
    _noSelfloop = false;
    _isunweighted = false;
    _N = G.get_num_nodes();
 
    vector<double> tmp(_N, 0.0);
    _deg = tmp;
    for (int i = 0; i < _N; i++) {
        _deg[i]= G.wdegree(i);
    };

    vector<int> tmp2(_N); // deg_rank[k] is the id of the node with the kth largest degree. 
    iota(tmp2.begin(), tmp2.end(), 0);
    _deg_rank = tmp2;
    sort(
        _deg_rank.begin(),
        _deg_rank.end(),
        [&](int x, int y){return _deg[x] > _deg[y];}
    );
}

void KM_config::_generate_randomised_network(Graph& G_rand, mt19937_64& mtrnd){
	_Chung_Lu_Algorithm(_deg, _deg_rank, G_rand, _noSelfloop, _isunweighted, mtrnd);
}

void KM_config::_calc_Q(
    const Graph& G,
    const vector<int>& c,
    const vector<bool>& x,
    double& Q,
    vector<double>& q)
{
    int N = G.get_num_nodes();
    int K = *max_element(c.begin(), c.end()) + 1;
    q.assign(K, 0.0);
    vector<double> Dc(K);
    vector<double> Dp(K);
    fill(Dc.begin(), Dc.end(), 0.0);
    fill(Dp.begin(), Dp.end(), 0.0);

    double double_M = 0.0;
    for (int i = 0; i < N; i++) {
	int sz = G.degree(i);
	double di = 0;
        for (int j = 0; j < sz; j++) {
	    int nei = -1;
	    double wj = -1;
	    G.get_weight(i, j, nei, wj);
            q[c[i]] += wj * !!(c[i] == c[nei]) * !!(x[i] | x[nei]);
	    di+=wj;
        }
        Dc[c[i]] += !!(x[i]) * di;
        Dp[c[i]] += !!(!x[i]) * di;
        double_M += di;
    }
    Q = 0;
    for (int k = 0; k < K; k++) {
        q[k] = (q[k] - (Dc[k] * Dc[k] + 2 * Dc[k] * Dp[k]) / double_M) / double_M;
        Q += q[k];
    }
}

/*-----------------------------
Private functions (internal use only)
-----------------------------*/

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
void KM_config::_Chung_Lu_Algorithm(const vector<double>& deg, const vector<int>& nodes, Graph& G, bool noSelfloop, bool isunweighted, mt19937_64& mtrnd)
{
    int N = deg.size();
    double M = accumulate(deg.begin(), deg.end(), 0.0);
    M /= 2;
    Graph tmp(N);
    G = tmp;
	
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
			double r = _udist(mtrnd);
			v = v + floor( log(r) / log(1-p) );
		}
		if(v < N ){
			double q = MIN(deg[ nodes[u] ] * deg[ nodes[v] ] / (2.0*M), 1);
			double w = 1;
			bool addEdge = false;
			if(isunweighted){
				double r = _udist(mtrnd);
				addEdge = r < q / p;	
			}else{
	    			poisson_distribution<int> distribution(q / p);
	    			w = distribution(mtrnd);
				addEdge = w>0; 	
			}
			if(addEdge){
				G.addEdge(nodes[u], nodes[v], w);
			
				if(u!=v){
					G.addEdge(nodes[v], nodes[u], w);
				}
			}
			p = q;
			v = v + 1;
		}
	}
    }
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


void KM_config::_km_config_label_switching(
    const Graph& G,
    const int num_of_runs,
    vector<int>& c,
    vector<bool>& x,
    double& Q,
    vector<double>& q,
    mt19937_64& mtrnd
    )
{

    Q = -1;
    for (int i = 0; i < num_of_runs; i++) {
        vector<int> ci;
        vector<bool> xi;
        vector<double> qi;
        double Qi = 0.0;

        _km_config_label_switching_core(G, ci, xi, mtrnd);

        _calc_Q(G, ci, xi, Qi, qi);
	
        if (Qi > Q) {
            c = ci;
            x = xi;
            Q = Qi;
            q = qi;
        }
    }
}


double KM_config::_calc_dQ_conf(double d_i_c,
    double d_i_p,
    double d_i,
    double D_c,
    double D_p,
    double selfloop,
    bool x,
    const double M)
{
    return 2 * (d_i_c + d_i_p * (!!(x)) - d_i * (D_c + D_p * !!(x)) / (2.0 * M)) + !!(x) * (selfloop - d_i * d_i / (2.0 * M));
}


void KM_config::_propose_new_label(
    const Graph& G,
    const vector<int>& c,
    const vector<bool>& x,
    const vector<double>& sum_of_deg_core,
    const vector<double>& sum_of_deg_peri,
    const double M,
    const int node_id,
    int& cprime,
    bool& xprime,
    double& dQ,
    mt19937_64& mtrnd)
{
    int N = G.get_num_nodes();
    int neighbourNum = G.degree(node_id);

    double deg = G.wdegree(node_id);

    vector<double> edges_to_core(N);
    vector<double> edges_to_peri(N);

    fill(edges_to_core.begin(), edges_to_core.end(), 0.0);
    fill(edges_to_peri.begin(), edges_to_peri.end(), 0.0);
    double selfloop = 0;
    for (int j = 0; j < neighbourNum; j++) {
        int nei = -1;
        double wj = -1;
        G.get_weight(node_id, j, nei, wj);
	
	if(node_id == nei){
		selfloop+= 1;
		continue;
	}
	
        edges_to_core[c[nei]] += wj * (double)!!(x[nei]);
        edges_to_peri[c[nei]] += wj * (double)!!(!x[nei]);
    }

    double D_core = sum_of_deg_core[c[node_id]] - deg * (double)!!(x[node_id]);
    double D_peri = sum_of_deg_peri[c[node_id]] - deg * (double)!!(!x[node_id]);
    double dQold = _calc_dQ_conf(edges_to_core[c[node_id]], edges_to_peri[c[node_id]], deg,
        D_core, D_peri, selfloop, x[node_id], M);

    dQ = 0;
    for (int j = 0; j < neighbourNum; j++) {
        int nei = -1;
        double wj = -1;
        G.get_weight(node_id, j, nei, wj);
        int cid = c[nei];

        D_core = sum_of_deg_core[cid] - deg * (double)!!( (c[node_id] == cid) & x[node_id]);
        D_peri = sum_of_deg_peri[cid] - deg * (double)!!( (c[node_id] == cid) & !x[node_id]);

        double Q_i_core = _calc_dQ_conf(edges_to_core[cid], edges_to_peri[cid],
            deg, D_core, D_peri, selfloop, true, M);
        double Q_i_peri = _calc_dQ_conf(edges_to_core[cid], edges_to_peri[cid],
            deg, D_core, D_peri, selfloop, false, M);
        Q_i_core -= dQold;
        Q_i_peri -= dQold;

        if (MAX(Q_i_core, Q_i_peri) < dQ)
            continue;

        if (Q_i_peri < Q_i_core) {
            xprime = true;
            cprime = cid;
            dQ = Q_i_core;
        }
        else if (Q_i_peri > Q_i_core) {
            xprime = false;
            cprime = cid;
            dQ = Q_i_peri;
        }
        else {
            cprime = cid;
            xprime = _udist(mtrnd) < 0.5;
            dQ = Q_i_core;
        }
    }
}




void KM_config::_km_config_label_switching_core(
    const Graph& G,
    vector<int>& c,
    vector<bool>& x,
    mt19937_64& mtrnd
    )
{
    /* Variable declarations */
    int N = G.get_num_nodes();
    vector<double> sum_of_deg_core(N);
    vector<double> sum_of_deg_peri(N);
    vector<int> order(N);
    vector<double> degs(N);
    double M = 0;
    bool isupdated = false;
    fill(sum_of_deg_core.begin(), sum_of_deg_core.end(), 0.0);
    fill(sum_of_deg_peri.begin(), sum_of_deg_peri.end(), 0.0);
    c.clear();
    x.clear();
    c.assign(N, 0);
    x.assign(N, true);
    for (int i = 0; i < N; i++) {
        order[i] = i;
        c[i] = i;
        double deg = G.wdegree(i);
	degs[i] = deg;
        sum_of_deg_core[i] += (double)!!(x[i]) * deg;
        M += deg;
    };
    M = M / 2;

    /* Label switching algorithm */
    do {
        isupdated = false;
        shuffle(order.begin(), order.end(), mtrnd);

        for (int scan_count = 0; scan_count < N; scan_count++) {
            int i = order[scan_count];

            int cprime = c[i]; // c'
            bool xprime = x[i]; // x'

            double dQ = 0;
            _propose_new_label(G, c, x, sum_of_deg_core, sum_of_deg_peri,
                M, i, cprime, xprime, dQ, mtrnd);

            if (dQ <= 0)
                continue;
		
            if ( (c[i] == cprime) & (x[i] == xprime) )
                continue;

            double deg = degs[i];
            if (x[i]) {
                sum_of_deg_core[c[i]] -= deg;
            }
            else {
                sum_of_deg_peri[c[i]] -= deg;
            }

            if (xprime) {
                sum_of_deg_core[cprime] += deg;
            }
            else {
                sum_of_deg_peri[cprime] += deg;
            }

            c[i] = cprime;
            x[i] = xprime;

            isupdated = true;
        }

    } while (isupdated == true);

    /* Remove empty core-periphery pairs */
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
