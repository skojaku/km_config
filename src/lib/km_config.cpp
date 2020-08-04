/*
*
* Implementation codes for the KM-config algorithm (C++ version)
*
*
* An algorithm for finding multiple core-periphery pairs in networks
*
*
* Core-periphery structure requires something else in the network
* Sadamori Kojaku and Naoki Masuda
* New J. Physics, 20, 1367-2630, 2018 
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

/* Function declarations */

// Private functions for km_config_switching
double calc_dQ(double d_i_c,
    double d_i_p,
    double d_i,
    double D_c,
    double D_p,
    double selfloop,
    bool x,
    const double M);


void propose_new_label(
    const vector<vector<int>>& A,
    const vector<vector<double>>& W,
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


void km_config_label_switching_core(
    const vector<vector<int>>& A,
    const vector<vector<double>>& W,
    vector<int>& c,
    vector<bool>& x,
    mt19937_64& mtrnd
    );


// Private functions for estimate_statistical_significance
double normcdf(double value);

/* Implementation codes */
void km_config_label_switching(
    const vector<vector<int>>& A,
    const vector<vector<double>>& W,
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

        km_config_label_switching_core(A, W, ci, xi, mtrnd);

        calc_Qconf(A, W, ci, xi, Qi, qi);

        if (Qi > Q) {
            c = ci;
            x = xi;
            Q = Qi;
            q = qi;
        }
    }
}


void estimate_statistical_significance(
    const vector<vector<int>>& A,
    const vector<vector<double>>& W,
    const vector<int>& c,
    const vector<bool>& x,
    const int num_of_runs,
    const int num_of_rand_nets,
    vector<double>& p_values)
{

    /* Initialise variables */
    bool noSelfloop = false;
    bool isunweighted = false;
    int K = *max_element(c.begin(), c.end()) + 1;
    int N = A.size();
    double Q;
    vector<double> q;
    vector<int> n(K);
    vector<double> deg(N);
    fill(deg.begin(), deg.end(), 0.0);
    fill(n.begin(), n.end(), 0);
    calc_Qconf(A, W, c, x, Q, q);
    for (int i = 0; i < N; i++) {
        deg[i] = accumulate(W[i].begin(), W[i].end(), 0.0);
        n[c[i]]++;
    };
    vector<int> deg_rank(N); // deg_rank[k] is the id of the node with the kth largest degree. 
    iota(deg_rank.begin(), deg_rank.end(), 0);
    sort(
        deg_rank.begin(),
        deg_rank.end(),
        [&](int x, int y){return deg[x] > deg[y];}
    );


    /* Generate \hat q^{(s)} and \hat n^{(s)} (1 \leq s \leq S) */
    int numthread;// create random number generator per each thread
    # pragma omp parallel
    {
    	numthread = omp_get_num_threads();
    }
    vector<mt19937_64> mtrnd_list(numthread);
    for(int i = 0; i < numthread; i++){
	mtrnd_list[i] = init_random_number_generator();
    }
 

    vector<int> nhat;
    vector<double> qhat;
    #ifdef _OPENMP
    #pragma omp parallel for shared(nhat, qhat, deg, mtrnd_list)
    #endif
    for (int it = 0; it < num_of_rand_nets; it++) {

        // Generate a randomised network using the configuration model.
        vector<vector<int>> A_rand;
        vector<vector<double>> W_rand;
        int tid = omp_get_thread_num();
        mt19937_64 mtrnd = mtrnd_list[tid];

	Chung_Lu_Algorithm(deg, deg_rank, A_rand, W_rand, noSelfloop, isunweighted, mtrnd);

        // Detect core-periphery pairs using the KM--config algorithm
        vector<int> c_rand;
        vector<bool> x_rand;
        vector<double> q_rand;
        double Q_rand;
        km_config_label_switching(A_rand, W_rand, num_of_runs, c_rand, x_rand, Q_rand, q_rand, mtrnd);

        // Save the quality and size of core-periphery pairs in the randomised network.
        int K_rand = q_rand.size();
        vector<int> nsr(K_rand);
        fill(nsr.begin(), nsr.end(), 0);
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
            double cum = normcdf(t);
            double w = exp(- (double)pow(n[k] - nhat[s], 2) / (2.0 * h * h * sig_nn)) + 1e-33;

            numer += cum * w;
            denom += w;
        }
        p_values[k] = 1.0 - numer / denom;
    }
}


void calc_Qconf(
    const vector<vector<int>>& A,
    const vector<vector<double>>& W,
    const vector<int>& c,
    const vector<bool>& x,
    double& Q,
    vector<double>& q)
{
    int N = A.size();
    int K = *max_element(c.begin(), c.end()) + 1;
    q.assign(K, 0.0);
    vector<double> Dc(K);
    vector<double> Dp(K);
    fill(Dc.begin(), Dc.end(), 0.0);
    fill(Dp.begin(), Dp.end(), 0.0);

    double double_M = 0.0;
    for (int i = 0; i < N; i++) {
	int Asize = A[i].size();
        for (int j = 0; j < Asize; j++) {
            int nei = A[i][j];
            q[c[i]] += W[i][j] * !!(c[i] == c[nei]) * !!(x[i] | x[nei]);
        }
        double di = accumulate(W[i].begin(), W[i].end(), 0.0);
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


double calc_dQ(double d_i_c,
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


void propose_new_label(
    const vector<vector<int>>& A,
    const vector<vector<double>>& W,
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
    int N = A.size();
    int neighbourNum = A[node_id].size();
    double deg = accumulate(W[node_id].begin(), W[node_id].end(), 0.0);
    vector<double> edges_to_core(N);
    vector<double> edges_to_peri(N);

    fill(edges_to_core.begin(), edges_to_core.end(), 0.0);
    fill(edges_to_peri.begin(), edges_to_peri.end(), 0.0);
    double selfloop = 0;
    for (int j = 0; j < neighbourNum; j++) {
        int nei = A[node_id][j];
	
	if(node_id == nei){
		selfloop+= 1;
		continue;
	}
	
        edges_to_core[c[nei]] += W[node_id][j] * (double)!!(x[nei]);
        edges_to_peri[c[nei]] += W[node_id][j] * (double)!!(!x[nei]);
    }

    double D_core = sum_of_deg_core[c[node_id]] - deg * (double)!!(x[node_id]);
    double D_peri = sum_of_deg_peri[c[node_id]] - deg * (double)!!(!x[node_id]);
    double dQold = calc_dQ(edges_to_core[c[node_id]], edges_to_peri[c[node_id]], deg,
        D_core, D_peri, selfloop, x[node_id], M);

    dQ = 0;
    for (int j = 0; j < neighbourNum; j++) {
        int nei = A[node_id][j];
        int cid = c[nei];

        D_core = sum_of_deg_core[cid] - deg * (double)!!( (c[node_id] == cid) & x[node_id]);
        D_peri = sum_of_deg_peri[cid] - deg * (double)!!( (c[node_id] == cid) & !x[node_id]);

        double Q_i_core = calc_dQ(edges_to_core[cid], edges_to_peri[cid],
            deg, D_core, D_peri, selfloop, true, M);
        double Q_i_peri = calc_dQ(edges_to_core[cid], edges_to_peri[cid],
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
            xprime = udist(mtrnd) < 0.5;
            dQ = Q_i_core;
        }
    }
}


double normcdf(double value)
{
    return 0.5 + 0.5 * erf(value * 1.0 / sqrt(2));
}


void km_config_label_switching_core(
    const vector<vector<int>>& A,
    const vector<vector<double>>& W,
    vector<int>& c,
    vector<bool>& x,
    mt19937_64& mtrnd
    )
{
    /* Variable declarations */
    int N = A.size();
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
        double deg = accumulate(W[i].begin(), W[i].end(), 0.0);
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
            propose_new_label(A, W, c, x, sum_of_deg_core, sum_of_deg_peri,
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

/*
void generate_randomised_net(
    vector<double>& deg,
    vector<vector<int>>& A,
    vector<vector<double>>& W)
{
    int N = deg.size();
    double M = accumulate(deg.begin(), deg.end(), 0.0);
    M /= 2;
    A.clear();
    W.clear();
    for (int i = 0; i < N; i++) {
        vector<int> tmp;
        A.push_back(tmp);
        vector<double> tmp2;
        W.push_back(tmp2);
    }

    double rnd = 0;
    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            poisson_distribution<int> distribution(deg[i] * deg[j] / (2.0 * (double)M));
            rnd = distribution(mtrnd);
            if (rnd <1 ) continue;
            A[i].push_back(j);
            W[i].push_back((double)rnd);
            
            if(i == j) continue;
                
            A[j].push_back(i);
            W[j].push_back((double)rnd);
        }
    }
}
*/
