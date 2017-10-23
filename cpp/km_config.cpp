/*
*
* Implementation code for the KM-config algorithm (C++ version)
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

/* Private function */
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

/* Private function */
void propose_new_label(const vector<vector<int> >& A,
    const vector<int>& c,
    const vector<bool>& x,
    const vector<double>& sum_of_deg_core,
    const vector<double>& sum_of_deg_peri,
    const double M,
    const int node_id,
    int& cprime,
    bool& xprime,
    double& dQ)
{
    int N = A.size();
    int neighbourNum = A[node_id].size();
    int deg = A[node_id].size();
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
	
        edges_to_core[c[nei]] += (double)!!(x[nei]);
        edges_to_peri[c[nei]] += (double)!!(!x[nei]);
    }

    double D_core = sum_of_deg_core[c[node_id]] - deg * (double)!!(x[node_id]);
    double D_peri = sum_of_deg_peri[c[node_id]] - deg * (double)!!(!x[node_id]);
    double dQold = calc_dQ(edges_to_core[c[node_id]], edges_to_peri[c[node_id]], A[node_id].size(),
        D_core, D_peri, selfloop, x[node_id], M);

    dQ = 0;
    for (int j = 0; j < neighbourNum; j++) {
        int nei = A[node_id][j];
        int cid = c[nei];

        D_core = sum_of_deg_core[cid] - deg * (double)!!(c[node_id] == cid & x[node_id]);
        D_peri = sum_of_deg_peri[cid] - deg * (double)!!(c[node_id] == cid & !x[node_id]);

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

/* Private function */
void km_config_label_switching_core(const vector<vector<int> >& A,
    vector<int>& c,
    vector<bool>& x)
{
    /* Variable declarations */
    int N = A.size();
    vector<double> sum_of_deg_core(N);
    vector<double> sum_of_deg_peri(N);
    vector<int> order(N);
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
        double deg = A[i].size();
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
            propose_new_label(A, c, x, sum_of_deg_core, sum_of_deg_peri,
                M, i, cprime, xprime, dQ);

            if (dQ <= 0)
                continue;

            if (c[i] == cprime & x[i] == xprime)
                continue;

            double deg = A[i].size();
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
        for (int j = 0; j < labs.size(); j++) {
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


/* Private function */
double normcdf(double value)
{
    return 0.5 + 0.5 * erf(value * M_SQRT1_2);
}

void km_config_label_switching(const vector<vector<int> >& A,
    const int num_of_runs,
    vector<int>& c,
    vector<bool>& x,
    double& Q,
    vector<double>& q)
{

    Q = -1;
    for (int i = 0; i < num_of_runs; i++) {
        vector<int> ci;
        vector<bool> xi;
        vector<double> qi;
        double Qi = 0.0;

        km_config_label_switching_core(A, ci, xi);

        calc_Qconf(A, ci, xi, Qi, qi);

        if (Qi > Q) {
            c = ci;
            x = xi;
            Q = Qi;
            q = qi;
        }
    }
}

/* Private function */
void generate_randomised_net(vector<double>& deg,
    vector<vector<int> >& A)
{
    int N = deg.size();
    int M = accumulate(deg.begin(), deg.end(), 0.0);
    M /= 2;
    A.clear();
    for (int i = 0; i < N; i++) {
        vector<int> tmp;
        A.push_back(tmp);
    }
    for (int i = 0; i < N; i++) {
        if ((udist(mtrnd)) <= deg[i] * deg[i] / (4.0 * (double)M)){
            A[i].push_back(i);
            A[i].push_back(i);
	}
        
	for (int j = i + 1; j < N; j++) {
            if ((udist(mtrnd)) > deg[i] * deg[j] / (2.0 * (double)M))
                continue;
            A[i].push_back(j);
            A[j].push_back(i);
        }
    }
}

void estimate_statistical_significance(const vector<vector<int> >& A,
    const vector<int>& c,
    const vector<bool>& x,
    const int num_of_runs,
    const int num_of_rand_nets,
    vector<double>& p_values)
{

    /* Initialise variables */
    int K = *max_element(c.begin(), c.end()) + 1;
    int N = A.size();
    double Q;
    vector<double> q;
    vector<int> n(K);
    vector<double> deg(N);
    fill(deg.begin(), deg.end(), 0.0);
    fill(n.begin(), n.end(), 0);
    calc_Qconf(A, c, x, Q, q);
    for (int i = 0; i < N; i++) {
        deg[i] = A[i].size();
        n[c[i]]++;
    };

    /* Generate \hat q^{(s)} and \hat n^{(s)} (1 \leq s \leq S) */
    vector<int> nhat;
    vector<double> qhat;
    #ifdef _OPENMP
    #pragma omp parallel for shared(nhat, qhat, deg)
    #endif
    for (int it = 0; it < num_of_rand_nets; it++) {

        // Generate a randomised network using the configuration model.
        vector<vector<int> > A_rand;
        generate_randomised_net(deg, A_rand);

        // Detect core-periphery pairs using the KM--config algorithm
        vector<int> c_rand;
        vector<bool> x_rand;
        vector<double> q_rand;
        double Q_rand;
        km_config_label_switching(A_rand, num_of_runs, c_rand, x_rand, Q_rand, q_rand);

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

void calc_Qconf(const vector<vector<int> >& A,
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
        for (int j = 0; j < A[i].size(); j++) {
            int nei = A[i][j];
            q[c[i]] += !!(c[i] == c[nei]) * !!(x[i] | x[nei]);
        }
        double di = A[i].size();
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
