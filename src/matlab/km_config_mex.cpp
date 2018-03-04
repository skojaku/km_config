/*
*
* MATLAB wrapper for the C++ codes of the KM-config algorithm (km_config.h and km_config.cpp)
*
*
* "Core-periphery structure requires something else in the network"
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
* 
*
* COMPILE:
* 
*   mex CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3'  GCC='g++' src/km_config_mex.cpp
* 
*/

#include "../lib/km_config.h"
#include "../lib/km_config.cpp"
#include "mex.h"

void mexFunction(int nlhs, mxArray* plhs[],
    int nrhs, const mxArray* prhs[])
{
    double* edgeList = mxGetPr(prhs[0]);
    int N = (int)mxGetPr(prhs[1])[0];
    int Enum = (int)mxGetPr(prhs[2])[0];
    int num_of_runs = (int)mxGetPr(prhs[3])[0];
    int num_of_rand_nets = (int)mxGetPr(prhs[4])[0];

    /* Parse Input */
    vector<vector<int>> A(N, vector<int>(0));
    vector<vector<double>> W(N, vector<double>(0));

    for (int i = 0; i < Enum; i++) {
        int rid = (edgeList[i] - 1);
        int cid = (edgeList[i + (int)Enum] - 1);
        double w = edgeList[i + 2*(int)Enum];
        A[rid].push_back(cid);
        W[rid].push_back(w);

        if (rid != cid) {
            A[cid].push_back(rid);
            W[cid].push_back(w);
        }
    }

    vector<int> c(N);
    vector<bool> x(N);
    double Q;
    vector<double> Qs;
    
    mt19937_64 mtrnd;
    random_device r;
    seed_seq seed{ r(), r(), r(), r(), r(), r(), r(), r() };
    mtrnd.seed(seed);
    km_config_label_switching(A, W, num_of_runs, c, x, Q, Qs, mtrnd);

    
    int K = Qs.size();
    vector<double> p_values;
    if(num_of_rand_nets>=1){
    	estimate_statistical_significance(A, W, c, x, num_of_runs, num_of_rand_nets, p_values);
    }else{
	p_values.assign(K,1.0);
    }

    plhs[0] = mxCreateDoubleMatrix((mwSize)N, (mwSize)1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix((mwSize)N, (mwSize)1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix((mwSize)1, (mwSize)1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix((mwSize)K, (mwSize)1, mxREAL);
    plhs[4] = mxCreateDoubleMatrix((mwSize)K, (mwSize)1, mxREAL);

    double* retc = mxGetPr(plhs[0]);
    double* retx = mxGetPr(plhs[1]);
    double* retQ = mxGetPr(plhs[2]);
    double* retQs = mxGetPr(plhs[3]);
    double* retPvals = mxGetPr(plhs[4]);

    retQ[0] = Q;
    for (int i = 0; i < N; i++) {
        retc[i] = c[i] + 1;
        retx[i] = !!(x[i]);
    }
    for (int k = 0; k < K; k++) {
        retQs[k] = Qs[k];
        retPvals[k] = p_values[k];
    }

}
