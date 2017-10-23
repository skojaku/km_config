/*
*
* MATLAB wrapper for the C++ codes of the KM-config algorithm (km_config.h and km_config.cpp)
*
*
* "???"
* Sadamori Kojaku and Naoki Masuda
* Preprint arXiv:???
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
*   mex CXXFLAGS='$CXXFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' CXXOPTIMFLAGS='-O3 -DNDEBUG' LDOPTIMFLAGS='-O3' km_config_mex.cpp
* 
*/

#include "../cpp/km_config.h"
#include "../cpp/km_config.cpp"
#include "mex.h"

void init_random_number_generator(){
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
	std::seed_seq seed(&seeds[0], &seeds[624]);
	mtrnd.seed(seed);
}


void mexFunction(int nlhs, mxArray* plhs[],
    int nrhs, const mxArray* prhs[])
{
    double* edgeList = mxGetPr(prhs[0]);
    int N = (int)mxGetPr(prhs[1])[0];
    int Enum = (int)mxGetPr(prhs[2])[0];
    int num_of_runs = (int)mxGetPr(prhs[3])[0];
    int num_of_rand_nets = (int)mxGetPr(prhs[4])[0];

    /* Parse Input */
    vector<vector<int> > A;
    for (int i = 0; i < N; i++) {
        vector<int> tmp;
        A.push_back(tmp);
    }

    for (int i = 0; i < Enum; i++) {
        int rid = (edgeList[i] - 1);
        int cid = (edgeList[i + (int)Enum] - 1);
        A[rid].push_back(cid);

        if (rid != cid) {
            A[cid].push_back(rid);
        }
    }

    vector<int> c(N);
    vector<bool> x(N);
    double Q;
    vector<double> Qs;
    init_random_number_generator();
    km_config_label_switching(A, num_of_runs, c, x, Q, Qs);

    
    int K = Qs.size();
    vector<double> p_values;
    if(num_of_rand_nets>=1){
    	estimate_statistical_significance(A, c, x, num_of_runs, num_of_rand_nets, p_values);
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
