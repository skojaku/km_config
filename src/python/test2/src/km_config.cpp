#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <algorithm>
#include <iostream>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>

#include <km_config.h>

using namespace std;
namespace py = pybind11;

void readEdgeTable(py::array_t<double> edges_array_t, vector<vector<int>>& A, vector<vector<double>>& W, int& N, int& M)
{

    vector<int> edgeList;
    vector<double> wList;
    N = 0;	
    auto edges = edges_array_t.data();
    auto r = edges_array_t.request();
    M = r.shape[0];

    for(int i =0; i< M; i++){
        int sid = int(edges[3*i]) - 1;
        int did = int(edges[3*i+1]) - 1;
	double w = edges[3*i+2];
	
        if (sid == did)
            continue;

        if (N < sid)
            N = sid;
        if (N < did)
            N = did;
        edgeList.push_back(sid);
        edgeList.push_back(did);
        wList.push_back(w);
    }
    N = N + 1;
   
    vector<vector<int>>tmp(N);
    A = tmp; 
    vector<vector<double>>tmp2(N);
    W = tmp2; 

    int wid = 0; 
    int edgeListsize = edgeList.size();
    for (int i = 0; i < edgeListsize; i += 2) {
        int sid = edgeList[i];
        int did = edgeList[i + 1];
	double w = wList[wid];
	wid++;
        A[sid].push_back(did);
        A[did].push_back(sid);
        W[sid].push_back(w);
        W[did].push_back(w);
    }
}

py::list detect(py::array_t<double> edges, int num_of_runs, double significance_level, int num_of_rand_nets){
        int N = 0;
        int M = 0;	
        vector<vector<int> > A;
        vector<vector<double>> W;
	readEdgeTable(edges, A, W, N, M);
        mt19937_64 mtrnd = init_random_number_generator();
	vector<int> c(N);
	vector<bool> x(N);
	double Q;
	vector<double> q;
	
	km_config_label_switching(A, W, num_of_runs, c, x, Q, q, mtrnd);
	
	/* Statistical test */
	int K = q.size();
	vector<double> p_values(K);
	fill(p_values.begin(), p_values.end(), 0.0);
	double corrected_significance_level =1.0 - pow(1.0 - significance_level, 1.0 / (double)K); // Sidak correction.
	if (significance_level < 1.0) {
	    estimate_statistical_significance(A, W, c, x, num_of_runs, num_of_rand_nets, p_values);
	}
	

	py::array_t<double> cids_array_t(N);
	auto cids = cids_array_t.mutable_data();
	
	py::array_t<double> xs_array_t(N);
	auto xs = xs_array_t.mutable_data();
	
	py::array_t<double> pvals_array_t(N);
	auto pvals = pvals_array_t.mutable_data();
	
	py::array_t<bool> sigs_array_t(N);
	auto sigs = sigs_array_t.mutable_data();	

	for(int i = 0; i < N; i++){
		cids[i] = c[i];
		xs[i] = x[i];
		pvals[i] = p_values[cids[i]];
		sigs[i] = p_values[cids[i]]<=corrected_significance_level;
	}

	py::list results(4);
	results[0] = cids_array_t;
	results[1] = xs_array_t;
	results[2] = pvals_array_t;
	results[3] = sigs_array_t;
	return results;
}

PYBIND11_MODULE(km_config, m){
	m.doc() = "Core-periphery detection in networks";
	m.def("detect", &detect, "Detect multiple core-periphery pairs in networks",
	py::arg("edges"),
	py::arg("num_of_runs") = 10,
	py::arg("significance_level") = 1.0, 
	py::arg("num_of_rand_nets") = 500
	);
}
