#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <algorithm>
#include <iostream>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>

#include <km_config.h>
#include <km_modmat.h>

using namespace std;
namespace py = pybind11;

void readEdgeTable(py::array_t<int> edges_array_t, py::array_t<double> w_array_t, Graph& G)
{

    vector<int> edgeList;
    vector<double> wList;
    int N = 0;	
    auto edges = edges_array_t.data();
    auto r = edges_array_t.request();
    int M = r.shape[0];
    auto ws = w_array_t.data();

    for(int i =0; i< M; i++){
        int sid = edges[2*i];
        int did = edges[2*i + 1];
	double w = ws[i];
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
  
    Graph tmp(N);
    G = tmp; 

    int wid = 0; 
    int edgeListsize = edgeList.size();
    for (int i = 0; i < edgeListsize; i += 2) {
        int sid = edgeList[i];
        int did = edgeList[i + 1];
	double w = wList[wid];
	G.addEdge(sid, did, w);
	wid++;
    }
}

/* KM algorithm based on the configuration model */
py::list detect_config(py::array_t<int> edges, py::array_t<double> ws, int num_of_runs, double significance_level, int num_of_rand_nets){
       	Graph G(0); 
	readEdgeTable(edges, ws, G);
	int N = G.get_num_nodes();

        //mt19937_64 mtrnd = init_random_number_generator();
	vector<double> q;

	KM_config km = KM_config(num_of_runs, significance_level);
	km.detect(G);
	
	vector<int>  c = km.get_c();
	vector<bool> x = km.get_x();
	vector<double> p_values = km.get_p_values();	
	int K = p_values.size();
	
	py::array_t<double> cids_array_t(N);
	auto cids = cids_array_t.mutable_data();
	
	py::array_t<double> xs_array_t(N);
	auto xs = xs_array_t.mutable_data();
	
	py::array_t<double> pvals_array_t(K);
	auto pvals = pvals_array_t.mutable_data();
	
	for(int i = 0; i < N; i++){
		cids[i] = c[i];
		xs[i] = x[i];
	}
	for(int i = 0; i < K; i++){
		pvals[i] = p_values[i];
	}

	py::list results(3);
	results[0] = cids_array_t;
	results[1] = xs_array_t;
	results[2] = pvals_array_t;
	return results;
}

/* KM algorithm based on an arbitary null models*/
py::list detect_modmat(py::array_t<int> edges, py::array_t<double> ws, int num_of_runs, double significance_level, int num_of_rand_nets){
       	Graph G(0); 
	readEdgeTable(edges, ws, G);
	int N = G.get_num_nodes();
	KM_modmat km = KM_modmat(num_of_runs, significance_level);

	km.detect(G);
	
	vector<int>  c = km.get_c();
	vector<bool> x = km.get_x();
	vector<double> p_values = km.get_p_values();	
	
	py::array_t<double> cids_array_t(N);
	auto cids = cids_array_t.mutable_data();
	
	py::array_t<double> xs_array_t(N);
	auto xs = xs_array_t.mutable_data();
	
	py::array_t<double> pvals_array_t(N);
	auto pvals = pvals_array_t.mutable_data();
	
	for(int i = 0; i < N; i++){
		cids[i] = c[i];
		xs[i] = x[i];
		pvals[i] = p_values[cids[i]];
	}

	py::list results(3);
	results[0] = cids_array_t;
	results[1] = xs_array_t;
	results[2] = pvals_array_t;
	return results;
}

PYBIND11_MODULE(_kmalgorithm, m){
	m.doc() = "Core-periphery detection in networks";
	m.def("detect_config", &detect_config, "Use the configuration model as null models",
		py::arg("edges"),
		py::arg("ws"),
		py::arg("num_of_runs") = 10,
		py::arg("significance_level") = 1.0, 
		py::arg("num_of_rand_nets") = 500
	);
	m.def("detect_modmat", &detect_modmat, "Use the user-provided modularity matrix",
	py::arg("edges"),
	py::arg("ws"),
	py::arg("num_of_runs") = 10,
	py::arg("significance_level") = 1.0, 
	py::arg("num_of_rand_nets") = 500
	);
}
