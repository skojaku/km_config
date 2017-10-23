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
#ifdef _OEPNMP
#include <omp.h>
#endif

#if !defined(MAX)
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#endif
using namespace std;

/* Global variables */
std::mt19937_64 mtrnd;
uniform_real_distribution<double> udist(0.0, 1.0);


/* 
*
* This function initialises mtrnd. 
*
*/
void init_random_number_generator();


/* ---- KM-config algorithm ----
* INPUT:
*
*   A[i][j] - Adjacency matrix. A[i][j] is the node id of the j-th neighbour of node i.
*
*   num_of_runs - number of times we run the KM-config algorithm.
*
*
* OUTPUT
*
*   c - N-dimensional vector. c[i] is the index of the core-periphery pair to which node i belongs.
*  
*   x - N-dimensional vector. If node i is a core node, x[i] = 1. If node i is a periphery node, x[i] = 0.
*  
*   Q - The quality value of the detected core-periphery structure.
*  
*   q - C-dimensional vector. q[i] represents the contribution of the i-th core-periphery pair to Q.
*/
void km_config_label_switching(const vector<vector<int> >& A,
    const int num_of_runs,
    vector<int>& c,
    vector<bool>& x,
    double& Q,
    vector<double>& q);


/* ---- Quaity function Qconf ----
* INPUT:
*
*   A[i][j] - Adjacency matrix. A[i][j] is the node id of the j-th neighbour of node i.
*
*   c - N-dimensional vector. c[i] is the index of the core-periphery pair to which node i belongs.
*  
*   x - N-dimensional vector. If node i is a core node, x[i] = 1. If node i is a periphery node, x[i] = 0.
*
*
* OUTPUT:
*
*   Q - The quality value of the detected core-periphery structure.
*  
*   q - C-dimensional vector. q[i] represents the contribution of the i-th core-periphery pair to Q.
*
*/
void calc_Qconf(const vector<vector<int> >& A,
    const vector<int>& c,
    const vector<bool>& x,
    double& Q,
    vector<double>& q);


/* ---- Statistical significance of core-periphery pairs ----
* INPUT
*
*   A[i][j] - Adjacency matrix. A[i][j] is the node id of the j-th neighbour of node i.
*
*   c - N-dimensional vector. c[i] is the index of the core-periphery pair to which node i belongs.
*  
*   x - N-dimensional vector. If node i is a core node, x[i] = 1. If node i is a periphery node, x[i] = 0.
*
*   num_of_runs - number of times we run the KM-config algorithm.
*
*   num_of_rand_nets - number of randomised networks.
*
*
* OUTPUT
*
*   p_values - C-dimensional vector. p_values[i] represents the p-value of the i-th core-periphery pair. 
*
*/
void estimate_statistical_significance(const vector<vector<int> >& A,
    const vector<int>& c,
    const vector<bool>& x,
    const int num_of_runs,
    const int num_of_rand_nets,
    vector<double>& p_values);
