/*
*
* Header file of the Kojaku-Masuda algorithm (C++ version)
*
*
* An algorithm for finding multiple core-periphery pairs in networks
* "???"
* Sadamori Kojaku and Naoki Masuda
* Preprint arXiv:???
* 
*
* Please do not distribute without contacting the authors above.
* If you find a bug in this code, please contact the authors.
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
* This function is for initialising mtrnd. 
* Implement this function for your environment in your .cpp file. 
*
* */
void init_random_number_generator();

/* ---- KM algorithm ----
* INPUT:
*
*  A[i][j] - Adjacency matrix. A[i][j] is the node id of the jth neighbour of node i.
*
*  num_of_runs - number of runs of the KM algorithm.
*
*
* OUTPUT
*
*  c - Nx1 vector. c[i] is the index of core-periphery pair to which each node i belongs, and
*      N is the number of nodes. 
*
*  x - Nx1 vector. x[i] = true if node i is a core node and x[i] = false if it is a peripheral node.
*
*  Q - Quality value of core-periphery pairs.  
*
*  q - Cx1 vector of qualities of core-periphery pairs. q[i] is the quality of the ith core-periphery pair, and
*      C is the number of core-periphery pairs. 
*
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
*   A[i][j] - Adjacency matrix. A[i][j] is the node id of the jth neighbour of node i.
*
*   c - Nx1 vector. c[i] is the index of core-periphery pair to which each node i belongs, and
*       N is the number of nodes. 
*
*   x - Nx1 vector. x[i] = true if node i is a core node and x[i] = false if it is a peripheral node.
*
*
* OUTPUT:
*
*   Q - Quality of core-periphery pairs.  
*
*   q - Cx1 vector of qualities of core-periphery pairs. q[i] is the quality of the ith core-periphery pair, and
*       C is the number of core-periphery pairs. 
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
*   A[i][j] - Adjacency matrix. A[i][j] is the node id of the jth neighbour of node i.
*
*   c - Nx1 vector. c[i] is the index of core-periphery pair to which each node i belongs, and
*       N is the number of nodes. 
*
*   x - Nx1 vector. x[i] = true if node i is a core node and x[i] = false if it is a peripheral node.
*
*   num_of_runs - number of runs of the KM algorithm.
*
*   num_of_rand_nets - number of randomised networks.
*
*
* OUTPUT
*
*   p_values - Cx1 vector. p_values[i] is the p-value of the ith core-periphery pair. 
*
*/
void estimate_statistical_significance(const vector<vector<int> >& A,
    const vector<int>& c,
    const vector<bool>& x,
    const int num_of_runs,
    const int num_of_rand_nets,
    vector<double>& p_values);
