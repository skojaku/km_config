#include <iostream>
#include <vector>

using namespace std;
typedef pair<int, double> Edge; // Edge 

class Graph{
public:
    vector<vector<Edge>> _edges;
    
    // constracter
    Graph(int num_nodes);
	
    // Getters
    int get_num_nodes() const;
    int get_num_edges() const;
    int degree(int nid) const;
    double wdegree(int nid) const;
    void get_weight(int nid, int j, int& nei, double& w) const;
    void print() const;
    vector<vector<double>> to_matrix() const;
    
    void addEdge(int u, int v, double w);
};
    
Graph::Graph(int num_nodes){
	vector<vector<Edge>> tmp(num_nodes, vector<Edge>(0));	
	_edges = tmp;
};

// Getter -----------
int Graph::get_num_nodes() const{
	return _edges.size();
}

int Graph::get_num_edges() const{
	int N = get_num_nodes();
	int M = 0;
    	for (int i = 0; i < N; i++) {
		 M+= _edges[i].size();
	}
	return M/2;
}

// get the id and weight of the jth neighbour of node nid
void Graph::get_weight(int nid, int j, int& nei, double& w) const{
	nei = _edges[nid][j].first;
	w = _edges[nid][j].second;
}

// get weighted degree
double Graph::wdegree(int nid) const{
	int sz = _edges[nid].size();
	double deg = 0;
    	for (int j = 0; j < sz; j++) {
        	deg+= _edges[nid][j].second;
	}
	return deg;
}

// get degree
int Graph::degree(int nid) const{
	return _edges[nid].size();
}

// add 
void Graph::addEdge(int u, int v, double w){
	Edge ed1 = make_pair(v, w);
	_edges[u].push_back(ed1);
	
	if(u==v) return;
	
	//Edge ed2 = make_pair(u, w );
	//_edges[v].push_back(ed2);
}

// add 
void Graph::print() const{
	int N = get_num_nodes();
	for(int i =0; i < N;i++){
		int sz = _edges[i].size();
		for(int j =0; j < sz;j++){
			cout<<i<<" "<<_edges[i][j].first<<" "<<_edges[i][j].second<<endl;
		}
	
	}
}

// add 
vector<vector<double>> Graph::to_matrix() const{
	int N = get_num_nodes();
	vector<vector<double>> M(N, vector<double>(N, 0));
		
	for(int i =0; i < N;i++){
		int sz = _edges[i].size();
		for(int j =0; j < sz;j++){	
			M[i][_edges[i][j].first] = _edges[i][j].second;
		}
	}
	return M;
}
