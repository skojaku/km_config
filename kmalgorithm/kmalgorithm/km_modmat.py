import _kmalgorithm as km 
import numpy as np
import networkx as nx

def km_modmat(G, significance_level = 0.05, num_of_runs = 10, num_of_rand_nets = 500):
	if isinstance(G,nx.classes.graph.Graph) == False:
		print("Pass Networkx.classes.graph.Graph object")
		return

	#Gi = nx.convert_node_labels_to_integers(G)
	node2id = dict(zip(G.nodes, range(len(G.nodes))))
	id2node= dict((v,k) for k,v in node2id.items())

	nx.relabel_nodes(G, node2id,False)
	edges = G.edges(data="weight")	

	node_pairs = np.array([ [edge[0], edge[1]] for edge in edges ]).astype(int)
	w = np.array([ edge[2] for edge in edges ]).astype(float)
	
	if all(np.isnan(w)):
		nx.set_edge_attributes(G, values =1, name='weight')
		w[:] = 1.0

	cppairs = km.detect_modmat(edges=node_pairs, ws=w, significance_level = significance_level, num_of_runs = num_of_runs, num_of_rand_nets = num_of_rand_nets)
	
	nx.relabel_nodes(G,id2node,False)
	return {'pair_id': dict(zip( range(len(cppairs[0])), cppairs[0].astype(int))),\
		'core_node': dict(zip( range(len(cppairs[0])), cppairs[1].astype(bool))),\
		'cp_pair_significance':cppairs[2]}
	#return {'pair_id':, 'is_core_node':cppairs[1], 'cp_pair_significance':cppairs[2]}	
	cppairs = km.detect_modmat(edges.astype(int), w.astype(float), significance_level = significance_level)
	return cppairs
