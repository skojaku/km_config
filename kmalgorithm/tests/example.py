import csv
import numpy as np
import pandas as pd
import networkx as nx

#import km_config as kmconfig
from kmalgorithm import km_config, km_modmat

#linkfilename='example_edge_list.txt'
#df = pd.read_csv(linkfilename, sep=' ')
G=nx.karate_club_graph()
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
print(w, node_pairs)
cppairs = km_config(G, significance_level = 0.05)

print(cppairs)
