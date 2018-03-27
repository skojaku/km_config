import csv
import numpy as np
import km_config as kmconfig


linkfilename='example_edge_list.txt'
edges = np.genfromtxt(linkfilename, delimiter=' ', skip_header = 0)
cppairs = kmconfig.detect(edges.astype(float), significance_level = 0.05)
print(cppairs)
