import csv
import numpy as np
import pandas as pd
#import km_config as kmconfig
from kmalgorithm import km_config, km_modmat

linkfilename='example_edge_list.txt'
df = pd.read_csv(linkfilename, sep=' ')

edges = np.array(df[["source", "target"]].values-1)
ws = np.array(df[["value"]].values)

cppairs = km_modmat(edges.astype(int), ws, significance_level = 1)

print(cppairs)
