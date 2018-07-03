import _kmalgorithm as km 

def km_config(edges, w, significance_level):
	cppairs = km.detect_config(edges.astype(int), w.astype(float), significance_level = significance_level)
	return cppairs
