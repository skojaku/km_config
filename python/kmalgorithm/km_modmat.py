import _kmalgorithm as km 

def km_modmat(edges, w, significance_level):
	cppairs = km.detect_modmat(edges.astype(int), w.astype(float), significance_level = significance_level)
	return cppairs
