# Makefile
.PHONY: all

MEXFLAGS := -fopenmp
MEXOFLAGS := -O3

#MEXCOMPILER := gcc
MEXCOMPILER := g++-4.7

all: km_config_mex.mexa64

km_config_mex.mexa64: km_config_mex.cpp ../lib/km_config.h ../lib/km_config.cpp
	mex CXXFLAGS='$$CXXFLAGS $(MEXFLAGS)' LDFLAGS='$$LDFLAGS $(MEXFLAGS)' CXXOPTIMFLAGS='$(MEXOFLAGS) -DNDEBUG' LDOPTIMFLAGS='$(MEXOFLAGS)' GCC='$(MEXCOMPILER)' km_config_mex.cpp

.PHONY: clean
clean:
	$(RM) *.mexa64 
