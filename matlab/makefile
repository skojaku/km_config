# Makefile
.PHONY: all

MEXFLAGS := -fopenmp

MEXOFLAGS := -O3

all: km_config_mex.mexa64

km_config_mex.mexa64: km_config_mex.cpp ../cpp/km_config.h ../cpp/km_config.cpp
	mex CXXFLAGS='$$CXXFLAGS $(MEXFLAGS)' LDFLAGS='$$LDFLAGS $(MEXFLAGS)' CXXOPTIMFLAGS='$(MEXOFLAGS) -DNDEBUG' LDOPTIMFLAGS='$(MEXOFLAGS)' km_config_mex.cpp

.PHONY: clean
clean:
	$(RM) *.mexa64 
