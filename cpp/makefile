# Makefile
.PHONY: all

CC := g++

#CFLAGS := -O3 -std=c++11 # use this option in case openmp is not working 
CFLAGS := -O3 -std=c++11 -fopenmp

all: km_config

km_config:km_config.cpp km_config.h km_config_cl.cpp
	$(CC) $(CFLAGS) -o km_config km_config_cl.cpp 

.PHONY: clean
clean:
	$(RM) km_config
