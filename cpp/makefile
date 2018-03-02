# Makefile
.PHONY: all

CC := g++

#CFLAGS := -O3 -std=c++11 # use this option in case openmp does not work 
CFLAGS := -O3 -std=c++11 -fopenmp

all: km_config

km_config:src/km_config.cpp src/km_config.h src/km_config_cl.cpp src/config_model.h
	$(CC) $(CFLAGS) -o km_config src/km_config_cl.cpp 

.PHONY: clean
clean:
	$(RM) km_config
