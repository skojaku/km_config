# Makefile
.PHONY: all

#CC := gcc
CC := g++

#CFLAGS := -O3 -std=c++11 # use this option in case openmp does not work 
CFLAGS := -O3 -std=c++11 -fopenmp

all: km_config

km_config:src/* kmalgorithm/*
	sudo rm -rf build kmalgorithm.egg* && sudo python3 setup.py build install

.PHONY: clean
clean:
	$(RM) km_config
