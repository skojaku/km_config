# Makefile
.PHONY: all

#CC := gcc
CC := g++

#CFLAGS := -O3 -std=c++11 # use this option in case openmp does not work 
CFLAGS := -O3 -std=c++11 -fopenmp

all: km_config

km_config:include/* src/*
	sudo rm -rf build km_config.egg* && sudo python3 setup.py build install

.PHONY: clean
clean:
	$(RM) km_config
