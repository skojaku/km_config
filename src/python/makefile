# Makefile
.PHONY: all

all: 
	mkdir build 
	cd build && cmake ../ && make && mv *.so ../
	
.PHONY: clean
clean:
	rm -rf build *.so
