# Makefile
.PHONY: all

cpp:
	make -C src/cpp

matlab:
	make -C src/matlab

python:
	make -C src/python

.PHONY: clean
clean:
	make -C src/cpp clean
	make -C src/matlab clean
	make -C src/python clean
