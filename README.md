# km_config
C++, Matlab and Python code for the KM--config algorithm.

Please cite  
  Kojaku, S. and Masuda, N. "Core-periphery structure requires somthing else in the network". Preprint arXiv:1710.07076 (2017).


## Files
Directory src/lib contains the main code written in C++ 
  * src/lib/km_config.h is the header file.
  * src/lib/km_config.cpp is the implementation file of the header file, km_config.h.

Directory src/cpp contains C++ code for command-line cliend 
  * src/cpp/km_config_cl.cpp is the code for the command-line client.
  * src/cpp/example.m is a usage example.
  * src/cpp/example_edge_list.txt is an edge list of a network consisting of two idealised core-periphery pairs.
  * src/cpp/makefile is the makefile for the command-line client. 

Directory src/matlab contains Matlab wrapper 
  * src/matlab/km_config_mex.cpp is the MATLAB wrapper for the C++ code (km_config.h and km_config.cpp)
  * src/matlab/km_config.m is the code for the MATLAB client.
  * src/matlab/example.m is a usage example.
  * src/matlab/makefile is the makefile for the MATLAB code. 

Directory src/python contains Python wrapper 
  * src/python/km_config.cpp is the Python wrapper for the C++ code (km_config.h and km_config.cpp)
  * src/python/example.py is a usage example.
  * src/python/makefile is the makefile for the Python code. 
  * src/python/CMakeLists.txt is the cmake file for the Python code. 
  * src/python/example_edge_list.txt is an edge list of a network consisting of two idealised core-periphery pairs.


## C++ 

### Compile:

  To compile, type
        
  ```bash
  make cpp
  ```
       
  This will produce an executable file "km_config" in the src/cpp directory.
 
 
### Usage:

  ``` bash
  km_config [input-file] [output-file] [options]
  ```

  km_config seeks multiple core-periphery pairs in the network given by [input-file] and saves the detected core-periphery pairs in [output-file].
  
**[input-file]**  
 * The file should contain a list of edges (space-separated).  
 * The first and second columns represent the IDs of the two nodes forming an edge.
 * The third column represents the weight of the edge between the two nodes. 
 * If the third column is not provided, then the weight is set to 1.  
 * The node's ID is assumed to start from 1. You can change this by [options], e.g., -i 0.
  
**[output_file]**  
 * The first column represents the node's ID.
 * The second column represents the index of the core-periphery pair to which each node belongs.
 * The third column indicates whether each node is a core node (= 1) or a peripheral node (= 0).
 * The fourth column indicates whether each node belongs to a significant core-periphery pair (= 1) or not (= 0).
  
  
**[options]**  
* -r R  
 Run the KM-config algorithm R times. (Default: 10)  
* -a ALPHA  
  Set the significance level before the Šidák correction to ALPHA. (Default: 1). If this option is not set, the statistical test is not carried out.
* -l N  
Set the number of randomised networks to N. (Default: 500)
* -i I
The node's ID starts from I. (Default: 1)
* -d "D"  
Change the delimiter for [input-file] and [output-file] to D. (Default: space)  


### Example (src/cpp/example.sh):
    
```bash
#To find core-periphery pairs, type
./km_config example_edge_list.txt result.txt
    
#To find significant core-periphery pairs at a significance level of 0.05, type
./km_config example_edge_list.txt result.txt -a 0.05 
```


## Matlab 
      
### Compile:

  To compile, type
        
  ```bash
  make matlab
  ```
    
  This will produce a mex file "km_config_mex.mexa64" in the src/matlab/ directory. 
  Copy "matlab/km_config_mex.mexa64" and "src/matlab/km_config.m" to your working directory.


  You may have the following message:
  
  ```bash
  Warning: You are using gcc version '5.4.0'. The version of gcc is not supported. 
  The version currently supported with MEX is '4.9.x'. 
  ```
  
  This means you are required to change the version of g++ compiler. 
  To remedy this, modify a line in ''./src/matlab/makefile'' as follows: 
  
  ```bash
  MEXCOMPILER := g++-(the version compatible with your mex compiler, e.g., g++-4.9) 
  ```
 
### Uage:

    [c, x, Q, q, p_vals] = km_config(A, num_of_runs, alpha, num_of_rand_nets);
 
 
  **INPUT:** 
 
  * `A` - N times N adjacency matrix, where N is the number of nodes. A(i, j) = w if nodes and j are adjacent, where w is the edge weight. Otherwise A(i, j) = 0. Matrix A should be symmetric, i.e., A(i, j) = A(j, i).
      
  * (optional) `num_of_runs` - Number of runs. (Default: 10) 
      
  * (optional) `alpha` - Statistical significance level before the Šidák correction. (Default: 1) If alpha is not set, the statistical test is not carried out. 
      
  * (optional) `num_of_rand_nets` - Number of randomised networks. (Default: 500) 


  **OUTPUT:**

  * `c` - N-dimensional column vector. c(i) is the index of the core-periphery pair to which node i belongs.
          If node i belongs to an insignificant core-periphery pair, then c(i) = NaN.
      
  * `x` - N-dimensional column vector. If node i is a core node, x(i) = 1. If node i is a periphery node, x(i) = 0.
          If node i belongs to an insignificant core-periphery pair, then x(i) = NaN.
      
  * `Q` - The quality value of the detected core-periphery pair.
      
  * `q` - C-dimensional column vector. q(i) is the contribution of the i-th core-periphery pair to Q.
          C is the number of core-periphery pairs.
      
  * `p_vals` - C-dimensional column vector. p_vals(i) is the statistical significance of the i-th core-periphery pair. 

  
### Example (src/matlab/example.m)
    
```matlab
% Construct a network with two core-periphery pairs 
A = [
	ones(5),ones(5);
	ones(5),zeros(5);
];
A = A - diag(diag(A));
A = kron(eye(2), A);

% Display adjacency matrix
disp('---- Adjacency matrix ----') 
full(A)

% Run KM-config algorithm
[c, x, Q, q] = km_config(A) % without statistical test

[c, x, Q, q, p_values] = km_config(A, 10, 0.05) % with statistical test
```


## Python 
      
### Compile:

  To compile, type
        
  ```bash
  make python 
  ```
    
This creates a shared library ''src/python/km_config.cpython-36m-x86_64-linux-gnu.so'' callable from python. 
Copy the shared library to your working directory. 
 
### Uage:

```python
import km_config as kmconfig
cppairs = kmconfig.detect(edges, num_of_runs = 10, significance_level = 0.05, num_of_rand_nets = 500)
```
 
 
  **INPUT:** 
 
  * `edges` - Mx3 Numpy array, where M is the number of edges. The first and second columns indicate the IDs of nodes connected by an edge. The third column indicates the weight of the edge.
      
  * (optional) `num_of_runs` - Number of runs. (Default: 10) 
      
  * (optional) `significance_level` - Statistical significance level before the Šidák correction. (Default: 1) If alpha is not set, the statistical test is not carried out. 
      
  * (optional) `num_of_rand_nets` - Number of randomised networks. (Default: 500) 


  **OUTPUT:**

  * `cppairs` - List of length 4. 
    * cppairs[0] - Numpy array of length N, where N is the number of nodes. communities[0][i] indicates the index of the core-periphery pair to which node i belongs.
    * cppairs[1] - Numpy array of length N. cppairs[1][i] = 1 or = 0 indicates a core or a peripheral node, respectively.
    * cppairs[1] - Numpy array of length N. cppairs[2][i] indicates the p-value of the core-periphery pair to which node i belongs.
    * cppairs[2] - Numpy array of length N. cppairs[3][i] = True or False indicates that node i belongs to the significant or insignificant cppairs, respectively.

  
### Example (src/python/example.m)
    
```python
import csv
import numpy as np
import km_config as kmconfig

linkfilename='example_edge_list.txt'
edges = np.genfromtxt(linkfilename, delimiter=' ', skip_header = 0)
cppairs = kmconfig.detect(edges, significance_level = 0.05)

print(cppairs)
```


### REQUIREMENT: 
      
  MATLAB 2012 or later.
 * Python3.4 or later
 * Cmake2.8 or later

---
Last updated: 17 October 2017
