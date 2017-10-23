%
% MATLAB client for the KM-config algorithm
%
%
% An algorithm for finding multiple core-periphery pairs in networks
% "???"
% Sadamori Kojaku and Naoki Masuda
% Preprint arXiv:???
% 
%
% Please do not distribute without contacting the authors.
%
%
% AUTHOR - Sadamori Kojaku
%
%
% DATE - 11 Oct 2017
function [cids, x, Q, q, p_vals] = km_config(A, varargin)

num_of_runs = 10;
num_of_rand_nets = 500;
alpha = 1; 
if nargin >= 2; num_of_runs = varargin{1}; end
if nargin >= 3; alpha = varargin{2}; end
if nargin == 4; num_of_rand_nets = varargin{3}; end

if alpha > 1-1e-4; num_of_rand_nets = 0; end

[node_ids1, node_ids2] = find(triu(A, 1));
N = size(A, 1);
[cids, x, Q, q, p_vals] = km_config_mex([node_ids1, node_ids2], N, length(node_ids1), num_of_runs, num_of_rand_nets);
	
if alpha > 1-1e-4; return; end % 
	
K = max(cids);
significant = p_vals <= ( 1 - (1 - alpha) ^(1 / size(K, 2)) ); 
x(~significant(cids)) = NaN;
cids(~significant(cids)) = NaN;
p_vals = p_vals(significant);
q = q(significant);
Q = sum(q);

end
