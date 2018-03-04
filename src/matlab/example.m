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
