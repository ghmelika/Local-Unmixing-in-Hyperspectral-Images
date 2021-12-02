function [alpha] = fraction_map(r,W,E)

% inputs:
%        r is the data set (L x N matrix)--> L bands, N pixles
%        W is the weight matrix (L x L matrix)(output of weight_matrix.m)
%        E is the endmember matrix (L x P matrix)(output of weighted EIA_ATGP.m)

% output:
%        alpha is the portion of each EM in all pixels (P x N)

M = E;


[L, N] = size(r); 
alpha = zeros(P,N);
alpha = (inv(M'*W*M))*M'*W*r;


end