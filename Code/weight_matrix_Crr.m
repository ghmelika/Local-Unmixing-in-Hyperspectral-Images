function [W,Vi] = weight_matrix(r,K,nbg,cat)
% inputs: 
%        r   is the data set (L x N matrix)--> L bands, N pixles
%        K   is the final number of categories...
%        ...(output of homogeneous_observation_categories.m)
%        nbg is a vector containing number of bands in each category (1xk)
%        cat is a matrix containing band numbers in each category (?xk)

% outputs:
%         W  is the weight observation matrix(in this case : weighted band)
%         Vi is submatrix of cofactor matrix for each category
%%
% r: is LxN matrix with the hyperspectral data set 
[L N] = size(r); 
%% Parameters
if (nargin < 2)
    error('Insufficient parameters');
end
%%

ri = []; % observation vector for each category
sigmai2 = zeros(K,1); % variance component for each category
Sum = 0;
for i = 1:K
    ri(1+Sum:nbg(1:i),1:N) = r(cat(sort(find(cat(:,i) ~= 0))));
    Sum = Sum + nbg(1,i);
    sigmai2(i,1) = diag(cov(ri(1+Sum:nbg(1:i),1:N))); 
end
 # Vi ???? ===> Crr ???? ===> W
 
 %%
[r k] = size(nbg); 
Crr = zeros(b,b); % equation 5 in "paper"
for i = 1:k
    sum2 = sigmai2(i,1)*Vi(i);   
    Crr = Crr + sum2;
end

W = inv(Crr);

return