function [nbg,cat,K] = category(w,k)
% [nVar_norm,Co] = Cov_Mat(w,k)

% inputs: 
%        w is the noise estimates for every pixel (LxN)
%        k is the number of equal intervals as the number of categories

% outputs:
%         nbg is a vector containing number of bands in each category (1xk)
%         cat is a matrix containing band numbers in each category (?xk)
%         K  is final number of categories

%         Co is the noise covariance matrix (LxL)

%% homogeneous observation categories determination 
%% Parameters
if (nargin < 1)
    error('Insufficient parameters');
end
if (nargin < 2 )
    k = 10;
end

%% noise variance estimation

% w is the output of 'estNoise.m'
% w is the noise estimates for every pixel (LxN)

% %% Covariance calculation
% [L, N] = size(w);
% % Remove mean from data
% u = mean(w.').';
% for i = 1:N
%     w(:,i) = w(:,i) - u;
% end
% 
% Co = (w*w.')/(N-1); % noise covariance matrix (LxL)
% 
% %% noise variance calculation
% nVar = diag(Co); % noise variance for each band (Lx1)
nVar = noisevar;
%% Normalization
minVal = min(nVar(:));
maxVal = max(nVar(:));

% normalized noise variance for each band (Lx1)
nVar_norm = nVar - minVal;
if (maxVal == minVal)
    normalizeData = zeros(size(nVar));
else
    nVar_norm = nVar_norm ./ (maxVal-minVal);
end
%% homogeneous observation categories
% for n categories in a baseline between 0 to 1:
% separating bands according to their noise variance into...
% homogeneous observation categories
cat = [];
for i = 1:k
    cat(1:length(find(nVar_norm <= i/k & nVar_norm >= ((i-1)/k))),i) = find(nVar_norm <= i/k & nVar_norm >= ((i-1)/k));
end

%% conditions for eliminating some obs. groups and merging with others 

% threshold for number of bands in each category: 
% number of bands in each group
nbg = zeros(1,k);
for i = 1:k
    if max(cat(:,i)) ~= 0
    nbg(1,i) = max(size(find(cat(:,i) ~= 0)));
    end
end

% eliminate the groups including no band or ...
% groups with bands less than threshold
T = 0.02*L; % threshold : 2% of all bands

% for first and last category:
% first:
if nbg(1,1) < T
    cat(nbg(1,2) + 1 : nbg(1,2) + nbg(1,1),1) = cat(1:nbg(1,1),1);
    cat(:,1) = [];
    k = k-1;
else
    cat(1,1) = cat(1,1);
end

for i=1:k
    nbg(1,i) = max(length(find(cat(:,i) ~= 0)));
end

for i = 2:k-1
    if nbg(1,i) < T 
        if nbg(1,i-1) < nbg(1,i+1)
            cat(nbg(1,i-1) + 1 : nbg(1,i-1) + nbg(1,i),i) = cat(1:nbg(1,i),i);
            cat(:,i) = [];
        else
            cat(nbg(1,i+1) + 1 : nbg(1,i+1) + nbg(1,i),i) = cat(1:nbg(1,i),i);
            cat(:,i) = [];
        end
    else 
        cat(:,i) = cat(:,i);
    end
end
[r c] = size(cat); % new number of categories

for i=1:c
    nbg(1,i) = max(size(find(cat(:,i) ~= 0)));
end

% last:
if nbg(1,c) < T
    cat(nbg(1,c-1) + 1 : nbg(1,c-1) + nbg(1,c),c) = cat(1:nbg(1,c),c);
    cat(:,c) = [];
else
    cat(1,c) = cat(1,c);
end
[row col] = size(cat); % final number of categories
K = col; % final number of categories
for i=1:col
    nbg(1,i) = max(size(find(cat(:,i) ~= 0)));
end

return