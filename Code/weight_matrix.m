function [W,Vi] = weight_matrix(r,K,nbg,cat,E,e)
% inputs: 
%        r   is the data set (L x N matrix)--> L bands, N pixles
%        K   is the final number of categories...
%        ...(output of homogeneous_observation_categories.m)
%        nbg is a vector containing number of bands in each category (1xK)
%        cat is a matrix containing band numbers in each category (?xK)
%        E   is set of induced endmembers [L x p] : output of EIA_ATGP.m
%        e   is a threshold for F
% outputs:
%         W  is the weight observation matrix(in this case : weighted band)
%         Vi is submatrix of cofactor matrix for each category
%%
% r: is LxN matrix with the hyperspectral data set 
[L, N] = size(r); 
%% Parameters
if (nargin < 5)
    error('Insufficient parameters');
end

if (nargin < 6)
    e = 10^-6;
end
%% iniatl observation covariance matrix in one segment (LxL)

M = mean(r.').';
for i = 1:N
    r(:,i) = r(:,i) - M;
end
Crr_in = (r*r.')/(N-1); 

%% weight matrix : inverse of covariance matrix 

W_in = inv(Crr_in); % (L x L)

%% M : matrix including endmember vectors in all "L" bands which is...
%      output of EIA_ATGP. its the initial value of endmembers (not weighted ATGP)

M = E; % (L x p) : p is number of endmembers

%% comuting P_ort which inverts observation space to residual space

I = eye(L,L);
P_ort = I - M*inv(M'*W_in*M)*M'*W_in; % (L x L)

%% v_in : computing residual of N pixels of the image...
%         (in one segment)in all "L" bands

v_in = P_ort*r; % (L x N)

%% estimation of variance component : sigma_i_cap

% 1. v_i : computing residual vector of each category (bi x N)
%        (N pixels in specified bands of each category)

v_i = [];

for i = 1:K
    v_i(1:nbg(1,i),1:N,i) = v_in(cat(find(cat(:,i)~=0)),1:N); 
end

% 2. W_i : computing weight matrix of each category (bi x bi)

W_i = [];
for i = 1:K
    W_i(1:nbg(1,i),1:nbg(1,i),i) = W_in(cat(find(cat(:,i)~=0)),cat(find(cat(:,i)~=0))); 
end

% 3. M_i : computing matrix including EMs vectors in each category (bi x p)

[~, p] = size(M); % (L x p)

for i = 1:K
    M_i(1:nbg(1,i),1:p,i) = M(cat(:,i),1:p); 
end

% (3.5). change v_i and W_i formats for including all N pixels of segment

% reshape of v_i : (bi.N x 1)
for i = 1:K
    vi(1:N*nbg(1,i),i) = reshape(v_i(1:nbg(1,i),1:N,i),[N*nbg(1,i),1]); 
end

% repeat of W_i for all N pixels : (bi.N x bi.N)
for i = 1:K
    for j = 0:N-1
        Wi(1+(j*nbg(1,i)):(j+1)*nbg(1,i) , 1+(j*nbg(1,i)):(j+1)*nbg(1,i) , i) = W_i(1:nbg(1,i),1:nbg(1,i),i);
    end
end
%% determining cofactor matrices : Vi    
Vi = zeros(L,L);
i = 1; 
Vi(1:nbg(1,i),1:nbg(1,i),i) = Crr_i(:,:,i);
s = nbg(1,1);
for i = 1:K
    Crr_i(1:nbg(1,i),1:nbg(1,i),i) = inv(W_i(1:nbg(1,i),1:nbg(1,i),i));
end

for i = 2:K
    Vi(s+1:s+nbg(1,i),s+1:s+nbg(1,i),i) = Crr_i(:,:,i);
    s = s+nbg(1,i);
end
%% %%%%%%%%%%%%%%%%%%%%%%%

% 4. sigma_i_cap : variance component of each category

sigma2_i_cap = zeros(1,K);

for i =1:K
    sigma2_i_cap(1,i) = ((vi(1:N*nbg(1,i),i))'*Wi(1+(j*nbg(1,i)):(j+1)*nbg(1,i) , 1+(j*nbg(1,i)):(j+1)*nbg(1,i) , i)*vi(1:N*nbg(1,i),i))/(nbg(1,i)-trace((inv(M'*W*M))*(M_i(1:nbg(1,i),1:p,i))'*W_i(1:nbg(1,i),1:nbg(1,i),i)*M_i(1:nbg(1,i),1:p,i)));
end

%% update time %%

% for i = 0 :

df = trace(P_ort);

% reshape of residual matrix(including all bands) to residual vector(NL x 1)
v = reshape(v_in,[N*L,1]); 

% repeat of W_i for all N pixels : (NL x NL)
for i = 0:N-1
    W(1+(i*L):(i+1)*L , 1+(i*L):(i+1)*L) = W_in;
end

% determining sigma2_i_cap
sigma0_cap2 = (v'*W*v)/df;
sums = 0;
for i = 1:K
    sums = (sigma2_i_cap(1,i)-1)^2 + sums;
end
sums = sums + sigma0_cap2;
F = sqrt(sums/K+1);

%% after this step we need to check if F>e or not
%  if F>e we have to repeat LSU but with some changes:

% computing Crr and W(weight matrix) after update
sums = 0;
for i = 1:K
    sums = sigma2_i_cap(1,i)*Vi(:,:,i);
end
sums = sums + sigma0_cap2;
Crr_updt = sigma2_i_cap(1,i);
    
W_updt = inv(Crr_updt);

% now : for update we must use following W_in in all above formulas 
W_in = W_updt;
% so all the other variables will be updated and the process is repeated
end

%%  final W:
% final W is the one which derived from the step that F<e
W = W_in;