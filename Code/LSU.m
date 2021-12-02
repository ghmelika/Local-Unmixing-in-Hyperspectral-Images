function [alpha,W_final,F,P_rdn,si2cap] = LU(r,M,W,K,p,Vi,e,nbg,cat)
% inputs: 
%        r  is the data set (L x N matrix)--> L bands, N pixles
%        M  is the matrix including endmember vectors in each band (bxp)...
%        ... (output of EIA_FIPPI.m named E)
%        W  is the observation weight matrix(weighted bands) (bxb) ...
%        ... (output of weight_matrix_Crr.m)
%        K  is the final number of categories ...
%        ... (output of homogeneous_observation_categories.m)
%        p  is positive integer number of endmembers in the scene ...
%        ... (output of EIA_HFC.m)
%        Vi is sum_matrix of cofactor matrix (bxb)...
%        ... (output of weight_matrix_Xrr.m)
%        e  is a threshold as the variance component convergence criterion 
%        nbg is a vector containing number of bands in each category (1xk)
%        cat is a matrix containing band numbers in each category (?xK)

% outputs:
%        alpha  is frequency of endmembers in each pixel
%        F      is statistic index for evaluating convergence of variance component
%        si2cap is variance component
%        P_rdn  is Redundancy Numbers 

%%
%        each pixel is a linear mixture of p endmembers
%        signatures R = M x s, where s = gamma x alfa
%        gamma is a illumination perturbation factor and
%        alfa are the abundance fractions of each endmember.
%        for a given R, we need to decide the M and s

%% Parameters
if (nargin < 2)
    error('Insufficient parameters');
end
if (nargin < 3)
    e = 10^(-6);
end
%% LSU
% r: is LxN matrix with the hyperspectral data set 
[b N] = size(r);
I = eye(b,b); % b: number of bands
alpha = []; % (p,N,K) : for each category
v = [];
sigma_cap2 = [];
Crr = zeros(b,b); % equation 5 in "paper"
SUM = 0; 
while F > e
    
    
    for i = 1:K
        % endmembers abundance matrix for each category
        alpha(:,:,i) = inv((M(1:nbg(1,i),1:p))'*W(1:nbg(i),1:nbg(i))*...
            M(1:nbg(1,i),1:p))*(M(1:nbg(1,i),1:p))'*W(1:nbg(i),1:nbg(i))*r(cat(:,i),N);
    
        % residual for each category
        v(cat(:,i),N) = (M(1:nbg(1,i),1:p))*alpha(:,:,i)-r(cat(:,i),N); 
    
        % variance of each category after LS adjustment (bxb)
        sigma_cap2(1:nbg(1,i),N,i) = ((v(cat(:,i),N))'*W(1:nbg(1,i),1:nbg(1,i))*(v(cat(:,i),N)))/...
            (nbg(1,i)-trace(inv(M'*W*M)*M(cat(:,i),1:p)'*W(1:nbg(1,i),1:nbg(1,i))*M(cat(:,i),1:p)));
        %sigma(i) = trace(sigma_cap2());
    
        sum1 = sigma_cap2(1:nbg(1,i),N,i)*Vi(:,:,i);   
        Crr = Crr + sum1;
        SUM1 = (sigma_cap2(:,:,i)-1)^2;
        SUM = SUM + SUM1; 
        
    end
    
    F = sqrt(SUM/(K+1));
    %w = diag(sigma);
end

W = inv(Crr);
si2cap = I-M*(inv(M'*W*M))*M'*W; 
% P_rdn = (v'*W*v)/df;



return