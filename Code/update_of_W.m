% inputs: r , E , 

F = 1;
if F > e
        
    sums = 0;
    for i = 1:K
        sums = sigma2_i_cap(1,i)*Vi(:,:,i);
    end
    sums = sums + sigma0_cap2;
    Crr_updt = sigma2_i_cap(1,i);
    
    W_updt_updt = inv(Crr_updt);
    
    M = E; % (L x p) : p is number of endmembers
    
    I = eye(L,L);
    P_ort_updt = I - M*inv(M'*W_updt*M)*M'*W_updt; % (L x L)
    v_updt = P_ort_updt*r; % (L x N)
    
    v_i_updt = [];
    
    for i = 1:K
    v_i_updt(1:nbg(1,i),1:N,i) = v_updt(cat(find(cat(:,i)~=0)),1:N); 
    end
    
    % 2. W_updt_i_updt : computing W_updteight matrix of each category (bi x bi)
    
    W_updt_i_updt = [];
    
    for i = 1:K
        W_updt_i_updt(1:nbg(1,i),1:nbg(1,i),i) = W_updt_updt(cat(find(cat(:,i)~=0)),cat(find(cat(:,i)~=0))); 
    end
    
    
    %% constant
    % 3. M_i : computing matrix including EMs vectors in each category (bi x p)
    
    [~, p] = size(M); % (L x p)
    
    for i = 1:K
        M_i(1:nbg(1,i),1:p,i) = M(cat(:,i),1:p); 
    end
    
    
    % (3.5). change v_i and W_updt_i formats for including all N pixels of segment
    
    % reshape of v_i_updt : (bi.N x 1)
    
    for i = 1:K
        vi_updt(1:N*nbg(1,i),i) = reshape(v_i_updt(1:nbg(1,i),1:N,i),[N*nbg(1,i),1]); 
    end
    
    % repeat of W_updt_i for all N pixels : (bi.N x bi.N)
    
    for i = 1:K
        for j = 0:N-1
            W_updti_updt(1+(j*nbg(1,i)):(j+1)*nbg(1,i) , 1+(j*nbg(1,i)):(j+1)*nbg(1,i) , i) = W_updt_i_updt(1:nbg(1,i),1:nbg(1,i),i);
        end
    end
    
    % 4. sigma_i_cap : variance component of each category

    sigma2_i_cap = zeros(1,K);

    for i =1:K
        sigma2_i_cap(1,i) = ((vi_updt(1:N*nbg(1,i),i))'*W_updti_updt(1+(j*nbg(1,i)):(j+1)*nbg(1,i) , 1+(j*nbg(1,i)):(j+1)*nbg(1,i) , i)*vi_updt(1:N*nbg(1,i),i))/(nbg(1,i)-trace((inv(M'*W_updt*M))*(M_i(1:nbg(1,i),1:p,i))'*W_updt_i_updt(1:nbg(1,i),1:nbg(1,i),i)*M_i(1:nbg(1,i),1:p,i)));
    end
    
    % for i = 0 :

    df = trace(P_ort_updt);

    % reshape of residual matrix(including all bands) to residual vector(NL x 1)
    v = reshape(v_updt,[N*L,1]); 

    % repeat of W_updt_i for all N pixels : (NL x NL)
    for i = 0:N-1
        W_updt(1+(i*L):(i+1)*L , 1+(i*L):(i+1)*L) = W_updt;
    end
    
    sigma0_cap2 = (v'*W_updt*v)/df;
    sums = 0;
    for i = 1:K
        sums = (sigma2_i_cap(1,i)-1)^2 + sums;
    end
    sums = sums + sigma0_cap2;
    F = sqrt(sums/K+1);
end