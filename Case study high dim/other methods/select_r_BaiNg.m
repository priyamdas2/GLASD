function r_hat = select_r_BaiNg(Y, rmax)
    [n,p] = size(Y);
    Y = Y - mean(Y,1);
    
    Sigma = (Y' * Y)/n;
    Sigma = (Sigma + Sigma')/2;
    
    [V,D] = eig(Sigma);
    [lambda,idx] = sort(diag(D),'descend');
    
    IC = zeros(rmax+1,1);
    
    for r = 0:rmax
        if r == 0
            Sigma_r = zeros(p);
        else
            Vr = V(:,idx(1:r));
            Lr = diag(lambda(1:r));
            Sigma_r = Vr * Lr * Vr';
        end
        
        Res = Sigma - Sigma_r;
        IC(r+1) = log(mean(Res(:).^2)) + r*(p+n)/(p*n)*log(min(p,n));
    end
    
    [~,ind] = min(IC);
    r_hat = ind - 1;
end
