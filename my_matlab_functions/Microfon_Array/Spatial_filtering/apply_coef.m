function [y] = apply_coef(s_in, W, L, K)

    X = zeros(K*L,1);
    y = zeros(K*L,1);
 
    for iter = 1:length(s_in)
        X(K*L-K+1:K*L) = [];
        X = [s_in(iter, :).'; X];
        y(iter) = W.'*X;

    end
end

