function [y,W] = my_Frost(s_in, J, K, mu)
%UNTITLED2 Summary of this function goes here
%   J - filter order
%   K - Number of microphones
%   s_in - input signal s_in(:,K)
    X = zeros(K*J,1);
    y = zeros(K*J,1);
    C = zeros(K*J,J);
    W = zeros(K*J,1);
    I = diag(ones(1,K*J));
    P = zeros(K*J,K*J);
    F = zeros(K*J,1);                                    
    for j = 1:J
        C(:,j) = [zeros(1,(j-1)*K) ones(1,K) zeros(1,J*K-j*K)].';
    end
    FF = fir1(J-1,0.99,'low',chebwin(J,30)).';
    % figure(111);
    % freqz(bhi,1)
    F = C/(C.'*C)*FF;                           %!!!C* inv(C.'*C)*FF;
    P = I - C/(C.'*C)*C.';
    W = F;

    for iter = 1:length(s_in)
        X(K*J-K+1:K*J) = [];
        X = [s_in(iter, :).'; X];
        y(iter) = W.'*X;
        W = P*(W - mu*y(iter)*X) + F;
%         constr(iter) = sum(ceil(C.'*W*10000) ~= ceil(FF*10000));
    end
end

